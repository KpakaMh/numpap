import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy import stats
from dataclasses import dataclass
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_curve, auc

from core.cancer_catalog import CANCER_CATALOG, CancerProfile


@dataclass
class PatientState:
    time: float = 0.0
    conc_central: float = 0.0
    conc_periph: float = 0.0
    tumor_sensitive: float = 0.39    # xS : cellules proliférantes sensibles (g)
    tumor_resistant: float = 0.01    # xR : sous-clone résistant (g)
    tumor_transit1: float = 0.0      # x2 : cellules endommagées en transit
    tumor_transit2: float = 0.0      # x3
    tumor_transit3: float = 0.0      # x4
    immune_effector: float = 0.10    # E(t) : activité lymphocytaire cytotoxique (CD8+)
    tumor_weight: float = 0.40
    angiogenesis_index: float = 1.0
    invasion_score: float = 1.0
    progression_stage: int = 2
    stage_label: str = "Stade 2 : Expansion clonale"


STAGE_LABELS = {
    1: "Stade 1 : Initiation tumorale",
    2: "Stade 2 : Expansion clonale",
    3: "Stade 3 : Dysplasie & Carcinome in situ",
    4: "Stade 4 : Invasion locale (EMT & MMP)",
    5: "Stade 5 : Intravasation vasculaire",
    6: "Stade 6 : Dissémination circulatoire (CTC)",
    7: "Stade 7 : Extravasation & Métastases"
}


class PatientPKPDEngine:
    def __init__(self, cancer_key: str = "lung_nsclc", biomarker_b0: float = 1.0,
                 pdl1_tmb: float = 1.0, initial_res_frac: float = 0.02,
                 anti_vegf: float = 0.0, immunotherapy: float = 0.0,
                 vc: float = 5.0, vp: float = 15.0, k10: float = 0.15,
                 k12: float = 0.08, k21: float = 0.04, psi: float = 20.0, k1: float = 0.025):
        self.profile: CancerProfile = CANCER_CATALOG[cancer_key]
        self.vc = vc
        self.vp = vp
        self.k10 = k10
        self.k12 = k12
        self.k21 = k21
        self.psi = psi
        self.k1 = k1

        # Paramètres biologiques du cancer
        self.lambda0 = self.profile.lambda0
        self.lambda1 = self.profile.lambda1
        self.stromal_permeability = self.profile.stromal_permeability
        self.base_k2 = self.profile.base_potency_k2
        self.emt_invasiveness = self.profile.emt_invasiveness

        # Paramètres patient & combinaisons thérapeutiques
        self.biomarker_b0 = biomarker_b0
        self.pdl1_tmb = pdl1_tmb
        self.anti_vegf = anti_vegf          # Intensité blocage anti-VEGF (0.0 à 1.0)
        self.immunotherapy = immunotherapy  # Intensité blocage PD-1/PD-L1 (0.0 à 1.0)
        self.mutation_rate = 1.2e-4         # Taux de transition phénotypique vers la résistance

        w0 = self.profile.initial_mass
        res_frac = float(np.clip(initial_res_frac, 0.0, 0.5))
        self.state = PatientState(
            tumor_sensitive=w0 * (1.0 - res_frac),
            tumor_resistant=w0 * res_frac,
            tumor_weight=w0
        )
        self._update_progression_metrics()

    @property
    def effective_k2_sensitive(self) -> float:
        biomarker_hill = (2.0 * self.biomarker_b0) / (0.5 + self.biomarker_b0)
        return self.base_k2 * self.stromal_permeability * biomarker_hill

    @property
    def effective_k2_resistant(self) -> float:
        # Le sous-clone résistant conserve seulement 5 % de sensibilité à la thérapie primaire
        return 0.05 * self.effective_k2_sensitive

    def administer_bolus(self, dose_mg: float):
        self.state.conc_central += dose_mg / self.vc

    def _ode_system(self, t, y, dose_rate):
        Cc, Cp, xS, xR, x2, x3, x4, E = y

        # 1. Pharmacocinétique bi-compartimentale
        dCc = -(self.k10 + self.k12) * Cc + self.k21 * Cp + (dose_rate / self.vc)
        dCp = self.k12 * Cc - self.k21 * Cp

        # 2. Croissance tumorale (Simeoni) modulée par l'anti-angiogénique (anti-VEGF)
        w = max(1e-6, xS + xR + x2 + x3 + x4)
        eff_lambda1 = self.lambda1 * (1.0 - 0.55 * self.anti_vegf)
        denom = (1.0 + (self.lambda0 / max(1e-4, eff_lambda1) * w) ** self.psi) ** (1.0 / self.psi)

        growth_S = (self.lambda0 * xS) / denom
        growth_R = (0.88 * self.lambda0 * xR) / denom  # Léger coût adaptatif de la résistance

        # 3. Cytotoxicité directe (thérapie ciblée/chimio) + lyse immunitaire (CD8+)
        drug_kill_S = self.effective_k2_sensitive * Cc * xS
        drug_kill_R = self.effective_k2_resistant * Cc * xR
        imm_kill_S = 0.015 * E * xS
        imm_kill_R = 0.015 * E * xR

        # 4. Émergence de résistance sous pression thérapeutique
        switching = self.mutation_rate * (1.0 + 2.0 * Cc) * xS

        dxS = growth_S - drug_kill_S - imm_kill_S - switching
        dxR = growth_R - drug_kill_R - imm_kill_R + switching
        dx2 = (drug_kill_S + drug_kill_R + imm_kill_S + imm_kill_R) - self.k1 * x2
        dx3 = self.k1 * (x2 - x3)
        dx4 = self.k1 * (x3 - x4)

        # 5. Dynamique de l'effecteur immunitaire (activé par PD-1/PD-L1 + relargage antigénique x2)
        imm_stimulation = 0.08 * self.immunotherapy * (self.pdl1_tmb / (0.8 + self.pdl1_tmb)) * (1.0 + 1.5 * x2)
        dE = imm_stimulation - 0.05 * (E - 0.10)

        return [dCc, dCp, dxS, dxR, dx2, dx3, dx4, dE]

    def _update_progression_metrics(self):
        w = self.state.tumor_weight
        base_angio = 0.5 + 1.5 * (w / 4.0)
        self.state.angiogenesis_index = float(np.clip(base_angio * (1.0 - 0.65 * self.anti_vegf), 0.2, 3.0))

        if self.profile.is_hematologic:
            raw_score = 3.5 + 1.4 * np.log1p(w * 2.5) * self.emt_invasiveness
        else:
            raw_score = 1.5 + 1.8 * np.log1p(w * 2.0) * (0.6 + 0.4 * self.emt_invasiveness)

        self.state.invasion_score = float(np.clip(raw_score, 1.0, 7.0))
        stage = int(np.clip(np.round(self.state.invasion_score), 1, 7))
        self.state.progression_stage = stage
        self.state.stage_label = STAGE_LABELS[stage]

    def step(self, dt: float = 2.0, dose_rate: float = 0.0):
        y0 = [
            self.state.conc_central, self.state.conc_periph,
            self.state.tumor_sensitive, self.state.tumor_resistant,
            self.state.tumor_transit1, self.state.tumor_transit2,
            self.state.tumor_transit3, self.state.immune_effector
        ]
        sol = solve_ivp(self._ode_system, [0, dt], y0, args=(dose_rate,), method="RK45")
        y_next = sol.y[:, -1]

        self.state.time += dt
        self.state.conc_central = max(0.0, y_next[0])
        self.state.conc_periph = max(0.0, y_next[1])
        self.state.tumor_sensitive = max(0.0, y_next[2])
        self.state.tumor_resistant = max(0.0, y_next[3])
        self.state.tumor_transit1 = max(0.0, y_next[4])
        self.state.tumor_transit2 = max(0.0, y_next[5])
        self.state.tumor_transit3 = max(0.0, y_next[6])
        self.state.immune_effector = max(0.0, y_next[7])

        self.state.tumor_weight = (
            self.state.tumor_sensitive + self.state.tumor_resistant +
            self.state.tumor_transit1 + self.state.tumor_transit2 + self.state.tumor_transit3
        )
        self._update_progression_metrics()
        return self.state


def run_virtual_cohort_trial(cancer_key: str, n_patients: int = 200, bolus_mg: float = 120.0,
                             bolus_interval_h: int = 72, anti_vegf: float = 0.3,
                             immunotherapy: float = 0.5, duration_days: int = 28,
                             seed: int = 42) -> dict:
    """
    Génère une cohorte QSP biologiquement plausible, simule la trithérapie
    et entraîne un modèle Random Forest (Validation Croisée) pour identifier
    les biomarqueurs prédictifs de la réponse RECIST 1.1.
    """
    profile = CANCER_CATALOG[cancer_key]
    rng = np.random.default_rng(seed)

    # 1. Génération multi-paramètres avec corrélations biologiques et bornes physiologiques
    z_bio = rng.normal(0.0, 1.0, size=n_patients)
    z_imm = rng.normal(0.0, 1.0, size=n_patients)
    z_pk = rng.normal(0.0, 1.0, size=n_patients)

    biomarker_b0 = np.clip(np.exp(0.0 + 0.60 * z_bio), 0.05, 6.0)
    pdl1_tmb = np.clip(np.exp(0.0 + 0.65 * z_imm), 0.05, 8.0)
    k10_vec = np.clip(0.15 * np.exp(0.25 * z_pk), 0.04, 0.50)
    vc_vec = np.clip(5.0 * np.exp(0.20 * (0.5 * z_pk + 0.866 * rng.normal(0, 1, n_patients))), 2.5, 10.0)
    lambda0_vec = np.clip(profile.lambda0 * np.exp(0.25 * rng.normal(0, 1, n_patients)), 0.002, 0.05)
    stroma_vec = np.clip(
        profile.stromal_permeability * np.exp(0.20 * (0.35 * z_bio + rng.normal(0, 0.8, n_patients))),
        0.10, 1.0
    )
    res_frac_vec = np.clip(np.exp(np.log(0.02) + 0.65 * (-0.3 * z_bio + rng.normal(0, 1, n_patients))), 0.001, 0.30)
    w0_vec = np.clip(profile.initial_mass * np.exp(0.22 * rng.normal(0, 1, n_patients)), 0.10, 2.0)

    biomarker_hill = (2.0 * biomarker_b0) / (0.5 + biomarker_b0)
    k2_sens_vec = profile.base_potency_k2 * stroma_vec * biomarker_hill
    k2_res_vec = 0.05 * k2_sens_vec

    # 2. Intégration RK4 vectorisée (8 états x N patients)
    total_hours = duration_days * 24
    dt = 1.0
    y = np.zeros((8, n_patients))
    y[2, :] = w0_vec * (1.0 - res_frac_vec)
    y[3, :] = w0_vec * res_frac_vec
    y[7, :] = 0.10

    k12, k21, lambda1, psi, k1 = 0.08, 0.04, profile.lambda1, 20.0, 0.025
    eff_lambda1 = lambda1 * (1.0 - 0.55 * anti_vegf)
    mut_rate = 1.2e-4

    def deriv(state):
        Cc, Cp, xS, xR, x2, x3, x4, E = state
        dCc = -(k10_vec + k12) * Cc + k21 * Cp
        dCp = k12 * Cc - k21 * Cp

        w = np.maximum(1e-6, xS + xR + x2 + x3 + x4)
        denom = (1.0 + (lambda0_vec / max(1e-4, eff_lambda1) * w) ** psi) ** (1.0 / psi)

        growth_S = (lambda0_vec * xS) / denom
        growth_R = (0.88 * lambda0_vec * xR) / denom

        drug_kill_S = k2_sens_vec * Cc * xS
        drug_kill_R = k2_res_vec * Cc * xR
        imm_kill_S = 0.015 * E * xS
        imm_kill_R = 0.015 * E * xR

        switching = mut_rate * (1.0 + 2.0 * Cc) * xS

        dxS = growth_S - drug_kill_S - imm_kill_S - switching
        dxR = growth_R - drug_kill_R - imm_kill_R + switching
        dx2 = (drug_kill_S + drug_kill_R + imm_kill_S + imm_kill_R) - k1 * x2
        dx3 = k1 * (x2 - x3)
        dx4 = k1 * (x3 - x4)
        dE = 0.08 * immunotherapy * (pdl1_tmb / (0.8 + pdl1_tmb)) * (1.0 + 1.5 * x2) - 0.05 * (E - 0.10)
        return np.array([dCc, dCp, dxS, dxR, dx2, dx3, dx4, dE])

    time_grid = np.arange(0, total_hours + 1, 1)
    w_trajectories = np.zeros((len(time_grid), n_patients))
    w_trajectories[0, :] = w0_vec

    for h in range(1, total_hours + 1):
        if bolus_mg > 0 and ((h - 1) % bolus_interval_h == 0):
            y[0, :] += bolus_mg / vc_vec
        k_1 = deriv(y)
        k_2 = deriv(y + 0.5 * dt * k_1)
        k_3 = deriv(y + 0.5 * dt * k_2)
        k_4 = deriv(y + dt * k_3)
        y = np.maximum(0.0, y + (dt / 6.0) * (k_1 + 2 * k_2 + 2 * k_3 + k_4))
        w_trajectories[h, :] = np.sum(y[2:7, :], axis=0)

    # 3. Stratification Clinique RECIST 1.1
    w_final = w_trajectories[-1, :]
    resistant_burden_final = y[3, :]
    pct_change = ((w_final - w0_vec) / w0_vec) * 100.0

    recist_labels = np.where(
        w_final < 0.05, "CR (Réponse Complète)",
        np.where(pct_change <= -30.0, "PR (Réponse Partielle)",
                 np.where(pct_change >= 20.0, "PD (Progression)", "SD (Maladie Stable)"))
    )
    is_responder = ((w_final < 0.05) | (pct_change <= -30.0)).astype(int)

    # 4. Machine Learning : Découverte de Biomarqueurs Prédictifs (Random Forest + CV)
    feature_names = [
        f"Biomarqueur ({profile.biomarker_name.split('/')[0].strip()})",
        "Score Immunitaire PD-L1 / TMB",
        "Perméabilité Stromale",
        "Taux Prolifération (λ0)",
        "Fraction Résistante Basale (xR0)",
        "Clairance Patient (k10)",
        "Masse Tumorale Basale (W0)"
    ]
    X = np.column_stack([biomarker_b0, pdl1_tmb, stroma_vec, lambda0_vec, res_frac_vec, k10_vec, w0_vec])

    n_pos = int(np.sum(is_responder == 1))
    n_neg = int(np.sum(is_responder == 0))
    min_class = min(n_pos, n_neg)

    if min_class >= 5:
        n_splits = min(5, min_class)
        clf = RandomForestClassifier(n_estimators=150, max_depth=5, random_state=seed)
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        y_proba = cross_val_predict(clf, X, is_responder, cv=cv, method="predict_proba")[:, 1]
        fpr, tpr, _ = roc_curve(is_responder, y_proba)
        roc_auc = float(auc(fpr, tpr))
        clf.fit(X, is_responder)
        importances = clf.feature_importances_
        ml_status = f"Validation Croisée ({n_splits}-Fold) — AUC = {roc_auc:.3f}"
    else:
        fpr, tpr = np.array([0.0, 1.0]), np.array([0.0, 1.0])
        roc_auc = 0.5
        importances = np.zeros(len(feature_names))
        ml_status = "Cohorte quasi-homogène (ajustez la dose pour obtenir des répondeurs et non-répondeurs)"

    # Construction du DataFrame complet exportable en CSV
    cohort_df = pd.DataFrame({
        "Patient_ID": [f"VP_{i+1:03d}" for i in range(n_patients)],
        "Cancer_Type": profile.name,
        "Biomarker_Target_B0": np.round(biomarker_b0, 4),
        "Biomarker_PDL1_TMB": np.round(pdl1_tmb, 4),
        "Stromal_Permeability": np.round(stroma_vec, 4),
        "Proliferation_Lambda0": np.round(lambda0_vec, 5),
        "Initial_Resistant_Fraction": np.round(res_frac_vec, 4),
        "Clearance_k10": np.round(k10_vec, 4),
        "Baseline_Tumor_g": np.round(w0_vec, 4),
        "Final_Tumor_g": np.round(w_final, 4),
        "Final_Resistant_Clone_g": np.round(resistant_burden_final, 4),
        "Tumor_Change_Pct": np.round(pct_change, 2),
        "RECIST_1_1": recist_labels,
        "Responder_Binary": is_responder
    })

    return {
        "time_days": time_grid / 24.0,
        "w_median": np.median(w_trajectories, axis=1),
        "w_p05": np.percentile(w_trajectories, 5, axis=1),
        "w_p95": np.percentile(w_trajectories, 95, axis=1),
        "w_p25": np.percentile(w_trajectories, 25, axis=1),
        "w_p75": np.percentile(w_trajectories, 75, axis=1),
        "recist": recist_labels,
        "orr": float(np.mean(is_responder) * 100.0),
        "fpr": fpr,
        "tpr": tpr,
        "roc_auc": roc_auc,
        "feature_names": feature_names,
        "importances": importances,
        "ml_status": ml_status,
        "cohort_df": cohort_df
    }