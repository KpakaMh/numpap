# 🧬 QSP Oncology Digital Twin — Plateforme de Pharmacologie Quantitative & Médecine Computationnelle

![Python](https://img.shields.io/badge/Python-3.11+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Dash](https://img.shields.io/badge/Dash_%26_Plotly-Interactive_UI-008DE4?style=for-the-badge&logo=plotly&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy_%26_NumPy-ODE_Solvers-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white)
![Scikit-Learn](https://img.shields.io/badge/Scikit--Learn-Biomarker_ML-F7931E?style=for-the-badge&logo=scikitlearn&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-Reproducible-2496ED?style=for-the-badge&logo=docker&logoColor=white)

**QSP Oncology Digital Twin** (`Numpap`) est une plateforme interactive de **Pharmacologie Quantitative des Systèmes (QSP)** et de simulation d'essais cliniques *in silico*. Elle couple un système mécaniste d'**équations différentielles ordinaires (EDO)** à un pipeline de **Machine Learning** pour modéliser la progression tumorale, l'émergence de résistances clonales, l'effet de trithérapies et l'identification de biomarqueurs prédictifs sur les **11 localisations cancéreuses majeures (données GLOBOCAN 2024 / CIRC-IARC)**.

---

## 🎯 Objectifs Scientifiques & Fonctionnalités Clés

### 1. Jumeau Numérique Patient ($N = 1$) & Microenvironnement Dynamique
* **Catalogue Oncologique Multi-Cancers (GLOBOCAN 2024) :** 13 sous-types moléculaires couvrant les 11 cancers les plus fréquents au monde, paramétrés selon leur cinétique de croissance ($\lambda_0, \lambda_1$), leur perméabilité stromale ($\eta_{\text{stroma}}$), leurs mutations pilotes (`EGFR`, `HER2`, `KRAS`, `TP53`, `BRAF`, `HPV E6/E7`...) et leur tropisme métastatique.
* **Modélisation des 7 Stades de Progression Tumorale :** Suivi en temps réel depuis l'*Initiation tumorale* (Stade 1) jusqu'à l'*Invasion locale EMT/MMP* (Stade 4), l'*Intravasation* (Stade 5) et la *Dissémination métastatique* (Stade 7).
* **Combinaisons Thérapeutiques (Trithérapie *In Silico*) :**
  * **Thérapie ciblée / Chimiothérapie :** Bolus IV et perfusion continue, modulés par la barrière stromale desmoplastique et l'expression du biomarqueur cible $B_0$.
  * **Anti-angiogénique (Anti-VEGF) :** Inhibition de la néovascularisation tumorale réduisant la capacité de croissance asymptotique ($\lambda_1$).
  * **Immunothérapie (Inhibiteur PD-1 / PD-L1) :** Activation des lymphocytes T cytotoxiques effecteurs $E(t)$ selon le score PD-L1 / TMB, capables d'éliminer les sous-clones résistants.
* **Imagerie Vectorielle Anatomique (SVG/CSS) :** Localisation dynamique de la tumeur primaire sur la silhouette anatomique, épaisseur visuelle de la barrière stromale, spicules d'invasion EMT et apparition des foyers métastatiques secondaires.

### 2. Essai Clinique Virtuel (VPop) & Découverte de Biomarqueurs par Machine Learning
* **Génération de Cohortes Virtuelles ($N = 100$ à $400$ patients) :** Échantillonnage Monte-Carlo sous contraintes physiologiques intégrant la variabilité inter-individuelle PK ($k_{10}, V_c$), tumorale ($\lambda_0, W_0$), microenvironnementale ($\eta_{\text{stroma}}$) et clonale ($x_{R,0}$).
* **Solveur RK4 Vectorisé & Stratification RECIST 1.1 :** Intégration simultanée de l'ensemble de la cohorte et classification clinique automatique en **CR** (*Réponse Complète*), **PR** (*Réponse Partielle*), **SD** (*Maladie Stable*) et **PD** (*Progression*), avec calcul du taux de réponse objective (**ORR**).
* **Pipeline Machine Learning (`scikit-learn`) :** Entraînement d'un classifieur **Random Forest** avec validation croisée stratifiée (*Stratified K-Fold*), génération de la **courbe ROC-AUC** et extraction de l'**importance relative des biomarqueurs prédictifs** (*Feature Importance*).
* **Export de Données :** Téléchargement direct de la cohorte simulée au format `.CSV` pour des analyses biostatistiques complémentaires sous R ou Python.

---

## 🧮 Formulation Mathématique (Système à 8 EDO)

Le moteur couvre la pharmacocinétique bi-compartimentale ($C_c, C_p$), la croissance tumorale de Simeoni scindée en populations sensible ($x_S$) et résistante ($x_R$), les compartiments de transit apoptotique ($x_2, x_3, x_4$) et la réponse immunitaire cytotoxique ($E$) :

$$\begin{aligned} \frac{dC_c}{dt} &= -(k_{10} + k_{12}) C_c + k_{21} C_p + \frac{R_{\text{inf}}}{V_c}, \qquad \frac{dC_p}{dt} = k_{12} C_c - k_{21} C_p \\ \frac{dx_S}{dt} &= \frac{\lambda_0 x_S}{\left[1 + \left(\frac{\lambda_0}{\lambda_{1,\text{eff}}} W\right)^\psi\right]^{1/\psi}} - k_{2,S} C_c x_S - k_{\text{imm}} E x_S - \mu(C_c) x_S \\ \frac{dx_R}{dt} &= \frac{0.88\,\lambda_0 x_R}{\left[1 + \left(\frac{\lambda_0}{\lambda_{1,\text{eff}}} W\right)^\psi\right]^{1/\psi}} - k_{2,R} C_c x_R - k_{\text{imm}} E x_R + \mu(C_c) x_S \\ \frac{dx_2}{dt} &= \text{Kill}_{\text{total}} - k_1 x_2, \qquad \frac{dx_3}{dt} = k_1 (x_2 - x_3), \qquad \frac{dx_4}{dt} = k_1 (x_3 - x_4) \\ \frac{dE}{dt} &= \alpha_{\text{ICI}} \cdot \text{ICI} \cdot \left(\frac{\text{TMB}}{0.8 + \text{TMB}}\right)(1 + 1.5 x_2) - \delta_E (E - E_0) \end{aligned}$$

où la puissance thérapeutique effective intègre la barrière stromale et la saturation du biomarqueur cible $B_0$ :
$$k_{2,S} = k_{2,\text{base}} \cdot \eta_{\text{stroma}} \cdot \left(\frac{2 B_0}{0.5 + B_0}\right), \qquad W(t) = x_S + x_R + x_2 + x_3 + x_4$$

---

## 📊 Catalogue des 11 Pathologies Tumorales Intégrées (GLOBOCAN 2024)

| Rang | Localisation | Sous-types Modélisés | Altérations Moléculaires Clés | Spécificité Microenvironnementale |
| :---: | :--- | :--- | :--- | :--- |
| **1** | **Poumon** | Adénocarcinome (NSCLC) & Petites Cellules (SCLC) | `EGFR`, `ALK`, `KRAS`, `TP53`, `RB1` | Dissémination précoce (SCLC) |
| **2** | **Sein** | Sous-type `HER2+` & Triple Négatif (`TNBC`) | `ERBB2`, `BRCA1/2`, `PIK3CA`, `PTEN` | Forte transition EMT & MMP (TNBC) |
| **3** | **Côlon-Rectum** | Adénocarcinome colique | `APC`, `KRAS`, `TP53`, `dMMR/MSI` | Séquence adénome $\rightarrow$ carcinome |
| **4** | **Prostate** | Adénocarcinome acinaire | Signalisation `AR`, `TMPRSS2-ERG`, `PTEN` | Progression lente, tropisme osseux |
| **5** | **Estomac** | Adénocarcinome intestinal / diffus | `CDH1` (E-cadhérine), `TP53`, `HER2` | Perte d'adhérence épithéliale |
| **6** | **Thyroïde** | Carcinome papillaire | `BRAF V600E`, `RAS`, `RET/PTC` | Forte perméabilité, diffusion lymphatique |
| **7** | **Foie** | Carcinome hépatocellulaire (CHC) | `TERT`, `TP53`, `CTNNB1`, `VEGF` | Stroma fibrotique / cirrhotique |
| **8** | **Vessie** | Carcinome urothélial (TVIM) | `FGFR3`, `TP53`, `RB1`, `ERBB2` | Infiltration du muscle détrusor |
| **9** | **Col de l'utérus** | Carcinome épidermoïde `HPV+` | Oncoprotéines virales `E6/E7` (`p53`/`pRb`) | Oncogenèse virale HPV 16/18 |
| **10** | **Lymphome (LNH)** | Lymphome diffus à grandes cellules B (DLBCL) | `MYC`, `BCL2`, `BCL6` | Lignée mobile (perméabilité 95 %) |
| **11** | **Pancréas** | Adénocarcinome canalaire (PDAC) | `KRAS`, `TP53`, `CDKN2A`, `SMAD4` | Stroma desmoplastique dense ($\eta = 0.28$) |

---

## 🏗️ Architecture Modulaire du Projet

```text
Numpap/
├── core/
│   ├── __init__.py
│   ├── cancer_catalog.py        # Base de connaissances des 11 cancers GLOBOCAN 2024 & paramètres QSP
│   └── patient_pkpd_engine.py   # Solveur EDO individuel (8 compartiments) + Cohorte RK4 & Pipeline ML
├── services/
│   ├── __init__.py
│   └── patient_logger.py        # Journalisation des séries temporelles PK/PD, clonales et immunitaires
├── ui/
│   ├── __init__.py
│   └── svg_patient_tumor.py     # Rendu anatomique multi-organes & microenvironnement stromal (SVG/CSS)
├── app.py                       # Interface web interactive Dash (Mode Individuel & Essai Clinique Virtuel)
├── requirements.txt             # Dépendances Python
├── Dockerfile                   # Image Docker légère (python:3.11-slim)
├── docker-compose.yml           # Orchestration du service web
└── README.md
