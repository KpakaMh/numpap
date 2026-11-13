from dataclasses import dataclass
from typing import Tuple, List

@dataclass(frozen=True)
class CancerProfile:
    key: str
    rank_globocan: int
    new_cases_2024: str
    name: str
    subtype: str
    origin: str
    mutations: str
    biomarker_name: str
    invasiveness_desc: str
    metastatic_sites: str
    # Paramètres mécanistes QSP
    lambda0: float          # Croissance exponentielle initiale (1/h)
    lambda1: float          # Croissance linéaire asymptotique (g/h)
    stromal_permeability: float  # Pénétration du médicament dans le microenvironnement (0.15 à 1.0)
    base_potency_k2: float  # Sensibilité intrinsèque au traitement (L/mg/h)
    emt_invasiveness: float # Facteur d'activation EMT / MMP (0.2 à 2.5)
    initial_mass: float     # Masse tumorale initiale (g)
    # Coordonnées anatomiques sur le SVG patient (x, y) et sites métastatiques
    organ_coords: Tuple[int, int]
    metastasis_coords: List[Tuple[int, int]]
    is_hematologic: bool = False


CANCER_CATALOG = {
    "lung_nsclc": CancerProfile(
        key="lung_nsclc", rank_globocan=1, new_cases_2024="2 637 005",
        name="1. Cancer du Poumon (Adénocarcinome NSCLC)",
        subtype="Adénocarcinome bronchique / alvéolaire",
        origin="Cellules épithéliales bronchiques et alvéolaires",
        mutations="EGFR, ALK, KRAS G12C, TP53",
        biomarker_name="Score d'expression / mutation EGFR-ALK",
        invasiveness_desc="Infiltration des bronches, vaisseaux sanguins et plèvre.",
        metastatic_sites="Cerveau, foie, os, glandes surrénales",
        lambda0=0.014, lambda1=0.06, stromal_permeability=0.75,
        base_potency_k2=0.007, emt_invasiveness=1.5, initial_mass=0.4,
        organ_coords=(65, 82), metastasis_coords=[(50, 22), (45, 102), (35, 130)]
    ),
    "lung_sclc": CancerProfile(
        key="lung_sclc", rank_globocan=1, new_cases_2024="2 637 005",
        name="1b. Cancer du Poumon (À Petites Cellules - SCLC)",
        subtype="Carcinome neuroendocrine à petites cellules",
        origin="Cellules neuroendocrines épithéliales bronchiques",
        mutations="TP53, RB1 (inactivation quasi-systématique), MYC",
        biomarker_name="Expression DLL3 / Signature neuroendocrine",
        invasiveness_desc="Très agressif, franchissement précoce de la membrane basale.",
        metastatic_sites="Dissémination précoce : cerveau, foie, moelle osseuse",
        lambda0=0.024, lambda1=0.09, stromal_permeability=0.80,
        base_potency_k2=0.0085, emt_invasiveness=2.3, initial_mass=0.5,
        organ_coords=(65, 82), metastasis_coords=[(50, 22), (45, 102), (68, 125)]
    ),
    "breast_her2": CancerProfile(
        key="breast_her2", rank_globocan=2, new_cases_2024="2 434 087",
        name="2. Cancer du Sein (Sous-type HER2+)",
        subtype="Carcinome canalaire invasif HER2-amplifié",
        origin="Cellules épithéliales des canaux galactophores",
        mutations="Amplification ERBB2 (HER2), PIK3CA, TP53",
        biomarker_name="Amplification membranaire HER2 (Score IHC/FISH)",
        invasiveness_desc="Infiltration du stroma mammaire et envahissement lymphatique axillaire.",
        metastatic_sites="Ganglions axillaires, os, poumons, foie, cerveau",
        lambda0=0.013, lambda1=0.055, stromal_permeability=0.80,
        base_potency_k2=0.0095, emt_invasiveness=1.3, initial_mass=0.35,
        organ_coords=(36, 86), metastasis_coords=[(28, 76), (65, 82), (45, 102)]
    ),
    "breast_tnbc": CancerProfile(
        key="breast_tnbc", rank_globocan=2, new_cases_2024="2 434 087",
        name="2b. Cancer du Sein (Triple Négatif - TNBC)",
        subtype="Carcinome basal-like (ER-, PR-, HER2-)",
        origin="Cellules épithéliales canalaires / progéniteurs basaux",
        mutations="TP53, BRCA1/2, perte PTEN, voie PI3K/AKT",
        biomarker_name="Score PD-L1 (CPS) / Déficit recombinaison homologue (HRD)",
        invasiveness_desc="Forte transition épithélio-mésenchymateuse (EMT) et MMP élevées.",
        metastatic_sites="Poumons, cerveau, foie (dissémination hématogène rapide)",
        lambda0=0.019, lambda1=0.075, stromal_permeability=0.65,
        base_potency_k2=0.0055, emt_invasiveness=1.9, initial_mass=0.4,
        organ_coords=(36, 86), metastasis_coords=[(50, 22), (65, 82), (45, 102)]
    ),
    "colorectal": CancerProfile(
        key="colorectal", rank_globocan=3, new_cases_2024="2 041 007",
        name="3. Cancer Colorectal (Adénocarcinome)",
        subtype="Séquence Adénome → Dysplasie → Carcinome invasif",
        origin="Cellules épithéliales glandulaires de la muqueuse colique",
        mutations="APC, KRAS, TP53, SMAD4 ou instabilité microsatellitaire (dMMR/MSI)",
        biomarker_name="Statut MSI-H / Charge mutationnelle tumorale (TMB)",
        invasiveness_desc="Infiltration progressive de la sous-muqueuse, musculeuse et séreuse.",
        metastatic_sites="Foie (via veine porte) et poumons",
        lambda0=0.010, lambda1=0.045, stromal_permeability=0.70,
        base_potency_k2=0.0065, emt_invasiveness=1.1, initial_mass=0.45,
        organ_coords=(42, 122), metastasis_coords=[(45, 102), (65, 82)]
    ),
    "prostate": CancerProfile(
        key="prostate", rank_globocan=4, new_cases_2024="1 546 112",
        name="4. Cancer de la Prostate (Adénocarcinome Acinaire)",
        subtype="Adénocarcinome hormonodépendant (Score Gleason / ISUP)",
        origin="Cellules épithéliales glandulaires prostatiques",
        mutations="Signalisation Récepteur Androgènes (AR), fusion TMPRSS2-ERG, PTEN, BRCA2",
        biomarker_name="Activité transcriptionnelle AR / PSA cinétique",
        invasiveness_desc="Extension locale périprostatique, vésicules séminales et ganglions pelviens.",
        metastatic_sites="Os (métastases ostéocondensantes) et ganglions pelviens",
        lambda0=0.006, lambda1=0.030, stromal_permeability=0.75,
        base_potency_k2=0.0075, emt_invasiveness=0.7, initial_mass=0.3,
        organ_coords=(50, 136), metastasis_coords=[(35, 130), (65, 130)]
    ),
    "stomach": CancerProfile(
        key="stomach", rank_globocan=5, new_cases_2024="980 286",
        name="5. Cancer de l’Estomac (Adénocarcinome Gastrique)",
        subtype="Type intestinal ou diffus (cellules en bague à chaton)",
        origin="Cellules épithéliales de la muqueuse gastrique",
        mutations="CDH1 (perte E-cadhérine), TP53, HER2, inflammation H. pylori",
        biomarker_name="Expression Claudine 18.2 / Statut HER2-MSI",
        invasiveness_desc="Infiltration pariétale profonde (musculeuse et séreuse) avec perte d'adhérence.",
        metastatic_sites="Ganglions régionaux, péritoine, foie, ovaires",
        lambda0=0.013, lambda1=0.055, stromal_permeability=0.60,
        base_potency_k2=0.006, emt_invasiveness=1.5, initial_mass=0.4,
        organ_coords=(57, 105), metastasis_coords=[(45, 102), (50, 122)]
    ),
    "thyroid": CancerProfile(
        key="thyroid", rank_globocan=6, new_cases_2024="959 281",
        name="6. Cancer de la Thyroïde (Papillaire vs Anaplasique)",
        subtype="Carcinome papillaire différencié (profil BRAF)",
        origin="Cellules folliculaires thyroïdiennes",
        mutations="BRAF V600E, RAS, réarrangements RET/PTC",
        biomarker_name="Mutation BRAF V600E / Captation iodure (NIS)",
        invasiveness_desc="Dissémination lymphatique cervicale fréquente mais progression lente.",
        metastatic_sites="Ganglions cervicaux, poumons (formes avancées)",
        lambda0=0.005, lambda1=0.025, stromal_permeability=0.85,
        base_potency_k2=0.008, emt_invasiveness=0.5, initial_mass=0.25,
        organ_coords=(50, 44), metastasis_coords=[(42, 48), (65, 82)]
    ),
    "liver": CancerProfile(
        key="liver", rank_globocan=7, new_cases_2024="843 045",
        name="7. Cancer du Foie (Carcinome Hépatocellulaire - CHC)",
        subtype="CHC sur microenvironnement cirrhotique / fibrotique",
        origin="Hépatocytes transformés sur terrain inflammatoire chronique",
        mutations="Promoteur TERT, TP53, CTNNB1 (voie Wnt/β-caténine), VEGF élevé",
        biomarker_name="Signature angiogénique VEGF / Alpha-fœtoprotéine (AFP)",
        invasiveness_desc="Invasion vasculaire fréquente (thrombose tumorale de la veine porte).",
        metastatic_sites="Poumons, os, ganglions lymphatiques",
        lambda0=0.014, lambda1=0.06, stromal_permeability=0.55,
        base_potency_k2=0.0055, emt_invasiveness=1.4, initial_mass=0.5,
        organ_coords=(44, 102), metastasis_coords=[(65, 82), (35, 130)]
    ),
    "bladder": CancerProfile(
        key="bladder", rank_globocan=8, new_cases_2024="635 264",
        name="8. Cancer de la Vessie (Carcinome Urothélial)",
        subtype="Carcinome urothélial infiltrant le muscle (TVIM)",
        origin="Cellules urothéliales de la muqueuse vésicale",
        mutations="FGFR3, TP53, RB1, ERBB2, charge mutationnelle élevée",
        biomarker_name="Altération FGFR3 / Expression Nectine-4 & PD-L1",
        invasiveness_desc="Pénétration du muscle détrusor puis des tissus graisseux périvésicaux.",
        metastatic_sites="Ganglions pelviens, os, poumons, foie",
        lambda0=0.012, lambda1=0.05, stromal_permeability=0.70,
        base_potency_k2=0.0068, emt_invasiveness=1.2, initial_mass=0.35,
        organ_coords=(50, 130), metastasis_coords=[(35, 130), (45, 102), (65, 82)]
    ),
    "cervix": CancerProfile(
        key="cervix", rank_globocan=9, new_cases_2024="604 196",
        name="9. Cancer du Col de l’Utérus (Carcinome Épidermoïde HPV+)",
        subtype="Oncogenèse virale persistante HPV 16 / 18",
        origin="Cellules épithéliales de la zone de jonction du col utérin",
        mutations="Oncoprotéines virales E6 (dégradation p53) & E7 (inactivation pRb), PIK3CA",
        biomarker_name="Charge virale E6/E7 HPV16-18 / Expression PD-L1",
        invasiveness_desc="Lésion intraépithéliale → carcinome in situ → invasion du stroma cervical.",
        metastatic_sites="Structures pelviennes, ganglions iliaques/para-aortiques, poumons",
        lambda0=0.011, lambda1=0.048, stromal_permeability=0.75,
        base_potency_k2=0.0072, emt_invasiveness=1.1, initial_mass=0.35,
        organ_coords=(50, 126), metastasis_coords=[(40, 124), (65, 82)]
    ),
    "lymphoma": CancerProfile(
        key="lymphoma", rank_globocan=10, new_cases_2024="562 793",
        name="10. Lymphome Non Hodgkinien (DLBCL)",
        subtype="Lymphome diffus à grandes cellules B (lignée mobile)",
        origin="Lymphocytes B du centre germinatif (cellules hématopoïétiques)",
        mutations="Réarrangements MYC, BCL2, BCL6 (lymphome double-hit)",
        biomarker_name="Densité antigénique CD20 / Index prolifératif Ki-67",
        invasiveness_desc="Cellules intrinsèquement mobiles : dissémination ganglionnaire et extranodale.",
        metastatic_sites="Chaînes ganglionnaires multiples, rate, moelle osseuse, SNC",
        lambda0=0.020, lambda1=0.08, stromal_permeability=0.95,
        base_potency_k2=0.011, emt_invasiveness=1.8, initial_mass=0.45,
        organ_coords=(38, 72), metastasis_coords=[(62, 72), (58, 102), (42, 124)],
        is_hematologic=True
    ),
    "pancreas": CancerProfile(
        key="pancreas", rank_globocan=11, new_cases_2024="531 318",
        name="11. Cancer du Pancréas (Adénocarcinome Canalaire - PDAC)",
        subtype="Adénocarcinome canalaire avec stroma desmoplastique dense",
        origin="Cellules épithéliales canalaires pancréatiques",
        mutations="KRAS (>90%), TP53, CDKN2A, SMAD4",
        biomarker_name="Perméabilité stromale / Signature KRAS-SMAD4 & CA 19-9",
        invasiveness_desc="Stroma fibrotique dense bloquant la diffusion vasculaire ; invasion périneurale précoce.",
        metastatic_sites="Foie, péritoine, poumons",
        lambda0=0.018, lambda1=0.075, stromal_permeability=0.28,
        base_potency_k2=0.0045, emt_invasiveness=2.1, initial_mass=0.45,
        organ_coords=(54, 111), metastasis_coords=[(45, 102), (50, 122), (65, 82)]
    )
}