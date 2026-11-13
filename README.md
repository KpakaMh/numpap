# Numpap — Prototype Exploratoire de Modélisation QSP & PK/PD en Oncologie Translationnelle

> **Statut du projet & Contexte d'Utilisation (*Context of Use* — Cadre FDA MIDD) :**  
> **Numpap** est un **prototype exploratoire de pharmacologie quantitive des systèmes (QSP) et de modélisation PK/PD mécaniste**. Il est conçu pour l'exploration d'hypothèses biologiques, l'analyse des relations exposition-réponse, la vérification numérique des bilans de masse et l'analyse de sensibilité globale. **Il ne constitue pas un jumeau numérique clinique validé** et ne doit pas être utilisé à des fins de décision médicale ou de pronostic individuel.

---

## 1. Philosophie Scientifique et Stratégie de Validation

Plutôt que d'appliquer une structure générique uniforme à de multiples localisations tumorales, l'architecture de **Numpap** adopte une démarche incrémentale rigoureuse :

1. **Modèle de référence unique et documenté :** Calibration et vérification sur un couple pathologie-molécule bien caractérisé — l'**adénocarcinome pulmonaire non à petites cellules (CBNPC / NSCLC) présentant une mutation activatrice d'*EGFR* ($\text{Ex19del}$ ou $\text{L858R}$)** traité par **Osimertinib** (inhibiteur de tyrosine kinase de 3e génération) et son métabolite actif **AZ5104**.
2. **Harmonisation stricte entre simulation individuelle ($N = 1$) et cohorte virtuelle (VPop) :** Les simulations d'un patient unique et d'une population virtuelle exécutent exactement le même système d'équations différentielles ordinaires (EDO) et les mêmes hypothèses biophysiques.
3. **Dissociation entre biomarqueurs mécanistes exploratoires et classifications cliniques :** Les variables internes du modèle (masses viables/apoptotiques, fractions clonales, ratios immunitaires) sont strictement séparées des référentiels cliniques (stadification anatomique TNM/AJCC et évaluation radiologique multi-lésionnelle RECIST 1.1).

---

## 2. Architecture Mathématique et Mécanistique (`EGFRNSCLCQSPEngine`)

### 2.1 Pharmacocinétique Bicompartimentale, Métabolite Actif et Bilans de Masse

La pharmacocinétique de la molécule mère (Osimertinib) est décrite par un modèle bicompartimental physiologique ($V_1, V_2, CL, Q$) couplé à un compartiment de distribution et d'élimination pour son métabolite actif N-déméthylé **AZ5104** ($V_{\text{met}}, CL_{\text{met}}$) avec prise en compte du rapport stoechiométrique des masses molaires ($\text{MW}_{\text{met}} / \text{MW}_{\text{parent}}$) :

$$
\begin{aligned}
\frac{dC_1}{dt} &= -\frac{CL + Q}{V_1} C_1 + \frac{Q}{V_1} C_2 + \frac{R_{\text{in}}(t)}{V_1} \\
\frac{dC_2}{dt} &= \frac{Q}{V_2} (C_1 - C_2) \\
\frac{dC_{\text{met}}}{dt} &= \frac{1}{V_{\text{met}}} \left( f_m \cdot CL \cdot C_1 \cdot \frac{\text{MW}_{\text{met}}}{\text{MW}_{\text{parent}}} - CL_{\text{met}} \cdot C_{\text{met}} \right)
\end{aligned}
$$

* **Discontinuités exactes des bolus :** Les administrations intraveineuses ou per os instantanées sont traitées par découpage exact des intervalles d'intégration aux instants $t_i$ des doses, avec mise à jour de la condition initiale $C_1(t_i^+) = C_1(t_i^-) + \text{Dose}_i / V_1$.
* **Conservation des flux PK :** Le moteur intègre simultanément l'élimination parentale non métabolisée $(1 - f_m) CL \cdot C_1$, la fraction biotransformée $f_m CL \cdot C_1$ et l'élimination du métabolite $CL_{\text{met}} C_{\text{met}}$, garantissant des résidus de bilan de masse inférieurs à $10^{-8}\ \text{mg}$.

### 2.2 Cascade d'Exposition : Totale, Libre Plasmatique et Libre Intratumorale

Conformément à l'hypothèse du médicament libre (*Free Drug Hypothesis*), seule la fraction non liée aux protéines plasmatiques ($f_u$) traverse les barrières biologiques et interagit avec le domaine kinase d'*EGFR*. Le modèle distingue trois niveaux d'exposition :

1. **Concentration plasmatique totale :** $C_1$ et $C_{\text{met}}$ exprimées en $\text{mg/L}$.
2. **Concentration plasmatique libre effective ($\mu\text{M}$) :**
   $$C_{u,\text{plasma,eff}}(t) = \frac{f_{u,\text{parent}} \cdot C_1(t) \cdot 1000}{\text{MW}_{\text{parent}}} + \rho_{\text{met}} \cdot \frac{f_{u,\text{met}} \cdot C_{\text{met}}(t) \cdot 1000}{\text{MW}_{\text{met}}}$$
   où $\rho_{\text{met}}$ représente la puissance relative libre du métabolite AZ5104 par rapport à la molécule mère.
3. **Concentration libre intratumorale ($C_{u,\text{tumor}}$ en $\mu\text{M}$) :** Modulée par le coefficient de partition libre $K_{p,uu}$, la barrière stromale desmoplastique ($b_{\text{stroma}}$) et l'index de perfusion vasculaire fonctionnelle $V_{\text{vasc}}(t)$ :
   $$C_{u,\text{tumor}}(t) = C_{u,\text{plasma,eff}}(t) \cdot K_{p,uu,\text{base}} \cdot (1 - b_{\text{stroma}}) \cdot \text{clip}\!\left(\frac{V_{\text{vasc}}(t)}{0.35 + V_{\text{vasc}}(t)},\ 0.15,\ 1.15\right)$$

### 2.3 Dynamique Tumorale à 3 États ($x_S \leftrightarrow x_P \to x_R$) et Séparation Viable / Apoptotique

Pour dépasser la limite d'une simple transition binaire sensible-résistant, la masse tumorale viable ($W_{\text{viable}}$) est structurée en trois sous-populations phénotypiques et génétiques, complétées par une cascade de trois compartiments de transit apoptotique ($W_{\text{apop}} = x_{A1} + x_{A2} + x_{A3}$) :

* **$x_S(t)$ — Cellules sensibles (*EGFR* $\text{Ex19del}/\text{L858R}$) :** Prolifèrent au taux $\lambda_0$ et subissent une lyse saturable de Hill pilotée par $C_{u,\text{tumor}}$ et la sensibilité de la cible $B_0$ ($EC_{50,u,\text{eff}} = EC_{50,u} / B_0$).
* **$x_P(t)$ — Cellules persistantes tolérantes (*Drug-Tolerant Persisters*, DTP) :** Population à cycle lent induite de manière épigénétique sous pression thérapeutique ($k_{SP} \cdot \text{Hill}(C_{u,\text{tumor}})$) et capable de **réversion** vers l'état sensible ($k_{PS} \cdot (1 - \text{Hill}(C_{u,\text{tumor}}))$) lors du *washout* pharmacologique.
* **$x_R(t)$ — Clone génétiquement résistant (ex. mutation secondaire *EGFR* $\text{C797S}$, amplification *MET*) :** Émerge de manière irréversible à partir de $x_S$ ($\mu_{SR}$) et, à un taux accru, à partir du réservoir persistant $x_P$ ($\mu_{PR}$).

$$
\begin{aligned}
\frac{dx_S}{dt} &= g_S - k_{\text{kill},S} x_S - k_{\text{imm},S} x_S - \Phi_{S \to P} + \Phi_{P \to S} - \Phi_{S \to R} \\
\frac{dx_P}{dt} &= g_P - k_{\text{kill},P} x_P - k_{\text{imm},P} x_P + \Phi_{S \to P} - \Phi_{P \to S} - \Phi_{P \to R} \\
\frac{dx_R}{dt} &= g_R - k_{\text{kill},R} x_R - k_{\text{imm},R} x_R + \Phi_{S \to R} + \Phi_{P \to R} \\
\frac{dx_{A1}}{dt} &= \text{Kill}_{\text{total}} - k_1 x_{A1}, \qquad \frac{dx_{A2}}{dt} = k_1(x_{A1} - x_{A2}), \qquad \frac{dx_{A3}}{dt} = k_1(x_{A2} - x_{A3})
\end{aligned}
$$

La masse tumorale totale suivie en imagerie vaut $W_{\text{total}}(t) = W_{\text{viable}}(t) + W_{\text{apop}}(t)$, tandis que seule la fraction viable $W_{\text{viable}}(t) = x_S(t) + x_P(t) + x_R(t)$ contribue à la prolifération et à la sécrétion de facteurs angiogéniques.

### 2.4 Microenvironnement Tumoral (TME) et Mécanismes d'Association

Chaque classe thérapeutique agit sur son compartiment physiologique propre :

* **Compartiment vasculaire ($V_{\text{vasc}}$) et Anti-VEGF :** La masse viable stimule l'angiogenèse via la sécrétion de VEGF ($k_{\text{angio}}$), modulant la capacité de charge tumorale effective $K_{\text{cap}}(V_{\text{vasc}})$ ainsi que la délivrance intratumorale $K_{p,uu,\text{eff}}$. Le traitement anti-VEGF inhibe ce signal angiogénique.
* **Immunité anti-tumorale ($\text{CD8}^+$) vs Immunosuppression ($\text{T}_{\text{reg}}$ / TAM-M2) :**
  * L'hypoxie induite par une faible vascularisation favorise le recrutement de cellules immunosuppressives $I_{\text{sup}}$ ($\text{T}_{\text{reg}}$ et macrophages associés aux tumeurs M2), qui freinent l'activité cytolytique et l'expansion lymphocytaire via le facteur $(1 + 1.5 I_{\text{sup}})^{-1}$.
  * Les inhibiteurs de points de contrôle immunitaire (anti-PD-1 / anti-PD-L1) lèvent l'épuisement des lymphocytes effecteurs $T_{\text{CD8}}$ proportionnellement au score $\text{TMB}/\text{PD-L1}$ et à la libération d'antigènes tumoraux issue de la clairance apoptotique ($k_1 x_{A3}$). Lorsque l'immunothérapie est désactivée ($\text{anti-PD-1} = 0$), l'effet imputable à l'ICI est strictement nul.

---

## 3. Dissociation des Indicateurs Exploratoires et des Référentiels Cliniques

Afin d'éviter toute sur-interprétation clinique des sorties simulées :

* **Indicateurs exploratoires produits par le moteur :** Masse viable $W_{\text{viable}}$ ($\text{g}$), masse apoptotique $W_{\text{apop}}$ ($\text{g}$), fraction clonale résistante $f_R = x_R / W_{\text{viable}}$, fraction persistante $f_P = x_P / W_{\text{viable}}$, index vasculaire $V_{\text{vasc}}$, ratio immunitaire $T_{\text{CD8}} / I_{\text{sup}}$, et diamètre sphérique équivalent unidimensionnel $d_{\text{eq}} = 20 \left(\frac{3 W_{\text{total}}}{4\pi}\right)^{1/3}$ ($\text{mm}$).
* **Limites cliniques explicites :**
  * Le modèle **ne calcule pas** de stadification clinique TNM/AJCC (qui requiert une évaluation anatomique de l'extension locale, ganglionnaire et métastatique à distance).
  * La variation relative du diamètre équivalent ($\Delta d_{\text{eq}} \le -30\%$) est utilisée uniquement comme critère unidimensionnel exploratoire *in silico* et **ne remplace pas** une évaluation RECIST 1.1 radiologique sur lésions cibles et non-cibles.

---

## 4. Vérification Numérique, Analyse de Sensibilité et Machine Learning

### 4.1 Suite Automatisée de Vérification Formelle (`pytest`)

Avant toute confrontation à des données expérimentales, le moteur est soumis à une batterie de **12 tests de vérification** (conformes aux bonnes pratiques ASME V&V 40) couvrant :

1. **Sans médicament :** Conformité exacte de $W_{\text{viable}}(t)$ au modèle de croissance non traité de référence ($W_{\text{apop}} = 0$).
2. **Dose nulle ($D = 0, R_{\text{in}} = 0$) :** Absence totale d'effet pharmacologique attribuable au médicament.
3. **Bolus isolé :** Vérification de la trajectoire PK face à la solution analytique bi-exponentielle exacte ($A e^{-\alpha t} + B e^{-\beta t}$).
4. **Clairance augmentée ($CL \uparrow$) :** Diminution stricte de l'exposition ($\text{AUC}$) et de l'efficacité antitumorale, toutes choses égales par ailleurs.
5. **Résistance initiale nulle ($x_R(0) = 0, \mu = 0$) :** Absence de création artificielle de cellules résistantes.
6. **Clone résistant sans avantage sélectif :** Conservation exacte du ratio clonal en l'absence de pression de sélection.
7. **Immunothérapie désactivée ($\text{anti-PD-1} = 0$) :** Activation lymphocytaire $T_{\text{CD8}}(t) = 0$ sur toute la simulation.
8. **Anti-VEGF désactivé ($\text{anti-VEGF} = 0$) :** Retour exact à la dynamique vasculaire basale.
9. **Tous traitements désactivés :** Identité stricte avec la trajectoire tumorale non traitée.
10. **Absence de cellules viables initiales ($W_{\text{viable}}(0) = 0$) :** Absence de génération spontanée de masse tumorale.
11. **Conservation des bilans de masse PK :** Vérification des résidus parentaux et métaboliques ($< 10^{-8}\ \text{mg}$).
12. **Convergence numérique multi-solveurs :** Accord entre l'intégrateur explicite `RK45` et le solveur implicite pour systèmes raides `Radau` ($\text{rtol} = 10^{-9}, \text{atol} = 10^{-12}$).

### 4.2 Cohorte Virtuelle (VPop), Analyse de Sensibilité PRCC et ML

* **Analyse de sensibilité globale (PRCC) :** Le module calcule les **coefficients de corrélation partielle sur les rangs (PRCC)** entre les paramètres échantillonnés ($CL, f_u, b_{\text{stroma}}, B_0, \text{TMB}/\text{PD-L1}, \lambda_0, x_{R0}$) et la masse tumorale viable finale afin de quantifier le poids mécanique de chaque source d'incertitude.
* **Pipeline Machine Learning sans fuite de données :** La prédiction de la réponse unidimensionnelle utilise un `Pipeline(StandardScaler(), RandomForestClassifier())` ajusté **strictement à l'intérieur de chaque pli** de validation croisée stratifiée (`StratifiedKFold`). L'importance des variables (indice de Gini) est présentée comme une mesure d'association prédictive au sein du modèle et **non comme une preuve de causalité biologique**.

---

## 5. Structure du Dépôt et Installation

```text
Numpap/
├── app.py                                 # Interface interactive Dash (Patient N=1 & Cohorte VPop)
├── core/
│   ├── egfr_nsclc_qsp_engine.py           # Moteur QSP/PK-PD de référence (CBNPC EGFR / Osimertinib)
│   ├── patient_pkpd_engine.py             # Moteur PK/PD généraliste & simulateur de cohorte
│   ├── drug_catalog.py                    # Fiches d'identité PK/PD et seuils libres (EC50,u / Ctox,u)
│   └── cancer_catalog.py                  # Paramètres biologiques des profils tumoraux
├── services/
│   └── patient_logger.py                  # Journalisation des trajectoires temporelles
├── ui/
│   └── svg_patient_tumor.py               # Visualisation anatomique SVG du microenvironnement
├── tests/
│   └── test_consortium_verification.py    # Suite pytest des 12 tests de vérification formelle
└── README.md
```

### Prérequis et Exécution (Windows PowerShell)

```powershell
# 1. Créer et activer l'environnement virtuel
python -m venv .venv
.\.venv\Scripts\Activate.ps1

# 2. Installer les dépendances scientifiques et d'interface
pip install numpy scipy pandas scikit-learn dash plotly pytest

# 3. Exécuter la suite complète de vérification numérique et biologique
python -m pytest -v

# 4. Lancer l'interface interactive locale
python app.py
```

---

## 6. Feuille de Route Scientifique

* **Étape 1 (En cours) — Vérification formelle et mécanique :** Validation des bilans de masse parent/métabolite, des conditions limites et de la sensibilité PRCC sur le modèle de référence *EGFR*-NSCLC.
* **Étape 2 — Calibration sur données publiées :** Estimation des paramètres ($\lambda_0, E_{\max}, EC_{50,u}, k_{SP}, \mu_{PR}$) à partir de courbes de croissance de xénogreffes dérivées de patients (PDX *EGFR* $\text{Ex19del}/\text{L858R}$) et de profils PK/PD cliniques publiés pour l'Osimertinib (essais AURA / FLAURA).
* **Étape 3 — Validation externe et quantification d'incertitude :** Confrontation des prédictions de survie sans progression simulée et d'émergence de résistance à une cohorte clinique indépendante non utilisée lors de la calibration.
* **Étape 4 — Extension modulaire multi-pathologies :** Déploiement progressif de modules spécifiques aux histologies complexes (ex. barrière desmoplastique dense pour l'adénocarcinome canalaire pancréatique PDAC, hématopoïèse pour les hémopathies malignes).
