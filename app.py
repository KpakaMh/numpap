import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objs as go
import numpy as np

from core.cancer_catalog import CANCER_CATALOG
from core.patient_pkpd_engine import PatientPKPDEngine, run_virtual_cohort_trial
from services.patient_logger import PatientLogger
from ui.svg_patient_tumor import generate_patient_oncology_svg

engine = PatientPKPDEngine(cancer_key="lung_nsclc")
logger = PatientLogger()
logger.log(engine.state)
latest_cohort_df = None

app = dash.Dash(__name__, title="QSP Oncology Digital Twin & ML Platform")

cancer_options = [
    {"label": f"#{p.rank_globocan} — {p.name} ({p.new_cases_2024} cas/an)", "value": k}
    for k, p in CANCER_CATALOG.items()
]

CARD_STYLE = {
    "backgroundColor": "#1e293b", "padding": "15px", "borderRadius": "10px",
    "border": "1px solid #334155"
}

app.layout = html.Div(style={
    "backgroundColor": "#0b0f19", "color": "#f8fafc", "minHeight": "100vh",
    "padding": "16px", "fontFamily": "Segoe UI, sans-serif"
}, children=[

    # En-tête Global
    html.Div(style={
        "display": "flex", "justifyContent": "space-between", "alignItems": "center",
        "borderBottom": "1px solid #1e293b", "paddingBottom": "12px", "marginBottom": "14px"
    }, children=[
        html.Div([
            html.H1("Plateforme QSP & Jumeau Numérique d'Oncologie Translationnelle", style={"margin": 0, "fontSize": "21px", "color": "#38bdf8"}),
            html.P("Modélisation PK/PD Multi-Cancers (11 localisations GLOBOCAN 2024), Résistance Clonale, Trithérapie & Machine Learning", style={"margin": "4px 0 0 0", "color": "#94a3b8", "fontSize": "12.5px"})
        ]),
        html.Div(style={"width": "430px"}, children=[
            html.Label("Pathologie Tumorale (CIRC / GLOBOCAN 2024) :", style={"fontSize": "12px", "color": "#cbd5e1", "fontWeight": "bold"}),
            dcc.Dropdown(
                id="dropdown-cancer", options=cancer_options, value="lung_nsclc", clearable=False,
                style={"color": "#0f172a", "fontSize": "13px", "marginTop": "4px"}
            )
        ])
    ]),

    dcc.Tabs(id="main-tabs", value="tab-individual", colors={
        "border": "#334155", "primary": "#38bdf8", "background": "#1e293b"
    }, children=[

        # ==================== ONGLET 1 : JUMEAU INDIVIDUEL ====================
        dcc.Tab(label="1. Jumeau Numérique Patient (N = 1) — Résistance Clonale & Combinaisons", value="tab-individual",
                style={"color": "#94a3b8", "backgroundColor": "#1e293b", "padding": "10px"},
                selected_style={"color": "#38bdf8", "backgroundColor": "#0f172a", "fontWeight": "bold", "padding": "10px"},
                children=[
            html.Div(style={"display": "grid", "gridTemplateColumns": "340px 340px 1fr", "gap": "15px", "marginTop": "15px"}, children=[

                # Colonne 1 : Paramètres Thérapeutiques & Biologiques
                html.Div(style={"display": "flex", "flexDirection": "column", "gap": "12px"}, children=[
                    html.Div(style=CARD_STYLE, children=[
                        html.H3("Protocole & Microenvironnement", style={"marginTop": 0, "color": "#38bdf8", "fontSize": "14px"}),
                        html.Div(style={"display": "flex", "gap": "8px", "marginBottom": "12px"}, children=[
                            html.Button("Bolus 120 mg", id="btn-bolus", n_clicks=0, style={"backgroundColor": "#0284c7", "color": "white", "border": "none", "padding": "8px", "borderRadius": "6px", "cursor": "pointer", "fontWeight": "bold", "flex": 1}),
                            html.Button("Démarrer", id="btn-toggle", n_clicks=0, style={"backgroundColor": "#10b981", "color": "white", "border": "none", "padding": "8px", "borderRadius": "6px", "cursor": "pointer", "fontWeight": "bold", "flex": 1}),
                            html.Button("Reset", id="btn-reset", n_clicks=0, style={"backgroundColor": "#ef4444", "color": "white", "border": "none", "padding": "8px", "borderRadius": "6px", "cursor": "pointer", "fontWeight": "bold"})
                        ]),

                        html.Label("Biomarqueur Cible B0 (Sensibilité thérapie primaire)", style={"color": "#cbd5e1", "fontSize": "11.5px"}),
                        dcc.Slider(id="slider-biomarker", min=0.1, max=3.0, step=0.1, value=1.0, marks={0.1: "0.1", 1.0: "1.0", 3.0: "3.0"}),

                        html.Label("Fraction Clonale Résistante Initiale xR0 (%)", style={"color": "#fca5a5", "fontSize": "11.5px"}),
                        dcc.Slider(id="slider-res-frac", min=0, max=20, step=2, value=2, marks={0: "0%", 10: "10%", 20: "20%"}),

                        html.Label("Anti-Angiogénique (Anti-VEGF)", style={"color": "#38bdf8", "fontSize": "11.5px"}),
                        dcc.Slider(id="slider-vegf", min=0.0, max=1.0, step=0.2, value=0.0, marks={0: "Off", 0.6: "Modéré", 1.0: "Max"}),

                        html.Label("Immunothérapie (Inhibiteur PD-1 / PD-L1)", style={"color": "#10b981", "fontSize": "11.5px"}),
                        dcc.Slider(id="slider-ici", min=0.0, max=1.0, step=0.2, value=0.0, marks={0: "Off", 0.6: "Standard", 1.0: "Max"}),

                        html.Label("Perfusion IV Continue (mg/h)", style={"color": "#cbd5e1", "fontSize": "11.5px"}),
                        dcc.Slider(id="slider-infusion", min=0, max=40, step=5, value=0, marks={0: "0", 20: "20", 40: "40"}),
                    ]),

                    html.Div(id="oncology-fiche-box", style=CARD_STYLE)
                ]),

                # Colonne 2 : Imagerie Anatomique & Progression
                html.Div(style=CARD_STYLE, children=[
                    html.H3("Anatomie & Microenvironnement Tumoral", style={"marginTop": 0, "color": "#38bdf8", "fontSize": "14px", "textAlign": "center"}),
                    html.Iframe(id="svg-patient-frame", style={"width": "100%", "height": "300px", "border": "none"}),
                    html.Div(id="stage-indicator-box", style={"marginTop": "10px", "padding": "10px", "backgroundColor": "#0f172a", "borderRadius": "8px", "border": "1px solid #334155"})
                ]),

                # Colonne 3 : Courbes PK & Dynamique Clonale / Immunitaire
                html.Div(style=CARD_STYLE, children=[
                    dcc.Graph(id="graph-pk", style={"height": "225px"}),
                    dcc.Graph(id="graph-pd", style={"height": "275px"})
                ])
            ])
        ]),

        # ==================== ONGLET 2 : COHORTE VIRTUELLE & MACHINE LEARNING ====================
        dcc.Tab(label="2. Cohorte Virtuelle (VPop) & Découverte de Biomarqueurs par Machine Learning", value="tab-cohort",
                style={"color": "#94a3b8", "backgroundColor": "#1e293b", "padding": "10px"},
                selected_style={"color": "#38bdf8", "backgroundColor": "#0f172a", "fontWeight": "bold", "padding": "10px"},
                children=[
            html.Div(style={"display": "grid", "gridTemplateColumns": "315px 1fr", "gap": "15px", "marginTop": "15px"}, children=[

                # Panneau de Contrôle de l'Essai In Silico
                html.Div(style=CARD_STYLE, children=[
                    html.H3("Design de l'Essai & Trithérapie", style={"marginTop": 0, "color": "#38bdf8", "fontSize": "14px"}),

                    html.Label("Taille de la Cohorte (N patients virtuels)", style={"fontSize": "11.5px", "color": "#cbd5e1"}),
                    dcc.Slider(id="cohort-n", min=100, max=400, step=50, value=200, marks={100: "100", 250: "250", 400: "400"}),

                    html.Label("Dose Bolus IV (mg) / 72h", style={"fontSize": "11.5px", "color": "#cbd5e1"}),
                    dcc.Slider(id="cohort-dose", min=0, max=250, step=25, value=125, marks={0: "0", 125: "125", 250: "250"}),

                    html.Label("Co-traitement Anti-VEGF", style={"fontSize": "11.5px", "color": "#38bdf8"}),
                    dcc.Slider(id="cohort-vegf", min=0.0, max=1.0, step=0.2, value=0.4, marks={0: "Off", 0.4: "0.4", 1.0: "1.0"}),

                    html.Label("Co-traitement Immunothérapie (PD-1/PD-L1)", style={"fontSize": "11.5px", "color": "#10b981"}),
                    dcc.Slider(id="cohort-ici", min=0.0, max=1.0, step=0.2, value=0.6, marks={0: "Off", 0.6: "0.6", 1.0: "1.0"}),

                    html.Label("Durée de suivi clinique (jours)", style={"fontSize": "11.5px", "color": "#cbd5e1"}),
                    dcc.Slider(id="cohort-days", min=14, max=42, step=7, value=28, marks={14: "14j", 28: "28j", 42: "42j"}),

                    html.Button("Simuler Cohorte & Entraîner ML", id="btn-run-cohort", n_clicks=0, style={
                        "width": "100%", "backgroundColor": "#10b981", "color": "white", "border": "none",
                        "padding": "10px", "borderRadius": "6px", "cursor": "pointer", "fontWeight": "bold", "marginTop": "10px"
                    }),

                    html.Button("Exporter Cohorte (.CSV)", id="btn-export-csv", n_clicks=0, style={
                        "width": "100%", "backgroundColor": "#334155", "color": "#38bdf8", "border": "1px solid #38bdf8",
                        "padding": "8px", "borderRadius": "6px", "cursor": "pointer", "fontWeight": "bold", "marginTop": "8px"
                    }),
                    dcc.Download(id="download-cohort-csv"),

                    html.Div(id="cohort-kpi-box", style={"marginTop": "12px", "padding": "10px", "backgroundColor": "#0f172a", "borderRadius": "8px", "border": "1px solid #334155"})
                ]),

                # Grille des 4 Graphiques : Population, RECIST, ROC-AUC, Feature Importance
                html.Div(style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "14px"}, children=[
                    html.Div(style=CARD_STYLE, children=[
                        dcc.Graph(id="graph-cohort-pop", style={"height": "245px"})
                    ]),
                    html.Div(style=CARD_STYLE, children=[
                        dcc.Graph(id="graph-cohort-recist", style={"height": "245px"})
                    ]),
                    html.Div(style=CARD_STYLE, children=[
                        dcc.Graph(id="graph-ml-roc", style={"height": "255px"})
                    ]),
                    html.Div(style=CARD_STYLE, children=[
                        dcc.Graph(id="graph-ml-importance", style={"height": "255px"})
                    ])
                ])
            ])
        ])
    ]),

    dcc.Interval(id="timer", interval=1000, n_intervals=0, disabled=True)
])


# --- CALLBACKS ONGLET 1 ---
@app.callback(
    Output("timer", "disabled"),
    Output("btn-toggle", "children"),
    Output("btn-toggle", "style"),
    Input("btn-toggle", "n_clicks"),
    State("timer", "disabled"),
    prevent_initial_call=True
)
def toggle_sim(n_clicks, is_disabled):
    new_state = not is_disabled
    txt = "Pause" if not new_state else "Démarrer"
    col = "#f59e0b" if not new_state else "#10b981"
    return new_state, txt, {"backgroundColor": col, "color": "white", "border": "none", "padding": "8px", "borderRadius": "6px", "cursor": "pointer", "fontWeight": "bold", "flex": 1}


@app.callback(
    Output("svg-patient-frame", "srcDoc"),
    Output("oncology-fiche-box", "children"),
    Output("stage-indicator-box", "children"),
    Output("graph-pk", "figure"),
    Output("graph-pd", "figure"),
    Input("timer", "n_intervals"),
    Input("btn-bolus", "n_clicks"),
    Input("btn-reset", "n_clicks"),
    Input("dropdown-cancer", "value"),
    Input("slider-biomarker", "value"),
    Input("slider-vegf", "value"),
    Input("slider-ici", "value"),
    State("slider-res-frac", "value"),
    State("slider-infusion", "value"),
    State("timer", "disabled")
)
def update_individual_twin(n_intervals, bolus_clicks, reset_clicks, cancer_key,
                           biomarker_b0, anti_vegf, ici_val, res_pct, infusion_rate, is_disabled):
    global engine, logger
    ctx = dash.callback_context
    trigger = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else ""

    if trigger in ("btn-reset", "dropdown-cancer"):
        engine = PatientPKPDEngine(
            cancer_key=cancer_key,
            biomarker_b0=float(biomarker_b0 or 1.0),
            initial_res_frac=float(res_pct or 2) / 100.0,
            anti_vegf=float(anti_vegf or 0.0),
            immunotherapy=float(ici_val or 0.0)
        )
        logger = PatientLogger()
        logger.log(engine.state)
    elif trigger == "btn-bolus":
        engine.administer_bolus(120.0)
        logger.log(engine.state)
    elif trigger in ("slider-biomarker", "slider-vegf", "slider-ici"):
        engine.biomarker_b0 = float(biomarker_b0 or 1.0)
        engine.anti_vegf = float(anti_vegf or 0.0)
        engine.immunotherapy = float(ici_val or 0.0)
    elif not is_disabled:
        engine.biomarker_b0 = float(biomarker_b0 or 1.0)
        engine.anti_vegf = float(anti_vegf or 0.0)
        engine.immunotherapy = float(ici_val or 0.0)
        engine.step(dt=2.0, dose_rate=float(infusion_rate or 0))
        logger.log(engine.state)

    df = logger.to_dataframe()
    st = engine.state
    prof = engine.profile

    svg_doc = generate_patient_oncology_svg(prof, st.tumor_weight, st.conc_central, st.angiogenesis_index, st.progression_stage)

    fiche_children = [
        html.H4(f"Fiche Biologique — Rang #{prof.rank_globocan}", style={"margin": "0 0 6px 0", "color": "#38bdf8", "fontSize": "13.5px"}),
        html.P([html.B("Origine : "), prof.origin], style={"fontSize": "11.5px", "margin": "3px 0"}),
        html.P([html.B("Mutations : "), prof.mutations], style={"fontSize": "11.5px", "margin": "3px 0", "color": "#fca5a5"}),
        html.P([html.B("Biomarqueur : "), prof.biomarker_name], style={"fontSize": "11.5px", "margin": "3px 0", "color": "#38bdf8"}),
        html.P([html.B("Métastases : "), prof.metastatic_sites], style={"fontSize": "11.5px", "margin": "3px 0", "color": "#f59e0b"}),
    ]

    res_share = (st.tumor_resistant / max(1e-6, st.tumor_sensitive + st.tumor_resistant)) * 100.0
    stage_children = [
        html.Div(st.stage_label, style={"fontWeight": "bold", "fontSize": "12.5px", "color": "#f59e0b" if st.progression_stage >= 4 else "#10b981"}),
        html.Div(f"Fraction résistante xR : {res_share:.1f}%  |  Effecteur immunitaire E(t) : {st.immune_effector:.2f}", style={"fontSize": "11px", "color": "#94a3b8", "marginTop": "4px"})
    ]

    fig_pk = go.Figure()
    fig_pk.add_trace(go.Scatter(x=df["time"], y=df["conc_central"], mode="lines", name="C_centrale (mg/L)", line=dict(color="#38bdf8", width=2.5)))
    fig_pk.add_trace(go.Scatter(x=df["time"], y=df["conc_periph"], mode="lines", name="C_périphérique (mg/L)", line=dict(color="#a855f7", width=1.5, dash="dash")))
    fig_pk.update_layout(
        title="Pharmacocinétique (PK Bi-Compartimentale)",
        paper_bgcolor="#1e293b", plot_bgcolor="#0f172a", font=dict(color="#cbd5e1", size=11),
        margin=dict(l=40, r=15, t=30, b=26), xaxis=dict(title="Temps (h)", gridcolor="#334155"),
        yaxis=dict(title="mg/L", gridcolor="#334155"), uirevision="constant"
    )

    fig_pd = go.Figure()
    fig_pd.add_trace(go.Scatter(x=df["time"], y=df["tumor_weight"], mode="lines", name="Masse Totale w(t)", line=dict(color="#f8fafc", width=2.5)))
    fig_pd.add_trace(go.Scatter(x=df["time"], y=df["tumor_sensitive"], mode="lines", name="Clone Sensible xS", line=dict(color="#10b981", width=2)))
    fig_pd.add_trace(go.Scatter(x=df["time"], y=df["tumor_resistant"], mode="lines", name="Clone Résistant xR", line=dict(color="#ef4444", width=2)))
    fig_pd.add_trace(go.Scatter(x=df["time"], y=df["immune_effector"], mode="lines", name="Activité Immunitaire E(t)", line=dict(color="#38bdf8", width=1.5, dash="dot")))
    fig_pd.update_layout(
        title=f"Dynamique Clonale & Immunitaire — {prof.subtype}",
        paper_bgcolor="#1e293b", plot_bgcolor="#0f172a", font=dict(color="#cbd5e1", size=11),
        margin=dict(l=40, r=15, t=30, b=26), xaxis=dict(title="Temps (h)", gridcolor="#334155"),
        yaxis=dict(title="Masse (g) / Indice", gridcolor="#334155"), uirevision="constant"
    )

    return svg_doc, fiche_children, stage_children, fig_pk, fig_pd


# --- CALLBACKS ONGLET 2 : COHORTE VIRTUELLE & ML ---
@app.callback(
    Output("cohort-kpi-box", "children"),
    Output("graph-cohort-pop", "figure"),
    Output("graph-cohort-recist", "figure"),
    Output("graph-ml-roc", "figure"),
    Output("graph-ml-importance", "figure"),
    Input("btn-run-cohort", "n_clicks"),
    Input("dropdown-cancer", "value"),
    State("cohort-n", "value"),
    State("cohort-dose", "value"),
    State("cohort-vegf", "value"),
    State("cohort-ici", "value"),
    State("cohort-days", "value")
)
def update_virtual_cohort(n_clicks, cancer_key, n_patients, dose_mg, anti_vegf, ici_val, duration_days):
    global latest_cohort_df
    res = run_virtual_cohort_trial(
        cancer_key=cancer_key,
        n_patients=int(n_patients or 200),
        bolus_mg=float(dose_mg or 125),
        anti_vegf=float(anti_vegf or 0.4),
        immunotherapy=float(ici_val or 0.6),
        duration_days=int(duration_days or 28)
    )
    latest_cohort_df = res["cohort_df"]

    kpi_children = [
        html.Div(f"Taux de Réponse Objective (ORR) : {res['orr']:.1f} %", style={"fontWeight": "bold", "color": "#10b981", "fontSize": "13.5px"}),
        html.Div(f"Modèle ML (Random Forest) : ROC-AUC = {res['roc_auc']:.3f}", style={"fontWeight": "bold", "color": "#38bdf8", "fontSize": "12.5px", "marginTop": "5px"}),
        html.Div(res["ml_status"], style={"fontSize": "11px", "color": "#94a3b8", "marginTop": "4px"})
    ]

    # 1. Trajectoires Populationnelles
    fig_pop = go.Figure()
    t = res["time_days"]
    fig_pop.add_trace(go.Scatter(x=np.concatenate([t, t[::-1]]), y=np.concatenate([res["w_p95"], res["w_p05"][::-1]]),
                                 fill="toself", fillcolor="rgba(56, 189, 248, 0.15)", line=dict(color="rgba(255,255,255,0)"),
                                 name="IC 5% – 95%"))
    fig_pop.add_trace(go.Scatter(x=np.concatenate([t, t[::-1]]), y=np.concatenate([res["w_p75"], res["w_p25"][::-1]]),
                                 fill="toself", fillcolor="rgba(56, 189, 248, 0.30)", line=dict(color="rgba(255,255,255,0)"),
                                 name="IQR 25% – 75%"))
    fig_pop.add_trace(go.Scatter(x=t, y=res["w_median"], mode="lines", name="Médiane Cohorte",
                                 line=dict(color="#38bdf8", width=3)))
    fig_pop.update_layout(
        title=f"Dynamique Tumorale Populationnelle (N = {n_patients})",
        paper_bgcolor="#1e293b", plot_bgcolor="#0f172a", font=dict(color="#cbd5e1", size=11),
        margin=dict(l=40, r=15, t=32, b=28), xaxis=dict(title="Jours", gridcolor="#334155"),
        yaxis=dict(title="Masse Tumorale (g)", gridcolor="#334155")
    )

    # 2. Répartition RECIST 1.1
    categories = ["CR (Réponse Complète)", "PR (Réponse Partielle)", "SD (Maladie Stable)", "PD (Progression)"]
    colors = ["#10b981", "#38bdf8", "#f59e0b", "#ef4444"]
    counts = [int(np.sum(res["recist"] == c)) for c in categories]

    fig_recist = go.Figure(data=[go.Bar(x=["CR", "PR", "SD", "PD"], y=counts, marker_color=colors, text=counts, textposition="auto")])
    fig_recist.update_layout(
        title="Stratification Clinique RECIST 1.1",
        paper_bgcolor="#1e293b", plot_bgcolor="#0f172a", font=dict(color="#cbd5e1", size=11),
        margin=dict(l=35, r=15, t=32, b=28), yaxis=dict(title="Patients", gridcolor="#334155")
    )

    # 3. Courbe ROC (Validation Croisée)
    fig_roc = go.Figure()
    fig_roc.add_trace(go.Scatter(x=res["fpr"], y=res["tpr"], mode="lines",
                                 name=f"Random Forest (AUC = {res['roc_auc']:.3f})",
                                 line=dict(color="#10b981", width=3)))
    fig_roc.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", name="Aléatoire (AUC = 0.50)",
                                 line=dict(color="#64748b", width=1.5, dash="dash")))
    fig_roc.update_layout(
        title=f"Prédiction ML de la Réponse (Courbe ROC-AUC = {res['roc_auc']:.3f})",
        paper_bgcolor="#1e293b", plot_bgcolor="#0f172a", font=dict(color="#cbd5e1", size=11),
        margin=dict(l=40, r=15, t=32, b=28), xaxis=dict(title="Taux Faux Positifs (1 - Spécificité)", gridcolor="#334155"),
        yaxis=dict(title="Taux Vrais Positifs (Sensibilité)", gridcolor="#334155")
    )

    # 4. Importance des Biomarqueurs (Feature Importance)
    order = np.argsort(res["importances"])
    sorted_names = [res["feature_names"][i] for i in order]
    sorted_imp = res["importances"][order]

    fig_imp = go.Figure(data=[go.Bar(
        x=sorted_imp, y=sorted_names, orientation="h",
        marker_color="#38bdf8", text=[f"{v:.2f}" for v in sorted_imp], textposition="auto"
    )])
    fig_imp.update_layout(
        title="Importance des Biomarqueurs Prédictifs (Random Forest)",
        paper_bgcolor="#1e293b", plot_bgcolor="#0f172a", font=dict(color="#cbd5e1", size=11),
        margin=dict(l=175, r=20, t=32, b=28), xaxis=dict(title="Importance relative", gridcolor="#334155")
    )

    return kpi_children, fig_pop, fig_recist, fig_roc, fig_imp


@app.callback(
    Output("download-cohort-csv", "data"),
    Input("btn-export-csv", "n_clicks"),
    prevent_initial_call=True
)
def export_cohort_csv(n_clicks):
    if latest_cohort_df is not None:
        return dcc.send_data_frame(latest_cohort_df.to_csv, "virtual_patient_cohort_qsp.csv", index=False)
    return dash.no_update


if __name__ == "__main__":
    app.run(debug=False, port=8051)