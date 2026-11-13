import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objs as go

from core.patient_pkpd_engine import PatientPKPDEngine
from services.patient_logger import PatientLogger
from ui.svg_patient_tumor import generate_patient_oncology_svg

engine = PatientPKPDEngine()
logger = PatientLogger()
logger.log(engine.state)

app = dash.Dash(__name__, title="Patient Digital Twin - PK/PD Oncology")

app.layout = html.Div(style={'backgroundColor': '#0b0f19', 'color': '#f8fafc', 'minHeight': '100vh', 'padding': '20px', 'fontFamily': 'Segoe UI, sans-serif'}, children=[
    
    # En-tête
    html.Div(style={'display': 'flex', 'justifyContent': 'space-between', 'alignItems': 'center', 'borderBottom': '1px solid #1e293b', 'paddingBottom': '15px'}, children=[
        html.Div([
            html.H1("Jumeau Virtuel Patient (PK/PD & Oncologie)", style={'margin': '0', 'fontSize': '24px', 'color': '#38bdf8'}),
            html.P("Simulateur In Silico de Réponse Tumorale & Cinetique Medicamenteuse (Modèle Simeoni / CRCL)", style={'margin': '5px 0 0 0', 'color': '#94a3b8', 'fontSize': '14px'})
        ]),
        html.Div([
            html.Button("Injecter Bolus (100mg)", id="btn-bolus", n_clicks=0, style={'backgroundColor': '#0284c7', 'color': 'white', 'border': 'none', 'padding': '10px 16px', 'borderRadius': '6px', 'cursor': 'pointer', 'fontWeight': 'bold', 'marginRight': '10px'}),
            html.Button("Démarrer", id="btn-toggle", n_clicks=0, style={'backgroundColor': '#10b981', 'color': 'white', 'border': 'none', 'padding': '10px 16px', 'borderRadius': '6px', 'cursor': 'pointer', 'fontWeight': 'bold', 'marginRight': '10px'}),
            html.Button("Réinitialiser", id="btn-reset", n_clicks=0, style={'backgroundColor': '#ef4444', 'color': 'white', 'border': 'none', 'padding': '10px 16px', 'borderRadius': '6px', 'cursor': 'pointer', 'fontWeight': 'bold'})
        ])
    ]),

    # Grille Principale
    html.Div(style={'display': 'grid', 'gridTemplateColumns': '320px 340px 1fr', 'gap': '20px', 'marginTop': '20px'}, children=[
        
        # Panneau Paramètres Cliniques
        html.Div(style={'backgroundColor': '#1e293b', 'padding': '20px', 'borderRadius': '10px', 'border': '1px solid #334155'}, children=[
            html.H3("Paramètres Posologiques & Patient", style={'marginTop': '0', 'color': '#38bdf8', 'fontSize': '16px'}),
            
            html.Label("Perfusion Continue (mg/h)", style={'color': '#cbd5e1', 'fontSize': '13px'}),
            dcc.Slider(id="slider-infusion", min=0, max=50, step=5, value=0, marks={0: '0', 25: '25', 50: '50'}, tooltip={"placement": "bottom", "always_visible": True}),
            html.Br(),
            
            html.Label("Puissance Médicament k2 (Potency)", style={'color': '#cbd5e1', 'fontSize': '13px'}),
            dcc.Slider(id="slider-k2", min=0.001, max=0.02, step=0.002, value=0.005, marks={0.001: 'Faible', 0.02: 'Forte'}, tooltip={"placement": "bottom", "always_visible": True}),
            html.Br(),
            
            html.Label("Clairance Patient k10 (1/h)", style={'color': '#cbd5e1', 'fontSize': '13px'}),
            dcc.Slider(id="slider-k10", min=0.05, max=0.4, step=0.05, value=0.15, marks={0.05: 'Lente', 0.4: 'Rapide'}, tooltip={"placement": "bottom", "always_visible": True}),
        ]),

        # Rendu Anatomie & Microenvironnement SVG
        html.Div(style={'backgroundColor': '#1e293b', 'padding': '15px', 'borderRadius': '10px', 'border': '1px solid #334155', 'textAlign': 'center'}, children=[
            html.H3("Imagerie Virtuelle & Tumeur", style={'marginTop': '0', 'color': '#38bdf8', 'fontSize': '16px'}),
            html.Iframe(id="svg-patient-frame", style={'width': '100%', 'height': '295px', 'border': 'none'})
        ]),

        # Graphiques Temps Réel
        html.Div(style={'backgroundColor': '#1e293b', 'padding': '20px', 'borderRadius': '10px', 'border': '1px solid #334155'}, children=[
            dcc.Graph(id="graph-pk", style={'height': '240px'}),
            dcc.Graph(id="graph-pd", style={'height': '240px'})
        ])
    ]),

    dcc.Interval(id="timer", interval=1000, n_intervals=0, disabled=True)
])


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
    btn_txt = "Pause" if not new_state else "Démarrer"
    btn_col = "#f59e0b" if not new_state else "#10b981"
    return new_state, btn_txt, {'backgroundColor': btn_col, 'color': 'white', 'border': 'none', 'padding': '10px 16px', 'borderRadius': '6px', 'cursor': 'pointer', 'fontWeight': 'bold', 'marginRight': '10px'}


@app.callback(
    Output("svg-patient-frame", "srcDoc"),
    Output("graph-pk", "figure"),
    Output("graph-pd", "figure"),
    Input("timer", "n_intervals"),
    Input("btn-bolus", "n_clicks"),
    Input("btn-reset", "n_clicks"),
    State("slider-infusion", "value"),
    State("slider-k2", "value"),
    State("slider-k10", "value"),
    State("timer", "disabled")
)
def update_patient(n_intervals, bolus_clicks, reset_clicks, infusion_rate, k2, k10, is_disabled):
    global engine, logger
    
    ctx = dash.callback_context
    trigger = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else ""

    if trigger == "btn-reset":
        engine = PatientPKPDEngine()
        logger = PatientLogger()
        logger.log(engine.state)
    elif trigger == "btn-bolus":
        engine.administer_bolus(100.0)
        logger.log(engine.state)
    elif not is_disabled:
        engine.k2 = k2 or 0.005
        engine.k10 = k10 or 0.15
        engine.step(dt=2.0, dose_rate=float(infusion_rate or 0))
        logger.log(engine.state)

    df = logger.to_dataframe()
    st = engine.state

    svg_doc = generate_patient_oncology_svg(st.tumor_weight, st.conc_central, st.angiogenesis_index)

    # Graph PK
    fig_pk = go.Figure()
    fig_pk.add_trace(go.Scatter(x=df["time"], y=df["conc_central"], mode='lines', name='C_centrale (mg/L)', line=dict(color='#38bdf8', width=2.5)))
    fig_pk.add_trace(go.Scatter(x=df["time"], y=df["conc_periph"], mode='lines', name='C_périphérique (mg/L)', line=dict(color='#a855f7', width=1.5, dash='dash')))
    fig_pk.update_layout(
        title="Profil Pharmacocinétique (PK - 2 Compartiments)",
        paper_bgcolor='#1e293b', plot_bgcolor='#0f172a', font=dict(color='#cbd5e1'),
        margin=dict(l=40, r=20, t=35, b=30), xaxis=dict(title="Temps (h)", gridcolor='#334155'),
        yaxis=dict(title="Concentration (mg/L)", gridcolor='#334155'), uirevision="constant"
    )

    # Graph PD
    fig_pd = go.Figure()
    fig_pd.add_trace(go.Scatter(x=df["time"], y=df["tumor_weight"], mode='lines', name='Masse Tumorale Totale (g)', line=dict(color='#ef4444', width=3)))
    fig_pd.add_trace(go.Scatter(x=df["time"], y=df["tumor_prolif"], mode='lines', name='Cellules Proliférantes x1', line=dict(color='#10b981', width=1.5)))
    fig_pd.add_trace(go.Scatter(x=df["time"], y=df["tumor_transit"], mode='lines', name='Cellules en Transit Nécrotique', line=dict(color='#f59e0b', width=1.5, dash='dot')))
    fig_pd.update_layout(
        title="Réponse Tumorale Dynamique (PD - Simeoni Transit Model)",
        paper_bgcolor='#1e293b', plot_bgcolor='#0f172a', font=dict(color='#cbd5e1'),
        margin=dict(l=40, r=20, t=35, b=30), xaxis=dict(title="Temps (h)", gridcolor='#334155'),
        yaxis=dict(title="Poids Tumeur (g)", gridcolor='#334155'), uirevision="constant"
    )

    return svg_doc, fig_pk, fig_pd


if __name__ == "__main__":
    app.run(debug=False, port=8051)