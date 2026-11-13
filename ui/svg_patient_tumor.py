from core.cancer_catalog import CancerProfile

def generate_patient_oncology_svg(profile: CancerProfile, tumor_weight: float,
                                  conc_c: float, angio_index: float, stage: int) -> str:
    r_tumor = min(48, max(9, int(15 * (tumor_weight ** 0.5))))
    necrotic_r = max(0, r_tumor - 14)

    # Épaisseur de la barrière stromale (inverse de la perméabilité)
    stroma_thickness = max(1, int((1.0 - profile.stromal_permeability) * 10))
    r_stroma = r_tumor + stroma_thickness

    blood_alpha = min(1.0, conc_c / 18.0)
    drug_color = f"rgba(56, 189, 248, {blood_alpha:.2f})"

    ox, oy = profile.organ_coords

    # Génération des métastases si Stade >= 5 (Intravasation / Dissémination / Métastases)
    metastasis_svg = ""
    if stage >= 5:
        active_mets = profile.metastasis_coords[:1] if stage == 5 else (
            profile.metastasis_coords[:2] if stage == 6 else profile.metastasis_coords
        )
        for mx, my in active_mets:
            metastasis_svg += f"""
            <line x1="{ox}" y1="{oy}" x2="{mx}" y2="{my}" stroke="#f59e0b" stroke-width="1" stroke-dasharray="2,2"/>
            <circle cx="{mx}" cy="{my}" r="3.2" fill="#f59e0b">
                <animate attributeName="opacity" values="0.4;1;0.4" dur="1.4s" repeatCount="indefinite"/>
            </circle>
            """

    # Spicules d'invasion EMT / MMP si Stade >= 4
    emt_spikes = ""
    if stage >= 4:
        for angle in range(0, 360, 45):
            emt_spikes += f'<line x1="0" y1="0" x2="0" y2="{-(r_stroma + 7)}" stroke="#f87171" stroke-width="1.5" transform="rotate({angle})"/>'

    return f"""<!DOCTYPE html>
<html>
<head>
<style>
    body {{ margin: 0; padding: 0; background-color: transparent; overflow: hidden; font-family: 'Segoe UI', sans-serif; }}
    @keyframes pulseVessel {{
        0% {{ stroke-dashoffset: 0; }}
        100% {{ stroke-dashoffset: -20; }}
    }}
    .blood-flow {{
        stroke-dasharray: 6, 4;
        animation: pulseVessel 1.4s linear infinite;
    }}
</style>
</head>
<body>
<svg width="100%" height="295" viewBox="0 0 330 300" xmlns="http://www.w3.org/2000/svg">
    <defs>
        <radialGradient id="necroticCore" cx="50%" cy="50%" r="50%">
            <stop offset="0%" stop-color="#0f172a"/>
            <stop offset="75%" stop-color="#4c1d95"/>
            <stop offset="100%" stop-color="#991b1b"/>
        </radialGradient>
        <radialGradient id="prolifRim" cx="50%" cy="50%" r="50%">
            <stop offset="55%" stop-color="#dc2626" stop-opacity="0.85"/>
            <stop offset="100%" stop-color="#f87171" stop-opacity="0.95"/>
        </radialGradient>
    </defs>

    <!-- Silhouette Anatomique Multi-Organes (Gauche) -->
    <g transform="translate(8, 8)">
        <circle cx="50" cy="22" r="13" fill="#1e293b" stroke="#64748b" stroke-width="1.5"/>
        <rect x="44" y="35" width="12" height="12" fill="#1e293b"/>
        <path d="M 26 48 Q 50 42 74 48 L 78 150 L 22 150 Z" fill="#1e293b" stroke="#64748b" stroke-width="1.5"/>

        <!-- Repères anatomiques subtils (Poumons, Foie, Tube digestif) -->
        <ellipse cx="38" cy="80" rx="10" ry="16" fill="#334155" opacity="0.45"/>
        <ellipse cx="62" cy="80" rx="10" ry="16" fill="#334155" opacity="0.45"/>
        <path d="M 35 98 Q 50 95 62 104 L 38 108 Z" fill="#475569" opacity="0.5"/>

        <!-- Réseau Vasculaire -->
        <path d="M 50 38 L 50 135 M 50 75 L 36 85 M 50 75 L 64 85" fill="none" stroke="#ef4444" stroke-width="1.6" opacity="0.7"/>
        <path class="blood-flow" d="M 50 38 L 50 135" fill="none" stroke="{drug_color}" stroke-width="2.5"/>

        <!-- Foyer Tumoral Primaire Dynamique -->
        <circle cx="{ox}" cy="{oy}" r="6" fill="#ef4444">
            <animate attributeName="r" values="4.5;7.5;4.5" dur="1.8s" repeatCount="indefinite"/>
        </circle>
        {metastasis_svg}

        <text x="50" y="172" fill="#cbd5e1" font-size="10" font-weight="bold" text-anchor="middle">Localisation Cible</text>
        <text x="50" y="186" fill="#f59e0b" font-size="9" text-anchor="middle">{'Métastases actives' if stage >= 5 else 'Lésion localisée'}</text>
    </g>

    <!-- Connexion Zoom Microenvironnement -->
    <line x1="{ox + 12}" y1="{oy + 8}" x2="155" y2="135" stroke="#38bdf8" stroke-width="1" stroke-dasharray="3,3" opacity="0.7"/>

    <!-- Microenvironnement Tumoral & Stroma (Droite) -->
    <g transform="translate(155, 14)">
        <rect x="0" y="0" width="165" height="265" rx="10" fill="#0f172a" stroke="#334155" stroke-width="1.5"/>
        <text x="82" y="19" fill="#38bdf8" font-size="11" font-weight="bold" text-anchor="middle">Microenvironnement</text>

        <!-- Vaisseaux Angiogéniques (VEGF) -->
        <path d="M 12 245 Q 42 190 16 135 T 42 35" fill="none" stroke="#b91c1c" stroke-width="{max(1, int(angio_index * 1.8))}"/>
        <path d="M 153 245 Q 123 185 149 125 T 123 35" fill="none" stroke="#b91c1c" stroke-width="{max(1, int(angio_index * 1.5))}"/>
        <path class="blood-flow" d="M 12 245 Q 42 190 16 135 T 42 35" fill="none" stroke="#38bdf8" stroke-width="2"/>

        <!-- Tumeur Multi-Couches -->
        <g transform="translate(82, 132)">
            {emt_spikes}
            <!-- Barrière Stromale Desmoplastique -->
            <circle cx="0" cy="0" r="{r_stroma}" fill="none" stroke="#94a3b8" stroke-width="{stroma_thickness}" stroke-dasharray="3,2" opacity="0.8"/>
            <!-- Couronne Proliférante x1 -->
            <circle cx="0" cy="0" r="{r_tumor}" fill="url(#prolifRim)" stroke="#fca5a5" stroke-width="1.2"/>
            <!-- Cœur Nécrotique / Hypoxique -->
            {f'<circle cx="0" cy="0" r="{necrotic_r}" fill="url(#necroticCore)"/>' if necrotic_r > 0 else ''}
        </g>

        <text x="82" y="218" fill="#cbd5e1" font-size="10" text-anchor="middle">Masse : {tumor_weight:.2f} g  |  Cc : {conc_c:.1f} mg/L</text>
        <text x="82" y="234" fill="#94a3b8" font-size="9" text-anchor="middle">Perméabilité stromale : {int(profile.stromal_permeability * 100)}%</text>
        <text x="82" y="250" fill="#f59e0b" font-size="9" font-weight="bold" text-anchor="middle">Stade {stage}/7</text>
    </g>
</svg>
</body>
</html>"""