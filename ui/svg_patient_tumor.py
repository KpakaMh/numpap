def generate_patient_oncology_svg(tumor_weight: float, conc_c: float, angio_index: float) -> str:
    # Calcul du rayon de la tumeur (normalisé)
    r_tumor = min(55, max(10, int(15 * (tumor_weight ** 0.5))))
    necrotic_r = max(0, r_tumor - 15)
    
    # Couleur du sang selon la concentration de drogue
    blood_alpha = min(1.0, conc_c / 20.0)
    drug_color = f"rgba(56, 189, 248, {blood_alpha:.2f})"
    
    html_content = f"""<!DOCTYPE html>
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
        animation: pulseVessel 1.5s linear infinite;
    }}
</style>
</head>
<body>
<svg width="100%" height="290" viewBox="0 0 320 300" xmlns="http://www.w3.org/2000/svg">
    <defs>
        <radialGradient id="necroticCore" cx="50%" cy="50%" r="50%">
            <stop offset="0%" stop-color="#1e1b4b"/>
            <stop offset="70%" stop-color="#4c1d95"/>
            <stop offset="100%" stop-color="#991b1b"/>
        </radialGradient>
        <radialGradient id="prolifRim" cx="50%" cy="50%" r="50%">
            <stop offset="60%" stop-color="#dc2626" stop-opacity="0.8"/>
            <stop offset="100%" stop-color="#f87171" stop-opacity="0.95"/>
        </radialGradient>
    </defs>

    <!-- Silhouette Humaine Stylisée (Gauche) -->
    <g transform="translate(10, 10)" opacity="0.85">
        <!-- Tête et Buste -->
        <circle cx="50" cy="25" r="14" fill="#334155" stroke="#94a3b8" stroke-width="1.5"/>
        <path d="M 30 45 Q 50 40 70 45 L 75 140 L 25 140 Z" fill="#1e293b" stroke="#64748b" stroke-width="1.5"/>
        
        <!-- Système Vasculaire Cible -->
        <path d="M 50 45 L 50 120 M 50 65 L 35 85 M 50 65 L 65 85" fill="none" stroke="#ef4444" stroke-width="2"/>
        <path class="blood-flow" d="M 50 45 L 50 120" fill="none" stroke="{drug_color}" stroke-width="2.5"/>

        <!-- Organe Cible (Poumon / Tumeur Primaire) -->
        <circle cx="65" cy="85" r="10" fill="#475569" opacity="0.6"/>
        <circle cx="65" cy="85" r="5" fill="#ef4444">
            <animate attributeName="r" values="4;6;4" dur="2s" repeatCount="indefinite"/>
        </circle>
        <text x="50" y="160" fill="#94a3b8" font-size="10" text-anchor="middle">Patient In Vivo</text>
    </g>

    <!-- Ligne de Connexion Microenvironnement -->
    <line x1="80" y1="90" x2="160" y2="135" stroke="#38bdf8" stroke-width="1" stroke-dasharray="3,3" opacity="0.6"/>

    <!-- Microenvironnement Tumoral (Droite) -->
    <g transform="translate(160, 20)">
        <!-- Cadre d'observation -->
        <rect x="0" y="0" width="150" height="250" rx="10" fill="#0f172a" stroke="#334155" stroke-width="1.5"/>
        <text x="75" y="20" fill="#38bdf8" font-size="11" font-weight="bold" text-anchor="middle">Microenvironnement</text>

        <!-- Vaisseaux Angiogéniques -->
        <path d="M 10 230 Q 40 180 15 130 T 40 40" fill="none" stroke="#b91c1c" stroke-width="{int(angio_index * 2)}"/>
        <path d="M 140 230 Q 110 170 135 110 T 110 40" fill="none" stroke="#b91c1c" stroke-width="{int(angio_index * 1.5)}"/>
        
        <!-- Perfusion Médicamenteuse dans les Vaisseaux -->
        <path class="blood-flow" d="M 10 230 Q 40 180 15 130 T 40 40" fill="none" stroke="#38bdf8" stroke-width="2"/>

        <!-- Tumeur 3D Multi-couches -->
        <g transform="translate(75, 135)">
            <!-- Zone Proliférante Externe -->
            <circle cx="0" cy="0" r="{r_tumor}" fill="url(#prolifRim)" stroke="#fca5a5" stroke-width="1.5"/>
            <!-- Cœur Hypoxique / Nécrotique -->
            {f'<circle cx="0" cy="0" r="{necrotic_r}" fill="url(#necroticCore)" stroke="#6b21a8" stroke-width="1"/>' if necrotic_r > 0 else ''}
        </g>

        <!-- Informations Textuelles -->
        <text x="75" y="215" fill="#cbd5e1" font-size="10" text-anchor="middle">Masse : {tumor_weight:.2f} g</text>
        <text x="75" y="230" fill="#38bdf8" font-size="10" text-anchor="middle">C_Cc : {conc_c:.1f} mg/L</text>
    </g>
</svg>
</body>
</html>"""
    return html_content