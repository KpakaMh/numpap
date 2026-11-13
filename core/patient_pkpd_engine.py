import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass

@dataclass
class PatientState:
    time: float = 0.0
    conc_central: float = 0.0      # C_c (mg/L)
    conc_periph: float = 0.0       # C_p (mg/L)
    tumor_prolif: float = 0.2      # x1 (g ou cm3)
    tumor_transit1: float = 0.0    # x2
    tumor_transit2: float = 0.0    # x3
    tumor_transit3: float = 0.0    # x4
    tumor_weight: float = 0.2      # w total (g)
    angiogenesis_index: float = 1.0 # Facteur de vascularisation

class PatientPKPDEngine:
    def __init__(self, 
                 vc: float = 5.0,        # Volume central (L)
                 vp: float = 15.0,       # Volume périphérique (L)
                 k10: float = 0.15,      # Élimination (1/h)
                 k12: float = 0.08,      # Inter-compartiment (1/h)
                 k21: float = 0.04,      # Inter-compartiment (1/h)
                 lambda0: float = 0.012, # Taux de croissance exponentielle (1/h)
                 lambda1: float = 0.05,  # Taux de croissance linéaire (g/h)
                 psi: float = 20.0,      # Constante de transition Simeoni
                 k1: float = 0.02,       # Constante de transit de mort cellulaire (1/h)
                 k2: float = 0.005,      # Puissance du médicament (L/mg/h)
                 initial_tumor_w: float = 0.2):
        
        self.vc = vc
        self.vp = vp
        self.k10 = k10
        self.k12 = k12
        self.k21 = k21
        self.lambda0 = lambda0
        self.lambda1 = lambda1
        self.psi = psi
        self.k1 = k1
        self.k2 = k2

        self.state = PatientState(tumor_prolif=initial_tumor_w, tumor_weight=initial_tumor_w)

    def administer_bolus(self, dose_mg: float):
        """Ajoute une dose bolus IV dans le compartiment central."""
        self.state.conc_central += dose_mg / self.vc

    def _ode_system(self, t, y, dose_rate):
        Cc, Cp, x1, x2, x3, x4 = y

        # PK
        dCc = -(self.k10 + self.k12) * Cc + self.k21 * Cp + (dose_rate / self.vc)
        dCp = self.k12 * Cc - self.k21 * Cp

        # PD (Simeoni et al.)
        w = max(1e-6, x1 + x2 + x3 + x4)
        growth_rate = (self.lambda0 * x1) / ((1.0 + (self.lambda0 / self.lambda1 * w) ** self.psi) ** (1.0 / self.psi))
        
        kill_effect = self.k2 * Cc * x1

        dx1 = growth_rate - kill_effect
        dx2 = kill_effect - self.k1 * x2
        dx3 = self.k1 * (x2 - x3)
        dx4 = self.k1 * (x3 - x4)

        return [dCc, dCp, dx1, dx2, dx3, dx4]

    def step(self, dt: float = 1.0, dose_rate: float = 0.0):
        """Avance la simulation de dt heures."""
        y0 = [
            self.state.conc_central,
            self.state.conc_periph,
            self.state.tumor_prolif,
            self.state.tumor_transit1,
            self.state.tumor_transit2,
            self.state.tumor_transit3
        ]

        sol = solve_ivp(
            self._ode_system,
            [0, dt],
            y0,
            args=(dose_rate,),
            method='RK45'
        )

        y_next = sol.y[:, -1]
        
        self.state.time += dt
        self.state.conc_central = max(0.0, y_next[0])
        self.state.conc_periph = max(0.0, y_next[1])
        self.state.tumor_prolif = max(0.0, y_next[2])
        self.state.tumor_transit1 = max(0.0, y_next[3])
        self.state.tumor_transit2 = max(0.0, y_next[4])
        self.state.tumor_transit3 = max(0.0, y_next[5])
        
        self.state.tumor_weight = (
            self.state.tumor_prolif + 
            self.state.tumor_transit1 + 
            self.state.tumor_transit2 + 
            self.state.tumor_transit3
        )
        
        # Angiogenèse stimulée par la masse tumorale
        self.state.angiogenesis_index = min(3.0, 0.5 + 1.5 * (self.state.tumor_weight / 5.0))
        return self.state