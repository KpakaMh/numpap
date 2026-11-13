import pandas as pd

class PatientLogger:
    def __init__(self):
        self.history = []

    def log(self, state):
        self.history.append({
            "time": state.time,
            "conc_central": state.conc_central,
            "conc_periph": state.conc_periph,
            "tumor_prolif": state.tumor_prolif,
            "tumor_transit": state.tumor_transit1 + state.tumor_transit2 + state.tumor_transit3,
            "tumor_weight": state.tumor_weight,
            "angiogenesis": state.angiogenesis_index
        })

    def to_dataframe(self):
        if not self.history:
            return pd.DataFrame(columns=["time", "conc_central", "conc_periph", "tumor_prolif", "tumor_transit", "tumor_weight", "angiogenesis"])
        return pd.DataFrame(self.history)