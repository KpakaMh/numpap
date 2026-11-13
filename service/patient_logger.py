import pandas as pd

class PatientLogger:
    def __init__(self):
        self.history = []

    def log(self, state):
        self.history.append({
            "time": state.time,
            "conc_central": state.conc_central,
            "conc_periph": state.conc_periph,
            "tumor_sensitive": state.tumor_sensitive,
            "tumor_resistant": state.tumor_resistant,
            "tumor_transit": state.tumor_transit1 + state.tumor_transit2 + state.tumor_transit3,
            "tumor_weight": state.tumor_weight,
            "immune_effector": state.immune_effector,
            "angiogenesis": state.angiogenesis_index
        })

    def to_dataframe(self):
        cols = [
            "time", "conc_central", "conc_periph", "tumor_sensitive",
            "tumor_resistant", "tumor_transit", "tumor_weight",
            "immune_effector", "angiogenesis"
        ]
        if not self.history:
            return pd.DataFrame(columns=cols)
        return pd.DataFrame(self.history)