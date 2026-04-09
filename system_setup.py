"""
System Setup for RENASCENT-Q Proxy
TET Collective — Due qubit (dimer microtubuli / NV-center proxy)
Autore: PhysSoliman / TET Collective
"""

import qutip as qt
import numpy as np

def create_two_qubit_operators():
    """Restituisce tutti gli operatori necessari per il sistema a due qubit."""
    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2))
    sx2 = qt.tensor(qt.qeye(2), qt.sigmax())
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz())
    sm1 = qt.tensor(qt.sigmam(), qt.qeye(2))
    sm2 = qt.tensor(qt.qeye(2), qt.sigmam())
    return sx1, sx2, sz1, sz2, sm1, sm2


def get_system_parameters(version="standard"):
    """Parametri del sistema - due versioni disponibili"""
    if version == "standard":
        return {
            "Delta": 1.0,      # splitting energetico (unità arbitrarie)
            "J": 0.5,          # accoppiamento Ising-like
            "g": 0.35,         # ampiezza drive SAW
            "omega": 2.0 * np.pi,  # frequenza SAW
            "gamma_relax": 0.08,   # tasso rilassamento
            "gamma_dephase": 0.04, # tasso dephasing
            "T": 300,          # temperatura (K)
            "beta": ((1 + np.sqrt(5))/2)**(-2),  # φ^{-2} ≈ 0.382
            "alpha": 0.08      # parametro weak-value per O_future
        }
    elif version == "modified":
        return {
            "Delta": 1.2,
            "J": 0.65,
            "g": 0.45,         # drive più forte
            "omega": 2.5 * np.pi,
            "gamma_relax": 0.12,
            "gamma_dephase": 0.06,
            "T": 300,
            "beta": ((1 + np.sqrt(5))/2)**(-2),
            "alpha": 0.12      # weak-value più marcato
        }