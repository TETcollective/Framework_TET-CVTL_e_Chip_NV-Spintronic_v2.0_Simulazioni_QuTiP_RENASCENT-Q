"""
Termine Retrocausale-Negentropico Custom (D_beta)
TET Collective — RENASCENT-Q
"""

import qutip as qt

def custom_d_beta_superoperator(beta=0.382, alpha=0.08):
    """Costruisce il superoperatore D_beta[ρ]"""
    # Stato Bell |Φ⁺⟩
    Phi_plus = (qt.basis(4, 0) + qt.basis(4, 3)).unit()
    proj_bell = Phi_plus * Phi_plus.dag()
    
    XX = qt.tensor(qt.sigmax(), qt.sigmax())
    O_future = proj_bell + alpha * XX
    
    # Costruzione dissipatore
    term1 = qt.spre(O_future) * qt.spost(O_future.dag())
    term2 = 0.5 * (qt.spre(O_future.dag() * O_future) + qt.spost(O_future.dag() * O_future))
    
    D = beta * (term1 - term2)
    D = D + D.dag()   # hermitiano coniugato
    
    return D