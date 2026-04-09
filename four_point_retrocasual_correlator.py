"""
RENASCENT-Q: Simulazione Correlatore Four-Point Retrocausale
Autore: PhysSoliman - TET Collective
Descrizione: Evoluzione del correlatore a quattro punti con e senza kick retrocausale-negentropico.
             Utile per la sezione risultati numerici del paper.
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

# ==================== Parametri ====================
beta = 0.382          # φ^{-2}
gamma = 0.07          # tasso decoerenza standard
tlist = np.linspace(0, 30, 500)

# ==================== Operatori su 4 qubit ====================
sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2), qt.qeye(2))
sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2), qt.qeye(2))
sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz(), qt.qeye(2))
sz4 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2), qt.sigmaz())

O1, O2, O3, O4 = sz1, sz2, sz3, sz4

# Stato iniziale fortemente correlato (GHZ-like generalizzato)
psi0 = (qt.basis(16, 0) + qt.basis(16, 15)).unit()

# Hamiltoniano di evoluzione debole
H = 0.4 * (O1*O2 + O2*O3 + O3*O4)

# Collassatori decoerenza standard
c_ops_std = [
    np.sqrt(gamma) * qt.tensor(qt.destroy(2), qt.qeye(2), qt.qeye(2), qt.qeye(2)),
    np.sqrt(gamma) * qt.tensor(qt.qeye(2), qt.destroy(2), qt.qeye(2), qt.qeye(2)),
    np.sqrt(gamma) * qt.tensor(qt.qeye(2), qt.qeye(2), qt.destroy(2), qt.qeye(2)),
    np.sqrt(gamma) * qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2), qt.destroy(2))
]

# Termine retrocausale-negentropico
O_neg = (O1 + O2 + O3 + O4) / 2.0
c_ops_retro = [np.sqrt(beta) * O_neg]

# ==================== Evoluzione ====================
result_std = qt.mesolve(H, psi0, tlist, c_ops_std, e_ops=[O1*O2*O3*O4])
result_retro = qt.mesolve(H, psi0, tlist, c_ops_std + c_ops_retro, e_ops=[O1*O2*O3*O4])

# ==================== Plot ====================
plt.figure(figsize=(11, 7))

plt.plot(tlist, np.abs(result_std.expect[0]), 'r--', linewidth=2.5, 
         label='Four-point standard (β = 0)')
plt.plot(tlist, np.abs(result_retro.expect[0]), 'b-', linewidth=3.2, 
         label=f'Four-point retrocausale (β = {beta})')

plt.xlabel('Tempo')
plt.ylabel(r'$|C^{(4)}(t_1,t_2,t_3,t_4)|$')
plt.title('Persistenza del correlatore four-point retrocausale nel framework RENASCENT-Q')
plt.legend()
plt.grid(True, alpha=0.35)
plt.ylim(0, 1.05)
plt.tight_layout()

# Salva la figura per Overleaf
plt.savefig('four_point_retrocausal_correlator.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Statistiche finali
print(f"Valore finale correlatore four-point (standard): {np.abs(result_std.expect[0][-1]):.4f}")
print(f"Valore finale correlatore four-point (retro)   : {np.abs(result_retro.expect[0][-1]):.4f}")
print(f"Guadagno di persistenza grazie a β: {np.abs(result_retro.expect[0][-1]) / np.abs(result_std.expect[0][-1]):.2f}x")