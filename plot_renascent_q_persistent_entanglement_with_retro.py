
"""
RENASCENT-Q: Simulazione QuTiP con termine retrocausale-negentropico custom
Proxy due qubit con modulazione SAW e decoerenza termica.
Versione corretta - 03 Aprile 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# ====================== Parametri fisici ======================
omega = 1.0                # frequenza qubit base (arbitrary units)
g = 0.5                    # accoppiamento strain/SAW
gamma_diss = 0.08          # rate dissipazione termica
gamma_neg = 0.18           # rate negentropico base (da ottimizzare)
phi = (1 + np.sqrt(5)) / 2
beta = 1.0 / phi**2        # ≈ 0.382   phi^{-2}

T = 2.0                    # tempo massimo simulazione
N_steps = 400
times = np.linspace(0, T, N_steps)

# ====================== Operatori ======================
# Due qubit
sz1 = tensor(sigmaz(), qeye(2))
sz2 = tensor(qeye(2), sigmaz())
sx1 = tensor(sigmax(), qeye(2))
sx2 = tensor(qeye(2), sigmax())
sy1 = tensor(sigmay(), qeye(2))
sy2 = tensor(qeye(2), sigmay())
sm1 = tensor(sigmam(), qeye(2))
sm2 = tensor(qeye(2), sigmam())

# Hamiltoniano
H0 = omega * (sz1 + sz2) / 2
H_int = g * (sx1*sx2 + sy1*sy2)          # accoppiamento trasverso favorevole a stati Bell-like
H = H0 + H_int

# Stato iniziale entangled (Bell-like)
psi0 = (tensor(basis(2,0), basis(2,0)) + tensor(basis(2,1), basis(2,1))).unit()
rho0 = ket2dm(psi0)

# Dissipatori standard (decoerenza termica)
c_ops = [np.sqrt(gamma_diss) * sm1, np.sqrt(gamma_diss) * sm2]

# ====================== Evoluzione baseline (β=0) ======================
result_base = mesolve(H, rho0, times, c_ops,
                      e_ops=[lambda t, state: entropy_vn(state), lambda t, state: concurrence(state)])

# ====================== Termine retrocausale-negentropico ======================
# Approssimazione Markoviana del termine custom (L_neg ≈ raising operator)
c_neg = np.sqrt(gamma_neg * beta) * (sm1.dag() + sm2.dag())
c_ops_retro = c_ops + [c_neg]

result_retro = mesolve(H, rho0, times, c_ops_retro,
                       e_ops=[lambda t, state: entropy_vn(state), lambda t, state: concurrence(state)])

# ====================== Plot ======================
fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

axs[0].plot(times, result_base.expect[0], label='β = 0 (baseline)', linewidth=2.2)
axs[0].plot(times, result_retro.expect[0], label=f'β = φ⁻² ≈ {beta:.3f}', linewidth=2.2)
axs[0].set_ylabel('Entropia di von Neumann S(ρ)')
axs[0].legend(fontsize=11)
axs[0].grid(True, alpha=0.3)

axs[1].plot(times, result_base.expect[1], label='β = 0', linewidth=2.2)
axs[1].plot(times, result_retro.expect[1], label=f'β = φ⁻²', linewidth=2.2)
axs[1].set_ylabel('Concurrence C(ρ)')
axs[1].set_xlabel('Tempo (unità arbitrarie)')
axs[1].legend(fontsize=11)
axs[1].grid(True, alpha=0.3)

plt.suptitle('RENASCENT-Q — Stabilizzazione negentropica retrocausale\n'
             'Proxy due qubit con modulazione SAW e decoerenza termica',
             fontsize=14, y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('renascent_q_persistent_entanglement.jpg', dpi=350, bbox_inches='tight')
plt.show()

# Info finali
print(f"β = φ⁻² ≈ {beta:.4f}")
print(f"Entropia finale (baseline): {result_base.expect[0][-1]:.4f}")
print(f"Entropia finale (con retro): {result_retro.expect[0][-1]:.4f}")
print(f"Concurrence finale (baseline): {result_base.expect[1][-1]:.4f}")
print(f"Concurrence finale (con retro): {result_retro.expect[1][-1]:.4f}")