"""
RENASCENT-Q — VERSIONE ULTRA RAFFINATA FINALE
Post-selezione TSVF + braiding topologico ultra-continuo e calibrato
Obiettivo: massimizzare concurrence persistente nel proxy a due qubit
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# ====================== Parametri ultra-ottimizzati ======================
omega = 1.0
g = 0.95                    # accoppiamento molto forte per mantenere entanglement
gamma_diss = 0.048          # dissipazione ridotta al minimo realistico
gamma_neg = 0.45            # rate negentropico alto
gamma_braid = 0.75          # braiding continuo forte ma stabile
phi = (1 + np.sqrt(5)) / 2
beta = 1.0 / phi**2         # ≈ 0.382

T = 3.5
N_steps = 700
times = np.linspace(0, T, N_steps)

# ====================== Operatori ======================
sz1 = tensor(sigmaz(), qeye(2))
sz2 = tensor(qeye(2), sigmaz())
sx1 = tensor(sigmax(), qeye(2))
sx2 = tensor(qeye(2), sigmax())
sm1 = tensor(sigmam(), qeye(2))
sm2 = tensor(qeye(2), sigmam())

H0 = omega * (sz1 + sz2) / 2
H_int = g * (sx1*sx2 + sx1*sx2)   # accoppiamento trasverso molto forte
H = H0 + H_int

# Stato iniziale fortemente entangled (corretto per essere un vettore ket a due qubit)
psi0 = (tensor(basis(2,0), basis(2,0)) + tensor(basis(2,1), basis(2,1))).unit()
rho0 = ket2dm(psi0)

c_diss = [np.sqrt(gamma_diss) * sm1, np.sqrt(gamma_diss) * sm2]

def weak_value_proxy(rho):
    sx_tot = sx1 + sx2
    expect_sx = expect(sx_tot, rho)
    return np.abs(expect_sx) * 2.3 + 0.65   # weak value fortemente anomalo

# ====================== Simulazione Monte Carlo ======================
n_trajectories = 400
np.random.seed(42)

entropy_base = np.zeros(N_steps)
conc_base = np.zeros(N_steps)
entropy_retro = np.zeros(N_steps)
conc_retro = np.zeros(N_steps)

for i in range(n_trajectories):
    # ----- Baseline -----
    rho = rho0.copy()
    for j in range(N_steps):
        if j > 0:
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]
        entropy_base[j] += entropy_vn(rho)
        conc_base[j] += concurrence(rho)

    # ----- Retro + TSVF + braiding ultra-continuo -----
    rho = rho0.copy()
    for j in range(N_steps):
        if j > 0:
            # Evoluzione dissipativa
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]

            # Post-selezione TSVF calibrata
            w_val = weak_value_proxy(rho)
            p_post = min(1.0, np.exp((w_val - 1.18) * 4.1)) * beta * 0.98

            if np.random.rand() < p_post:
                L_neg = (sm1.dag() + sm2.dag()).unit()
                rho_temp = L_neg * rho * L_neg.dag()
                if rho_temp.norm() > 1e-12:
                    rho = rho_temp / rho_temp.norm()

            # Braiding ultra-continuo e calibrato (applicato ogni passo con forza ottimale)
            braid_strength = gamma_braid * 0.017   # valore piccolo ma costante
            braid_op = braid_strength * (tensor(sigmax(), sigmax()) + tensor(sigmay(), sigmay())).unit()
            rho = (braid_op * rho * braid_op.dag() + (1 - braid_strength) * rho).unit()

        entropy_retro[j] += entropy_vn(rho)
        conc_retro[j] += concurrence(rho)

# Media
entropy_base /= n_trajectories
conc_base /= n_trajectories
entropy_retro /= n_trajectories
conc_retro /= n_trajectories

# ====================== Plot ======================
fig, axs = plt.subplots(2, 1, figsize=(13, 10), sharex=True)

axs[0].plot(times, entropy_base, label='β = 0 (baseline)', linewidth=2.8)
axs[0].plot(times, entropy_retro, label=f'β = φ⁻² + post-selezione TSVF + braiding ultra-continuo', linewidth=2.8)
axs[0].set_ylabel('Entropia di von Neumann S(ρ) (media)')
axs[0].legend(fontsize=12)
axs[0].grid(True, alpha=0.3)

axs[1].plot(times, conc_base, label='β = 0', linewidth=2.8)
axs[1].plot(times, conc_retro, label=f'β = φ⁻² + post-selezione TSVF + braiding ultra-continuo', linewidth=2.8)
axs[1].set_ylabel('Concurrence C(ρ) (media)')
axs[1].set_xlabel('Tempo (unità arbitrarie)')
axs[1].legend(fontsize=12)
axs[1].grid(True, alpha=0.3)

plt.suptitle('RENASCENT-Q — Stabilizzazione negentropica persistente\n'
             'Post-selezione TSVF + braiding topologico ultra-continuo (versione finale raffinata)',
             fontsize=16, y=0.96)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('renascent_q_persistent_entanglement_ultra_refined.jpg', dpi=450, bbox_inches='tight')
plt.show()

# Risultati finali
print(f"Simulazione completata con {n_trajectories} traiettorie")
print(f"β = φ⁻² ≈ {beta:.4f}")
print(f"Entropia finale baseline: {entropy_base[-1]:.4f}")
print(f"Entropia finale con retro + braiding: {entropy_retro[-1]:.4f}")
print(f"Concurrence finale baseline: {conc_base[-1]:.4f}")
print(f"Concurrence finale con retro + braiding: {conc_retro[-1]:.4f}")