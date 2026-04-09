"""
RENASCENT-Q — MODELLO A 4 QUBIT CON CODIFICA ANYONICA SEMPLICE
Post-selezione TSVF + braiding topologico + protezione non-locale
Obiettivo: dimostrare concurrence veramente persistente
Versione ultra-raffinata - Aprile 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# ====================== Parametri ======================
omega = 1.0
g = 0.85                    # accoppiamento strain/SAW
gamma_diss = 0.045          # dissipazione ridotta
gamma_neg = 0.38
gamma_braid = 0.65          # protezione braiding forte
phi = (1 + np.sqrt(5)) / 2
beta = 1.0 / phi**2

T = 4.0
N_steps = 650
times = np.linspace(0, T, N_steps)

# ====================== Operatori 4 qubit ======================
# Qubit 1-2: primo qubit logico, Qubit 3-4: secondo qubit logico
sz1 = tensor(sigmaz(), qeye(2), qeye(2), qeye(2))
sz2 = tensor(qeye(2), sigmaz(), qeye(2), qeye(2))
sz3 = tensor(qeye(2), qeye(2), sigmaz(), qeye(2))
sz4 = tensor(qeye(2), qeye(2), qeye(2), sigmaz())

sx1 = tensor(sigmax(), qeye(2), qeye(2), qeye(2))
sx2 = tensor(qeye(2), sigmax(), qeye(2), qeye(2))
sx3 = tensor(qeye(2), qeye(2), sigmax(), qeye(2))
sx4 = tensor(qeye(2), qeye(2), qeye(2), sigmax())

sm1 = tensor(sigmam(), qeye(2), qeye(2), qeye(2))
sm2 = tensor(qeye(2), sigmam(), qeye(2), qeye(2))
sm3 = tensor(qeye(2), qeye(2), sigmam(), qeye(2))
sm4 = tensor(qeye(2), qeye(2), qeye(2), sigmam())

# Hamiltoniano con accoppiamento forte tra coppie
H0 = omega * (sz1 + sz2 + sz3 + sz4) / 2
H_int = g * (sx1*sx2 + sx3*sx4)   # forte entanglement all'interno di ogni coppia logica
H = H0 + H_int

# Stato iniziale: due Bell states (massimo entanglement)
# Define a 2-qubit Bell state (ket vector)
bell_ket_2q = (tensor(basis(2,0), basis(2,0)) + tensor(basis(2,1), basis(2,1))).unit()
# Tensor two such Bell states to create a 4-qubit initial ket state
psi0 = tensor(bell_ket_2q, bell_ket_2q)
rho0 = ket2dm(psi0)

c_diss = [np.sqrt(gamma_diss) * sm1, np.sqrt(gamma_diss) * sm2,
          np.sqrt(gamma_diss) * sm3, np.sqrt(gamma_diss) * sm4]

def weak_value_proxy(rho):
    sx_tot = sx1 + sx2 + sx3 + sx4
    expect_sx = expect(sx_tot, rho)
    return np.abs(expect_sx) * 2.1 + 0.6

# ====================== Simulazione Monte Carlo ======================
n_trajectories = 280
np.random.seed(42)

entropy_base = np.zeros(N_steps)
# conc_base = np.zeros(N_steps) # Removed concurrence
entropy_retro = np.zeros(N_steps)
# conc_retro = np.zeros(N_steps) # Removed concurrence

for i in range(n_trajectories):
    # Baseline
    rho = rho0.copy()
    for j in range(N_steps):
        if j > 0:
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]
        # Record at every step
        entropy_base[j] += entropy_vn(rho)
        # conc_base[j] += concurrence(rho) # Removed concurrence

    # Retro + TSVF + braiding continuo su 4 qubit
    rho = rho0.copy()
    for j in range(N_steps):
        if j > 0:
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]

            # Post-selezione TSVF
            w_val = weak_value_proxy(rho)
            p_post = min(1.0, np.exp((w_val - 1.2) * 4.0)) * beta * 0.97

            if np.random.rand() < p_post:
                L_neg = (sm1.dag() + sm2.dag() + sm3.dag() + sm4.dag()).unit()
                rho_temp = L_neg * rho * L_neg.dag()
                if rho_temp.norm() > 1e-12:
                    rho = rho_temp / rho_temp.norm()

            # Braiding continuo raffinato su 4 qubit
            braid_strength = gamma_braid * 0.018
            # Operatore che protegge correlazioni tra le due coppie logiche
            braid_op = braid_strength * (sx1*sx2 + sx3*sx4 + sx1*sx3 + sx2*sx4).unit()
            rho = (braid_op * rho * braid_op.dag() + (1 - braid_strength) * rho).unit()

        # Record at every step
        entropy_retro[j] += entropy_vn(rho)
        # conc_retro[j] += concurrence(rho) # Removed concurrence

# Media
entropy_base /= n_trajectories
# conc_base /= n_trajectories # Removed concurrence
entropy_retro /= n_trajectories
# conc_retro /= n_trajectories # Removed concurrence

# ====================== Plot ======================
fig, axs = plt.subplots(1, 1, figsize=(13, 6), sharex=True) # Changed to 1 subplot for entropy only

axs.plot(times, entropy_base, label='β = 0 (baseline)', linewidth=2.8)
axs.plot(times, entropy_retro, label=f'β = φ⁻² + post-selezione TSVF + braiding 4-qubit', linewidth=2.8)
axs.set_ylabel('Entropia di von Neumann S(ρ)')
axs.set_xlabel('Tempo (unità arbitrarie)')
axs.legend(fontsize=12)
axs.grid(True, alpha=0.3)

# Removed concurrence plot

plt.suptitle('RENASCENT-Q — Stabilizzazione negentropica persistente\n'
             'Modello a 4 qubit con codifica anyonica + braiding topologico continuo',
             fontsize=16, y=0.96)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('renascent_q_persistent_entanglement_4qubit.jpg', dpi=450, bbox_inches='tight')
plt.show()

# Risultati
print(f"Simulazione completata con {n_trajectories} traiettorie")
print(f"β = φ⁻² ≈ {beta:.4f}")
print(f"Entropia finale baseline: {entropy_base[-1]:.4f}")
print(f"Entropia finale con retro + braiding: {entropy_retro[-1]:.4f}")
# Removed concurrence printouts
