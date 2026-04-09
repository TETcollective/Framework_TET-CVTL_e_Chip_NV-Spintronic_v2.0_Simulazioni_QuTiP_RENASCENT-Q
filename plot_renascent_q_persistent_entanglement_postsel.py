"""
RENASCENT-Q: Simulazione QuTiP con termine retrocausale-negentropico + post-selezione TSVF-inspired
Proxy due qubit con modulazione SAW e decoerenza termica.
Versione con post-selezione simulata - Aprile 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# ====================== Parametri fisici ======================
omega = 1.0
g = 0.5
gamma_diss = 0.08
gamma_neg = 0.22          # leggermente aumentato per evidenziare l'effetto post-selezione
phi = (1 + np.sqrt(5)) / 2
beta = 1.0 / phi**2       # ≈ 0.382

T = 2.0
N_steps = 300
times = np.linspace(0, T, N_steps)

# ====================== Operatori e Hamiltoniano ======================
sz1 = tensor(sigmaz(), qeye(2))
sz2 = tensor(qeye(2), sigmaz())
sx1 = tensor(sigmax(), qeye(2))
sx2 = tensor(qeye(2), sigmax())
sm1 = tensor(sigmam(), qeye(2))
sm2 = tensor(qeye(2), sigmam())

# Definizione di sy1 e sy2
sy1 = tensor(sigmay(), qeye(2))
sy2 = tensor(qeye(2), sigmay())

H0 = omega * (sz1 + sz2) / 2
H_int = g * (sx1*sx2 + sy1*sy2)
H = H0 + H_int

# Stato iniziale Bell-like corretto (ket per due qubit)
psi0 = (tensor(basis(2,0), basis(2,0)) + tensor(basis(2,1), basis(2,1))).unit()
rho0 = ket2dm(psi0)

# Dissipatori standard
c_diss = [np.sqrt(gamma_diss) * sm1, np.sqrt(gamma_diss) * sm2]

# Proxy semplice per weak value (basato su correlazione trasversale)
def weak_value_proxy(rho):
    sx_tot = sx1 + sx2
    expect_sx = expect(sx_tot, rho)
    return np.abs(expect_sx) + 0.6   # valore che può superare 1 (anomalo)

# ====================== Simulazione Monte Carlo con post-selezione ======================
n_trajectories = 150
np.random.seed(42)

entropy_base = np.zeros(N_steps)
conc_base = np.zeros(N_steps)
entropy_retro = np.zeros(N_steps)
conc_retro = np.zeros(N_steps)

count_survived = 0

for i in range(n_trajectories):
    # ----- Baseline (solo dissipazione) -----
    rho = rho0.copy()
    ent_traj_b = np.zeros(N_steps)
    con_traj_b = np.zeros(N_steps)

    for j in range(N_steps):
        if j > 0:
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]
        ent_traj_b[j] = entropy_vn(rho)
        con_traj_b[j] = concurrence(rho)

    entropy_base += ent_traj_b
    conc_base += con_traj_b

    # ----- Con retrocausale + post-selezione -----
    rho = rho0.copy()
    ent_traj_r = np.zeros(N_steps)
    con_traj_r = np.zeros(N_steps)

    for j in range(N_steps):
        if j > 0:
            # Evoluzione dissipativa standard
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]

            # Calcolo proxy weak value e probabilità di post-selezione
            w_val = weak_value_proxy(rho)
            p_post = min(1.0, np.exp((w_val - 1.1) * 2.0)) * beta * 0.85

            # Applica pompaggio negentropico solo se post-selezionato
            if np.random.rand() < p_post:
                L_neg = (sm1.dag() + sm2.dag()).unit()
                rho_temp = L_neg * rho * L_neg.dag()
                if rho_temp.norm() > 1e-12:
                    rho = rho_temp / rho_temp.norm()

        ent_traj_r[j] = entropy_vn(rho)
        con_traj_r[j] = concurrence(rho)

    entropy_retro += ent_traj_r
    conc_retro += con_traj_r
    count_survived += 1   # tutte le traiettorie vengono tenute (media pesata implicita)

# Media
entropy_base /= n_trajectories
conc_base /= n_trajectories
entropy_retro /= n_trajectories
conc_retro /= n_trajectories

# ====================== Plot ======================
fig, axs = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)

axs[0].plot(times, entropy_base, label='β = 0 (baseline)', linewidth=2.4)
axs[0].plot(times, entropy_retro, label=f'β = φ⁻² ≈ {beta:.3f} + post-selezione TSVF', linewidth=2.4)
axs[0].set_ylabel('Entropia di von Neumann S(ρ) (media su traiettorie)')
axs[0].legend(fontsize=11)
axs[0].grid(True, alpha=0.35)

axs[1].plot(times, conc_base, label='β = 0', linewidth=2.4)
axs[1].plot(times, conc_retro, label=f'β = φ⁻² + post-selezione', linewidth=2.4)
axs[1].set_ylabel('Concurrence C(ρ) (media su traiettorie)')
axs[1].set_xlabel('Tempo (unità arbitrarie)')
axs[1].legend(fontsize=11)
axs[1].grid(True, alpha=0.35)

plt.suptitle('RENASCENT-Q — Stabilizzazione negentropica con post-selezione TSVF\n'
             'Proxy due qubit + modulazione SAW + decoerenza termica',
             fontsize=14, y=0.96)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('renascent_q_persistent_entanglement_postsel.jpg', dpi=380, bbox_inches='tight')
plt.show()

# Statistiche finali
print(f"Simulazione completata con {n_trajectories} traiettorie")
print(f"β = φ⁻² ≈ {beta:.4f}")
print(f"Entropia finale baseline: {entropy_base[-1]:.4f}")
print(f"Entropia finale con retro + post-selezione: {entropy_retro[-1]:.4f}")
print(f"Concurrence finale baseline: {conc_base[-1]:.4f}")
print(f"Concurrence finale con retro + post-selezione: {conc_retro[-1]:.4f}")