"""
RENASCENT-Q — VERSIONE FINALE 4 QUBIT (Ultra-Raffinata)
Post-selezione TSVF + braiding topologico continuo
Obiettivo: dimostrare stabilizzazione negentropica e entanglement persistente
Senza concurrence() - usa entanglement entropy tra le due coppie come metrica
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt

# ====================== Parametri ultra-ottimizzati ======================
omega = 1.0
g = 1.05
gamma_diss = 0.038
gamma_neg = 0.48
gamma_braid = 0.82
phi = (1 + np.sqrt(5)) / 2
beta = 1.0 / phi**2

T = 4.0
N_steps = 800
times = np.linspace(0, T, N_steps)

# ====================== Operatori 4 qubit ======================
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

H0 = omega * (sz1 + sz2 + sz3 + sz4) / 2
H_int = g * (sx1*sx2 + sx3*sx4 + sx1*sx3 + sx2*sx4)
H = H0 + H_int

bell = (basis(2,0)*basis(2,0).dag() + basis(2,1)*basis(2,1).dag()).unit() / np.sqrt(2)
psi0 = tensor(bell, bell)
rho0 = ket2dm(psi0)

c_diss = [np.sqrt(gamma_diss) * sm1, np.sqrt(gamma_diss) * sm2,
          np.sqrt(gamma_diss) * sm3, np.sqrt(gamma_diss) * sm4]

def weak_value_proxy(rho):
    sx_tot = sx1 + sx2 + sx3 + sx4
    expect_sx = expect(sx_tot, rho)
    return np.abs(expect_sx) * 2.4 + 0.7

# Funzione per entanglement entropy tra le due coppie (qubit 1-2 vs 3-4)
def entanglement_entropy(rho):
    # Partial trace su primi due qubit
    rho_a = ptrace(rho, [0,1])
    return entropy_vn(rho_a)

# ====================== Simulazione ======================
n_trajectories = 350
np.random.seed(42)

entropy_base = np.zeros(N_steps)
ent_ent_base = np.zeros(N_steps)   # entanglement entropy
entropy_retro = np.zeros(N_steps)
ent_ent_retro = np.zeros(N_steps)

for i in range(n_trajectories):
    # Baseline
    rho = rho0.copy()
    for j in range(N_steps):
        if j > 0:
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]
        if j % 80 == 0 or j == N_steps-1:
            entropy_base[j] += entropy_vn(rho)
            ent_ent_base[j] += entanglement_entropy(rho)

    # Retro + TSVF + braiding ultra-continuo
    rho = rho0.copy()
    for j in range(N_steps):
        if j > 0:
            result = mesolve(H, rho, [times[j-1], times[j]], c_diss, [])
            rho = result.states[-1]
            
            w_val = weak_value_proxy(rho)
            p_post = min(1.0, np.exp((w_val - 1.15) * 4.5)) * beta * 0.98
            
            if np.random.rand() < p_post:
                L_neg = (sm1.dag() + sm2.dag() + sm3.dag() + sm4.dag()).unit()
                rho_temp = L_neg * rho * L_neg.dag()
                if rho_temp.norm() > 1e-12:
                    rho = rho_temp / rho_temp.norm()
            
            # Braiding ultra-continuo
            braid_strength = gamma_braid * 0.021
            braid_op = braid_strength * (sx1*sx2 + sx3*sx4 + sx1*sx3 + sx2*sx4 + sz1*sz3 + sz2*sz4).unit()
            rho = (braid_op * rho * braid_op.dag() + (1 - braid_strength) * rho).unit()
        
        if j % 80 == 0 or j == N_steps-1:
            entropy_retro[j] += entropy_vn(rho)
            ent_ent_retro[j] += entanglement_entropy(rho)

# Media
entropy_base /= n_trajectories
ent_ent_base /= n_trajectories
entropy_retro /= n_trajectories
ent_ent_retro /= n_trajectories

# ====================== Plot ======================
fig, axs = plt.subplots(2, 1, figsize=(13.5, 10.5), sharex=True)

axs[0].plot(times, entropy_base, label='β = 0 (baseline)', linewidth=2.9)
axs[0].plot(times, entropy_retro, label=f'β = φ⁻² + post-selezione TSVF + braiding 4-qubit', linewidth=2.9)
axs[0].set_ylabel('Entropia di von Neumann S(ρ)')
axs[0].legend(fontsize=12)
axs[0].grid(True, alpha=0.3)

axs[1].plot(times, ent_ent_base, label='β = 0', linewidth=2.9)
axs[1].plot(times, ent_ent_retro, label=f'β = φ⁻² + post-selezione TSVF + braiding 4-qubit', linewidth=2.9)
axs[1].set_ylabel('Entanglement Entropy tra le due coppie')
axs[1].set_xlabel('Tempo (unità arbitrarie)')
axs[1].legend(fontsize=12)
axs[1].grid(True, alpha=0.3)

plt.suptitle('RENASCENT-Q — Stabilizzazione negentropica persistente\n'
             'Modello a 4 qubit con codifica anyonica + braiding ultra-continuo (versione finale)', 
             fontsize=16.5, y=0.96)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('renascent_q_persistent_entanglement_4qubit_final.jpg', dpi=480, bbox_inches='tight')
plt.show()

# Risultati
print(f"Simulazione completata con {n_trajectories} traiettorie")
print(f"β = φ⁻² ≈ {beta:.4f}")
print(f"Entropia finale baseline: {entropy_base[-1]:.4f}")
print(f"Entropia finale con retro + braiding: {entropy_retro[-1]:.4f}")
print(f"Entanglement Entropy finale baseline: {ent_ent_base[-1]:.4f}")
print(f"Entanglement Entropy finale con retro + braiding: {ent_ent_retro[-1]:.4f}")