"""
RENASCENT-Q: Four-Point Correlator + Logarithmic Negativity + Vacuum Torque Extraction
Codice principale per i risultati numerici del paper
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

# ==================== Parametri ====================
beta = 0.382
gamma = 0.065
omega_saw = 1.8
tlist = np.linspace(0, 35, 450)

# ==================== Operatori su 4 qubit ====================
sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2), qt.qeye(2))
sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2), qt.qeye(2))
sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz(), qt.qeye(2))
sz4 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2), qt.sigmaz())

# Stato iniziale fortemente entangled
psi0 = (qt.basis(16, 0) + qt.basis(16, 15)).unit()

# Hamiltoniano con modulazione SAW
def H_t(t):
    H0 = 0.35 * (sz1*sz2 + sz2*sz3 + sz3*sz4)
    return H0 + 0.18 * np.sin(omega_saw * t) * (sz1 + sz4)

# Collassatori decoerenza
c_ops = [np.sqrt(gamma) * qt.tensor(qt.destroy(2), qt.qeye(2), qt.qeye(2), qt.qeye(2)) for _ in range(4)]

# Termine retrocausale-negentropico
O_neg = (sz1 + sz2 + sz3 + sz4) / 2.0
c_ops_retro = [np.sqrt(beta) * O_neg]

# Liste risultati
four_point_retro = []
log_neg_retro = []
torque = []

print("Esecuzione simulazione principale RENASCENT-Q...")

for t in tlist:
    rho = qt.mesolve(H_t(t), psi0, [t], c_ops + c_ops_retro, []).states[-1]
    
    # Correlatore four-point
    C4 = qt.expect(sz1*sz2*sz3*sz4, rho)
    four_point_retro.append(np.abs(C4))
    
    # Logarithmic Negativity
    E_N = qt.logarithmic_negativity(rho)
    log_neg_retro.append(E_N)
    
    # Estrazione Vacuum Torque
    tau_vac = np.real(C4) * np.sin(omega_saw * t)
    torque.append(tau_vac)

# ==================== Plot principale ====================
fig, axs = plt.subplots(3, 1, figsize=(12, 11), sharex=True)

# Pannello 1: Four-point correlator
axs[0].plot(tlist, four_point_retro, 'b-', linewidth=3)
axs[0].set_ylabel(r'$|C^{(4)}_{\rm retro}|$')
axs[0].set_title('Correlatore Four-Point Retrocausale')
axs[0].grid(True, alpha=0.3)

# Pannello 2: Logarithmic Negativity
axs[1].plot(tlist, log_neg_retro, 'b-', linewidth=3)
axs[1].set_ylabel(r'$E_{\mathcal{N}}(\rho)$')
axs[1].set_title('Logarithmic Negativity')
axs[1].grid(True, alpha=0.3)

# Pannello 3: Vacuum Torque
axs[2].plot(tlist, torque, 'm-', linewidth=3)
axs[2].set_xlabel('Tempo')
axs[2].set_ylabel(r'$\tau_{\rm vac}(t)$')
axs[2].set_title('Estrazione di Vacuum Torque')
axs[2].grid(True, alpha=0.3)

plt.suptitle('RENASCENT-Q: Correlatori Four-Point, Logarithmic Negativity e Vacuum Torque', 
             fontsize=14, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('renascent_fourpoint_logneg_torque.pdf', dpi=300, bbox_inches='tight')
plt.show()

print("Figura salvata come: renascent_fourpoint_logneg_torque.pdf")
print(f"Log Neg finale: {log_neg_retro[-1]:.4f}")
print(f"Vacuum torque medio: {np.mean(np.abs(torque)):.4f}")