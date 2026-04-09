"""
RENASCENT-Q: Delayed Choice Quantum Eraser Multiqubit + Logarithmic Negativity
Versione stabile per QuTiP 5.x (Colab)
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

print("QuTiP version:", qt.__version__)  # Per verificare la versione

beta = 0.382
tlist = np.linspace(0, 20, 250)

# Stati di base
zero = qt.basis(2, 0)
one  = qt.basis(2, 1)

# Stato entangled iniziale a 3 qubit
psi0 = (qt.tensor(zero, zero, zero) + qt.tensor(one, one, one)).unit()

# Operatori
sigma_x_signal = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
Pi_path1 = qt.tensor(qt.ket2dm(zero), qt.qeye(2), qt.qeye(2))

# Post-selezione (erasure)
phi_future = (zero + one).unit() / np.sqrt(2)

# Weak Value
A_w = (phi_future.dag() * Pi_path1 * psi0).tr() / (phi_future.dag() * psi0).tr()
visibility = np.abs(A_w.real)

print(f"Weak Value (which-path) : {A_w:.4f}")
print(f"Visibilità massima DCQE : {visibility:.4f}")

# Evoluzione
H = 0.8 * sigma_x_signal

c_ops = [np.sqrt(0.08) * qt.tensor(qt.destroy(2), qt.qeye(2), qt.qeye(2))]
c_ops_retro = [np.sqrt(beta) * qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2)) / np.sqrt(2)]

log_neg_std = []
log_neg_retro = []
expect_signal = []
expect_retro = []

for t in tlist:
    rho_std = qt.mesolve(H, psi0, [t], c_ops, []).states[-1]
    rho_retro = qt.mesolve(H, psi0, [t], c_ops + c_ops_retro, []).states[-1]
    
    expect_signal.append(qt.expect(sigma_x_signal, rho_std))
    expect_retro.append(qt.expect(sigma_x_signal, rho_retro))
    
    # Logarithmic Negativity (corretta per QuTiP 5.x)
    log_neg_std.append(qt.logarithmic_negativity(rho_std))
    log_neg_retro.append(qt.logarithmic_negativity(rho_retro))

# Plot
fig, axs = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

axs[0].plot(tlist, expect_signal, 'r--', linewidth=2, label='Standard (β=0)')
axs[0].plot(tlist, expect_retro, 'b-', linewidth=3, label=f'Retrocausale β={beta}')
axs[0].axhline(y=visibility, color='g', linestyle=':', linewidth=2.5, 
               label=f'Visibilità DCQE = {visibility:.3f}')
axs[0].set_ylabel(r'$\langle \sigma_x \rangle$ sul qubit segnale')
axs[0].legend()
axs[0].grid(True, alpha=0.3)

axs[1].plot(tlist, log_neg_std, 'r--', linewidth=2, label='Standard')
axs[1].plot(tlist, log_neg_retro, 'b-', linewidth=3, label=f'Retrocausale (β={beta})')
axs[1].set_xlabel('Tempo')
axs[1].set_ylabel(r'Logarithmic Negativity $E_{\mathcal{N}}(\rho)$')
axs[1].legend()
axs[1].grid(True, alpha=0.3)

plt.suptitle('Delayed Choice Quantum Eraser con Kick Retrocausale-Negentropico\n(RENASCENT-Q Framework)')
plt.tight_layout()
plt.savefig('dcqe_multiqubit_with_logneg.pdf', dpi=300, bbox_inches='tight')
plt.show()

print("Figura salvata correttamente come: dcqe_multiqubit_with_logneg.pdf")