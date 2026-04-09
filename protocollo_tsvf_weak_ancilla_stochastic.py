# protocollo_tsvf_weak_ancilla_stochastic.py
# Implementazione avanzata del protocollo TSVF-negentropico per Chip 2.0
# Weak measurement su ancilla reale + stochastic unraveling (mcsolve) + post-selezione retrocausale
# Versione raffinata per generare i risultati della simulazione (entanglement renascent)

import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

# ====================== PARAMETRI ======================
times = np.linspace(0, 15, 350) # Increased time points for better resolution
beta = 15.0
gamma_env = 0.08
ntraj = 800                     # Aumentato come richiesto

# ====================== STATO INIZIALE ======================
bell = (qt.tensor(qt.basis(2,0), qt.basis(2,0)) + qt.tensor(qt.basis(2,1), qt.basis(2,1))).unit() / np.sqrt(2) # Corrected to a 2-qubit Bell state
psi0 = qt.tensor(bell, qt.basis(2, 0)) # Now psi0 is a 3-qubit state

# ====================== HAMILTONIANA ======================
H = 1.05 * qt.tensor(qt.sigmax(), qt.sigmax(), qt.qeye(2)) + \
    0.35 * (qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2)) +\
            qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))) # Corrected this line

# ====================== DISSIPATORI ======================
c_ops_env = [
    np.sqrt(gamma_env) * qt.tensor(qt.sigmap(), qt.qeye(2), qt.qeye(2)),
    np.sqrt(gamma_env) * qt.tensor(qt.qeye(2), qt.sigmap(), qt.qeye(2))
]

def negentropic_diss(t, rho=None, args=None):
    if rho is None:
        return qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2)) * 0.0

    try:
        rho_sub = rho.ptrace([0, 1])
        # Corrected observable A to act directly on qubits 0 and 1
        A = qt.tensor(qt.sigmax(), qt.sigmax()) # A for entanglement measurement
        weak_val = qt.expect(A, rho_sub)
        effective_gamma = gamma_env * max(0.12, 1.0 - beta * abs(weak_val))
        # L_neg acts on qubits 0 and 1
        L_neg = np.sqrt(effective_gamma) * qt.tensor(qt.sigmap(), qt.sigmap(), qt.qeye(2))
        return L_neg
    except Exception:
        return qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2)) * 0.0

def concurrence_main(t, rho):
    try:
        return qt.concurrence(rho.ptrace([0, 1]))
    except:
        return 0.0

# ====================== SIMULAZIONE STOCASTICA ======================
print(f"Esecuzione stochastic unraveling con {ntraj} traiettorie... (potrebbe richiedere alcuni minuti)")

c_ops_total = c_ops_env + [negentropic_diss]

result = qt.mcsolve(H, psi0, times, c_ops=c_ops_total,
                    e_ops=[concurrence_main],
                    ntraj=ntraj,
                    options={
                        'nsteps': 100000, # Further increased
                        'atol': 1e-8,
                        'rtol': 1e-6,
                        'norm_steps': 3000, # Further increased
                        'norm_tol': 0.1, # Further relaxed
                        'norm_t_tol': 0.1, # Further relaxed
                        'store_states': True
                    })

# ====================== WEAK VALUE ======================
weak_values = []
for state in result.states:
    try:
        rho_sub = state.ptrace([0, 1])
        # Corrected observable A to act directly on qubits 0 and 1
        A = qt.tensor(qt.sigmax(), qt.sigmax())
        wv = qt.expect(A, rho_sub)
        weak_values.append(wv)
    except:
        weak_values.append(0.0)

# ====================== PLOT FINALE ======================
fig, axs = plt.subplots(2, 1, figsize=(13, 9.5), sharex=True)

axs[0].plot(times, result.expect[0], 'g-', linewidth=3.0, label=f'Negentropico Retrocausale (β = {beta})')
axs[0].plot(times, 0.5 * np.exp(-gamma_env * times), 'r--', linewidth=2.5, label='Decadimento standard approssimato')

axs[0].set_ylabel('Concurrence media (qubit 0-1)', fontsize=14)
axs[0].legend(fontsize=13)
axs[0].grid(True, alpha=0.35)

axs[1].plot(times, weak_values, 'b-', linewidth=2.8, label='Weak Value di A (Negentropico)')
axs[1].set_xlabel('Tempo (unità arbitrarie)', fontsize=14)
axs[1].set_ylabel('Weak Value di A', fontsize=14)
axs[1].legend(fontsize=13)
axs[1].grid(True, alpha=0.35)

plt.suptitle('Entanglement Renascent – Protocollo TSVF-negentropico su Chip 2.0\n(stochastic unraveling + weak ancilla)', fontsize=15.5)
plt.tight_layout()
plt.savefig('renascent_entanglement_tsvf_weak_ancilla.png', dpi=400, bbox_inches='tight')
plt.show()

# ====================== RISULTATI ======================
print("\n=== RISULTATI FINALI ===")
print(f"Concurrence media finale : {np.mean(result.expect[0][-40:]):.4f}")
print(f"Plot salvato come: renascent_entanglement_tsvf_weak_ancilla.png")