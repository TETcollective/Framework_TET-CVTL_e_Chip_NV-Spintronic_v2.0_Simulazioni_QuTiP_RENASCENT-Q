# protocollo_tsvf_4qubit_anyonic_moire.py
# Evoluzione temporale del weak value ⟨A⟩_w con braiding moiré esplicito
# Versione raffinata per Chip 2.0 - 4 qubit con imprint anyonico

import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

# ====================== PARAMETRI ======================
times = np.linspace(0, 15, 300)
beta = 1.0 # Changed beta from 18.0 to 1.0 based on context for dynamic variation
gamma_env = 0.07
ntraj = 600

# ====================== STATO INIZIALE ======================
# bell dovrebbe essere uno stato di Bell a 2 qubit per un sistema a 4 qubit complessivo
bell = (qt.tensor(qt.basis(2,0), qt.basis(2,0)) + qt.tensor(qt.basis(2,1), qt.basis(2,1))).unit() / np.sqrt(2)
# 4 qubit: 2 sistema + 2 per braiding moiré semplificato
psi0 = qt.tensor([bell, qt.basis(2,0), qt.basis(2,0)])

# ====================== HAMILTONIANA con braiding moiré ======================
# Interazione XX tra qubit 0-1 + termine moiré semplificato tra 2-3
H = 1.2 * (qt.tensor(qt.sigmax(), qt.sigmax(), qt.qeye(2), qt.qeye(2)) +
           qt.tensor(qt.qeye(2), qt.sigmax(), qt.sigmax(), qt.qeye(2))) + \
    0.5 * qt.tensor(qt.sigmaz(), qt.qeye(2), qt.sigmaz(), qt.qeye(2))

# ====================== DISSIPATORI AMBIENTALI ======================
c_ops_env = [
    np.sqrt(gamma_env) * qt.tensor(qt.sigmap(), qt.qeye(2), qt.qeye(2), qt.qeye(2)),
    np.sqrt(gamma_env) * qt.tensor(qt.qeye(2), qt.sigmap(), qt.qeye(2), qt.qeye(2))
]

# ====================== DISSIPATORE NEGENTROPICO ======================
def negentropic_diss(t, rho=None, args=None):
    if rho is None:
        # Return a zero operator for initialization by QobjEvo (4-qubit system)
        return qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2), qt.qeye(2)) * 0.0
    try:
        rho_sub = rho.ptrace([0, 1])  # solo i due qubit principali
        A = qt.tensor(qt.sigmax(), qt.sigmax()) # Changed A operator for weak value sensitivity
        weak_val = qt.expect(A, rho_sub)
        effective_gamma = gamma_env * max(0.15, 1.0 - beta * abs(weak_val))
        L_neg = np.sqrt(effective_gamma) * qt.tensor(qt.sigmap(), qt.sigmap(), qt.qeye(2), qt.qeye(2))
        return L_neg
    except:
        return qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2), qt.qeye(2)) * 0.0 # Return a zero operator on error

# ====================== WEAK VALUE OBSERVABLE ======================
def weak_value_A(t, rho):
    try:
        rho_sub = rho.ptrace([0, 1])
        A = qt.tensor(qt.sigmax(), qt.sigmax()) # Changed A operator for weak value sensitivity
        return qt.expect(A, rho_sub)
    except:
        return 0.0

# ====================== SIMULAZIONE ======================
print(f"Esecuzione simulazione con {ntraj} traiettorie e braiding moiré...")

c_ops_total = c_ops_env + [negentropic_diss]

result = qt.mcsolve(H, psi0, times, c_ops=c_ops_total,
                    e_ops=[concurrence_main := lambda t, rho: qt.concurrence(rho.ptrace([0,1])),
                           weak_value_A],
                    ntraj=ntraj,
                    options={
                        'nsteps': 100000,
                        'atol': 1e-8,
                        'rtol': 1e-6,
                        'norm_steps': 3000,
                        'norm_tol': 0.1,  # Lowered tolerance (less strict)
                        'norm_t_tol': 0.1, # Lowered tolerance (less strict)
                        'store_states': True
                    })

# ====================== PLOT ======================
fig, axs = plt.subplots(2, 1, figsize=(13, 9.5), sharex=True)

axs[0].plot(times, result.expect[0], 'g-', linewidth=3.0, label='Concurrence (Negentropico)')
axs[0].plot(times, 0.5 * np.exp(-gamma_env * times), 'r--', linewidth=2.5, label='Decadimento standard')
axs[0].set_ylabel('Concurrence (qubit 0-1)', fontsize=14)
axs[0].legend(fontsize=13)
axs[0].grid(True, alpha=0.35)

axs[1].plot(times, result.expect[1], 'b-', linewidth=2.8, label='Weak Value di A (Reale)')
axs[1].set_xlabel('Tempo (unità arbitrarie)', fontsize=14)
axs[1].set_ylabel(r'Weak Value $\langle A \rangle_w$', fontsize=14)
axs[1].legend(fontsize=13)
axs[1].grid(True, alpha=0.35)

plt.suptitle('Evoluzione temporale del weak value ⟨A⟩_w con braiding moiré', fontsize=15.5)
plt.tight_layout()
plt.savefig('weak_value_evolution_moire.png', dpi=400, bbox_inches='tight')
plt.show()

print("\n=== RISULTATI ===")
print(f"Concurrence media finale: {np.mean(result.expect[0][-40:]):.4f}")
print("Plot salvato come: weak_value_evolution_moire.png")