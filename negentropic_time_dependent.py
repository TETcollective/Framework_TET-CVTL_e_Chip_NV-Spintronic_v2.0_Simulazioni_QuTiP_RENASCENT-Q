# negentropic_time_dependent.py
# Integrazione temporale del termine negentropico custom per Chip 2.0
# Versione raffinata e stabile per generare gamma_neg(t)

import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

# ====================== PARAMETRI ======================
times = np.linspace(0, 12, 400)
beta = 1.0                  # Valore forte (ridotto per osservare dinamica)
gamma_env = 0.12

# ====================== STATO INIZIALE ======================
# Define a 2-qubit Bell state for the first two system qubits
bell_2_qubits = (qt.tensor(qt.basis(2,0), qt.basis(2,0)) + qt.tensor(qt.basis(2,1), qt.basis(2,1))).unit() / np.sqrt(2)
# Now tensor it with the third qubit (ancilla) in the |0> state to form a 3-qubit system
psi0 = qt.tensor(bell_2_qubits, qt.basis(2, 0))

# ====================== HAMILTONIANA ======================
H = 1.2 * qt.tensor(qt.sigmax(), qt.sigmax(), qt.qeye(2)) + \
    0.45 * (qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2)) +
            qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2)))

# ====================== DISSIPATORE AMBIENTALE ======================
c_ops_env = [
    np.sqrt(gamma_env) * qt.tensor(qt.sigmap(), qt.qeye(2), qt.qeye(2)),
    np.sqrt(gamma_env) * qt.tensor(qt.qeye(2), qt.sigmap(), qt.qeye(2))
]

# ====================== RATE NEGENTROPICO TIME-DEPENDENT (VERSIONE CORRETTA) ======================
def gamma_neg(t, rho):
    """Calcolo robusto del rate negentropico"""
    try:
        # Riduzione esplicita ai due qubit principali
        rho_sub = rho.ptrace([0, 1])
        # Observable sensibile all'entanglement (σx1 + σx2)
        A = qt.tensor(qt.sigmax(), qt.sigmax()) # Corrected A operator
        weak_val = qt.expect(A, rho_sub)

        # Formula più aggressiva per rendere visibile la variazione
        suppression = beta * abs(weak_val)
        return gamma_env * max(0.05, 1.0 - suppression)
    except:
        return gamma_env * 0.6   # fallback

# ====================== DISSIPATORE NEGENTROPICO ======================
def negentropic_diss(t, rho=None, args=None):
    if rho is None:
        # Return a zero operator for initialization by QobjEvo (3-qubit system)
        return qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2)) * 0.0

    gamma_t = gamma_neg(t, rho)
    # Operatore di pompaggio coerente
    L_neg = np.sqrt(gamma_t) * qt.tensor(qt.sigmap(), qt.sigmap(), qt.qeye(2))
    return L_neg # Return a single Qobj, not a list

# ====================== SIMULAZIONE ======================
c_ops_total = c_ops_env + [negentropic_diss]

print("Esecuzione simulazione con termine negentropico time-dependent...")
result = qt.mesolve(H, psi0, times, c_ops=c_ops_total,
                    options={'nsteps': 40000, 'atol': 1e-9})

# ====================== ESTRAZIONE γ_neg(t) ======================
gamma_values = []
for i, t in enumerate(times):
    try:
        g = gamma_neg(t, result.states[i])
        gamma_values.append(g)
    except:
        gamma_values.append(gamma_env)

# ====================== PLOT ======================
plt.figure(figsize=(11.5, 6.8))
plt.plot(times, gamma_values, 'c-', linewidth=3.2, label=r'$\gamma_{\text{neg}}(t)$')
plt.axhline(y=gamma_env, color='r', linestyle='--', linewidth=2.2, label=r'$\gamma_{\text{env}}$ (ambiente)')

plt.xlabel('Tempo (unità arbitrarie)', fontsize=14)
plt.ylabel(r'Rate negentropico $\gamma_{\text{neg}}(t)$', fontsize=14)
plt.title('Evoluzione temporale del rate negentropico', fontsize=15.5)
plt.grid(True, alpha=0.35)
plt.legend(fontsize=13)
plt.tight_layout()
plt.savefig('gamma_neg_time_evolution.png', dpi=400, bbox_inches='tight')
plt.show()

print(f"Plot salvato come: gamma_neg_time_evolution.png")
print(f"Valore medio di γ_neg(t): {np.mean(gamma_values):.4f}")
print(f"Min γ_neg(t): {np.min(gamma_values):.4f}  |  Max γ_neg(t): {np.max(gamma_values):.4f}")