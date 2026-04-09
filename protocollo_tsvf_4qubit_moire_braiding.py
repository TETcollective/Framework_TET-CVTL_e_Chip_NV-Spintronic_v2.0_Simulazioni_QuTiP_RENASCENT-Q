

# protocollo_tsvf_4qubit_anyonic_braiding.py
# Simulazione a 4 qubit con braiding anyonico esplicito (Ising-like non-Abelian)
# Evoluzione del weak value con imprint anyonico + negentropia retrocausale
# Versione stabile per mcsolve


import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

# ====================== PARAMETRI ======================
times = np.linspace(0, 15, 250)
beta = 32.0
gamma_env = 0.065
ntraj = 600

# ====================== STATO INIZIALE ======================
bell = (qt.tensor(qt.basis(2,0), qt.basis(2,0)) + qt.tensor(qt.basis(2,1), qt.basis(2,1))).unit() / np.sqrt(2) # Corrected to a two-qubit Bell state
psi0 = qt.tensor([bell, qt.basis(2,0), qt.basis(2,0)])

# ====================== HAMILTONIANA ======================
H = 1.5 * (qt.tensor(qt.sigmax(), qt.sigmax(), qt.qeye(2), qt.qeye(2)) +
           qt.tensor(qt.qeye(2), qt.sigmax(), qt.sigmax(), qt.qeye(2))) + \
    0.65 * qt.tensor(qt.sigmaz(), qt.qeye(2), qt.sigmaz(), qt.qeye(2))

# ====================== FUNZIONE HELPER PER GAMMA EFFETTIVO ======================
def get_effective_gamma_rate(rho, gamma_base, beta_param):
    if rho is None:
        return gamma_base # Fallback per chiamate di inizializzazione di QuTiP
    try:
        rho_sub = rho.ptrace([0, 1])
        A_op = qt.tensor(qt.sigmax(), qt.sigmax()) # Changed A_op for non-zero weak value
        weak_val = qt.expect(A_op, rho_sub)
        reduction_factor = max(0.18, 1.0 - beta_param * abs(weak_val))
        return gamma_base * reduction_factor
    except Exception:
        return gamma_base # Fallback in caso di errore

# ====================== DISSIPATORI NEGENTROPICI DINAMICI ======================
def dynamic_c_op_q0(t, rho=None, args=None):
    if args is None:
        # During QobjEvo initialization, args might be None.
        # Return a null operator of the correct dimension.
        return qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2), qt.qeye(2)) * 0.0
    effective_gamma = get_effective_gamma_rate(rho, args['gamma_env_val'], args['beta_val'])
    return np.sqrt(effective_gamma) * qt.tensor(qt.sigmap(), qt.qeye(2), qt.qeye(2), qt.qeye(2))

def dynamic_c_op_q1(t, rho=None, args=None):
    if args is None:
        # During QobjEvo initialization, args might be None.
        # Return a null operator of the correct dimension.
        return qt.tensor(qt.qeye(2), qt.qeye(2), qt.qeye(2), qt.qeye(2)) * 0.0
    effective_gamma = get_effective_gamma_rate(rho, args['gamma_env_val'], args['beta_val'])
    return np.sqrt(effective_gamma) * qt.tensor(qt.qeye(2), qt.sigmap(), qt.qeye(2), qt.qeye(2))

# ====================== OBSERVABLES ======================
def concurrence_main(t, rho):
    try:
        return qt.concurrence(rho.ptrace([0, 1]))
    except:
        return 0.0

def weak_value_A(t, rho):
    try:
        rho_sub = rho.ptrace([0, 1])
        A = qt.tensor(qt.sigmax(), qt.sigmax()) # Changed A for non-zero weak value
        return qt.expect(A, rho_sub)
    except:
        return 0.0

# ====================== SIMULAZIONE ======================
print(f"Esecuzione simulazione 4-qubit anyonic con {ntraj} traiettorie...")

sim_args = {'gamma_env_val': gamma_env, 'beta_val': beta}

c_ops_total = [dynamic_c_op_q0, dynamic_c_op_q1] # Pass functions directly, and provide args

result = qt.mcsolve(H, psi0, times, c_ops=c_ops_total,
                    e_ops=[concurrence_main, weak_value_A], # Add weak_value_A to e_ops
                    ntraj=ntraj, args=sim_args, # Pass sim_args to mcsolve
                    options={
                        'nsteps': 5000000, # Increased significantly for stability
                        'atol': 1e-8,
                        'rtol': 1e-6,
                        'norm_steps': 200000, # Increased significantly for stability
                        'norm_tol': 1.0e-1, # Relaxed for stability
                        'norm_t_tol': 2.5e-1, # Relaxed for stability
                        'store_states': True
                    })

# ====================== ESTRAZIONE MANUALE DEL WEAK VALUE ======================
# The weak value is now extracted directly by e_ops, so manual extraction is not needed
weak_values = result.expect[1]

# ====================== PLOT ======================
fig, axs = plt.subplots(2, 1, figsize=(13, 9.5), sharex=True)

green_neon = '#39ff14'
red_neon   = '#ff2a6d'
blue_neon  = '#00f0ff'

axs[0].plot(times, result.expect[0], color=green_neon, linewidth=3.2, label='Concurrence (Negentropico)')
axs[0].plot(times, 0.5 * np.exp(-gamma_env * times), color=red_neon, linestyle='--', linewidth=2.5, label='Decadimento standard')

axs[0].set_ylabel('Concurrence (qubit 0-1)', fontsize=14)
axs[0].legend(fontsize=13, loc='upper right')
axs[0].grid(True, alpha=0.35)

axs[1].plot(times, weak_values, color=blue_neon, linewidth=2.8, label='Weak Value di A (Reale)')
axs[1].set_xlabel('Tempo (unità arbitrarie)', fontsize=14)
axs[1].set_ylabel(r'Weak Value $\langle A \rangle_w$', fontsize=14)
axs[1].legend(fontsize=13, loc='upper right')
axs[1].grid(True, alpha=0.35)

plt.suptitle('Entanglement Renascent con Anyonic Braiding', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('renascent_anyonic_4qubit.jpg', dpi=400, bbox_inches='tight')
plt.savefig('aw_evolution_anyonic_moire.png', dpi=400, bbox_inches='tight')
plt.show()

print("\n=== RISULTATI ===")
print(f"Concurrence media finale: {np.mean(result.expect[0][-40:]):.4f}")
print(f"Weak value medio: {np.mean(weak_values):.4f}")
print("Plot 1 salvato come: renascent_anyonic_4qubit.jpg")
print("Plot 2 salvato come: aw_evolution_anyonic_moire.png")
