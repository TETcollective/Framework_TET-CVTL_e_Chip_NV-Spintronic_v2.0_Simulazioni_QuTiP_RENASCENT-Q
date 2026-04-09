 # protocollo_tsvf_negentropico_chip2.py # Simulazione del protocollo TSVF-negentropico per Chip 2.0 # Versione stabile per figura: Entanglement Renascent + Weak Value Evolution

import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

# ====================== PARAMETRI ======================
times = np.linspace(0, 12, 300)
beta = 20.0
gamma_env_base = 0.12 # Define base environmental gamma

# ====================== STATO INIZIALE ======================
# Corrected: Define bell as a 2-qubit Bell state, then tensor with the third qubit
bell = (qt.tensor(qt.basis(2,0), qt.basis(2,0)) + qt.tensor(qt.basis(2,1), qt.basis(2,1))).unit() / np.sqrt(2)
psi0 = qt.tensor(bell, qt.basis(2, 0)) # This makes psi0 a 3-qubit state

# ====================== HAMILTONIANA ======================
H = 1.1 * qt.tensor(qt.sigmax(), qt.sigmax(), qt.qeye(2)) + \
    0.4 * (qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2)) +
           qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2)))

# ====================== DISSIPATORI ======================

# Standard (environmental) dissipators: these are time and state independent.
# So we can define them as a list of Qobj directly.
c_ops_standard_list = [
    np.sqrt(gamma_env_base) * qt.tensor(qt.sigmap(), qt.qeye(2), qt.qeye(2)),
    np.sqrt(gamma_env_base) * qt.tensor(qt.qeye(2), qt.sigmap(), qt.qeye(2))
]

# Pre-compute the A operator for weak value calculation, as it's constant
A_operator_for_weak_value = qt.tensor(qt.sigmax(), qt.sigmax()) # CHANGED: Using σx ⊗ σx for non-zero expectation with Bell state

# Helper function for negentropic gamma calculation
def get_effective_gamma_for_negentropy(t, rho, beta_val, gamma_env_val):
    if rho is None: # For QobjEvo initialization, rho might be None
        return gamma_env_val # Return a default value

    try:
        rho_sub = rho.ptrace([0, 1])
        # Use the pre-computed A_operator_for_weak_value
        weak_val = qt.expect(A_operator_for_weak_value, rho_sub)
        effective_gamma = gamma_env_val * max(0.15, 1.0 - beta_val * abs(weak_val))
        return effective_gamma
    except Exception:
        return gamma_env_val # Fallback in case of errors

# Negentropic collapse operator 1 (returns a single Qobj)
def negentropic_diss_op1(t, rho=None, args=None):
    if rho is None or args is None:
        # Return a zero operator during QobjEvo initialization
        return qt.tensor(qt.sigmap(), qt.qeye(2), qt.qeye(2)) * 0.0

    beta_val = args.get('beta', beta) # Get beta from args, fallback to global beta
    gamma_env_val = args.get('gamma_env_base', gamma_env_base) # Get gamma_env_base from args, fallback to global
    eff_gamma = get_effective_gamma_for_negentropy(t, rho, beta_val, gamma_env_val)
    return np.sqrt(eff_gamma) * qt.tensor(qt.sigmap(), qt.qeye(2), qt.qeye(2))

# Negentropic collapse operator 2 (returns a single Qobj)
def negentropic_diss_op2(t, rho=None, args=None):
    if rho is None or args is None:
        # Return a zero operator during QobjEvo initialization
        return qt.tensor(qt.qeye(2), qt.sigmap(), qt.qeye(2)) * 0.0

    beta_val = args.get('beta', beta)
    gamma_env_val = args.get('gamma_env_base', gamma_env_base)
    eff_gamma = get_effective_gamma_for_negentropy(t, rho, beta_val, gamma_env_val)
    return np.sqrt(eff_gamma) * qt.tensor(qt.qeye(2), qt.sigmap(), qt.qeye(2))

c_ops_negentropic_list = [negentropic_diss_op1, negentropic_diss_op2]

def concurrence_main(t, rho):
    try:
        return qt.concurrence(rho.ptrace([0, 1]))
    except:
        return 0.0

# Define a monitor for weak value as e_ops for easier plotting
def weak_val_monitor_func(t, rho):
    try:
        rho_sub = rho.ptrace([0, 1])
        # Use the pre-computed A_operator_for_weak_value
        return qt.expect(A_operator_for_weak_value, rho_sub)
    except:
        return 0.0

# ====================== SIMULAZIONI ======================
result_std = qt.mesolve(H, psi0, times, c_ops=c_ops_standard_list, # Use the list of Qobj operators
                        e_ops=[concurrence_main],
                        options={'nsteps': 30000, 'atol': 1e-9})

result_neg = qt.mesolve(H, psi0, times, c_ops=c_ops_negentropic_list, # Use the list of functions
                        e_ops=[concurrence_main, weak_val_monitor_func], # Add weak_val_monitor_func
                        args={'beta': beta, 'gamma_env_base': gamma_env_base}, # Pass parameters
                        options={'nsteps': 30000, 'atol': 1e-9, 'store_states': True})

# ====================== WEAK VALUE ======================
# Now that weak_val_monitor_func is an e_op, we can get it directly from .expect
weak_values = result_neg.expect[1] # Assumes it's the second e_op

# Print range of weak_values for debugging
print(f"Weak Values - Min: {np.min(weak_values):.4f}, Max: {np.max(weak_values):.4f}")

# ====================== PLOT SEMPLICE E PULITO ======================
fig, axs = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

axs[0].plot(times, result_std.expect[0], 'r--', linewidth=2.5, label='Concurrence (Standard)')
axs[0].plot(times, result_neg.expect[0], 'g-', linewidth=3.0, label=f'Concurrence (Negentropico, β = {beta})')
axs[0].set_ylabel('Concurrence (qubit 0-1)', fontsize=13)
axs[0].legend(fontsize=13)
axs[0].grid(True, alpha=0.35)

axs[1].plot(times, weak_values, 'b-', linewidth=2.5, label='Weak Value of A (Negentropic)')
axs[1].set_xlabel('Time (arbitrary units)', fontsize=13)
axs[1].set_ylabel('Weak Value of A', fontsize=13)
axs[1].legend(fontsize=13)
axs[1].grid(True, alpha=0.35)

# Dynamically adjust y-limits for weak value plot if there's variation
w_v_min, w_v_max = np.min(weak_values), np.max(weak_values)
if abs(w_v_max - w_v_min) > 1e-9: # If there is significant variation
    buffer = (w_v_max - w_v_min) * 0.1 # Add 10% buffer
    axs[1].set_ylim(w_v_min - buffer, w_v_max + buffer)
else:
    # If values are essentially constant, ensure a small, visible range around the value
    if abs(w_v_min) < 1e-9: # If essentially zero
        axs[1].set_ylim(-0.1, 0.1) # Set a small symmetric range around zero
    else: # If constant non-zero
        axs[1].set_ylim(w_v_min - abs(w_v_min)*0.1, w_v_max + abs(w_v_max)*0.1)

plt.suptitle('Entanglement Renascent - Concurrence and Weak Value Evolution', fontsize=15)
plt.tight_layout()
plt.show()

# ====================== RISULTATI ======================
print("\n=== RISULTATI FINALI ===")
print(f"Concurrence finale (standard)      : {result_std.expect[0][-1]:.4f}")
print(f"Concurrence finale (negentropico)  : {result_neg.expect[0][-1]:.4f}")
guadagno = (result_neg.expect[0][-1] - result_std.expect[0][-1]) / max(result_std.expect[0][-1], 1e-6)
print(f"Guadagno relativo alla fine        : {guadagno:+.1%}")