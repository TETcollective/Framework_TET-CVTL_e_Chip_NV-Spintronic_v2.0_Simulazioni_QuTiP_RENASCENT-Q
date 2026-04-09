import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 13,
    'figure.figsize': (15, 11),
    'figure.dpi': 400,
    'axes.grid': True,
    'grid.alpha': 0.3
})

# ========================= PARAMETRI =========================
omega = 2 * np.pi * 1.0
Delta = 1.0
J = 0.8
g = 0.5
gamma = 0.12
gamma_braid = 0.85
beta = (np.sqrt(5) - 1)/2

# Operatori
sx, sz = qt.sigmax(), qt.sigmaz()
I = qt.qeye(2)
sx1 = qt.tensor(sx, I)
sx2 = qt.tensor(I, sx)
sz1 = qt.tensor(sz, I)
sz2 = qt.tensor(I, sz)
sp1 = qt.tensor(qt.sigmap(), I)
sp2 = qt.tensor(I, qt.sigmap())

H0 = 0.5 * Delta * (sz1 + sz2) + J * sz1 * sz2
H1 = g * (sx1 + sx2)

def H_t(t, args=None):
    return H0 + H1 * np.sin(omega * t)

c_ops_std = [np.sqrt(gamma)*sp1, np.sqrt(gamma)*sp2, np.sqrt(gamma/2)*sz1, np.sqrt(gamma/2)*sz2]

# Corrected to be a two-qubit state
phi_plus = (qt.tensor(qt.basis(2,0), qt.basis(2,0)) + qt.tensor(qt.basis(2,1), qt.basis(2,1))).unit()
O_future = phi_plus * phi_plus.dag() + 0.08 * (sx1 * sx2)

def L_retro(rho):
    term = beta * (O_future * rho * O_future.dag() - 0.5*(O_future.dag()*O_future*rho + rho*O_future.dag()*O_future))
    return term + term.dag()

# Corrected to be a two-qubit identity operator
U_braid = qt.tensor(I, I)

def L_braid(rho):
    return gamma_braid * (U_braid * rho * U_braid.dag() - rho)

def vacuum_torque(rho):
    Lz_approx = 0.5 * (sx1 * sx2 + sz1 * sz2)
    return qt.expect(Lz_approx, rho)

# Custom function for manual partial transpose using numpy
def manual_partial_transpose(rho_dm_qobj, subsystem_index):
    if not rho_dm_qobj.isoper:
        raise ValueError("Input must be a density matrix (operator).")
    if rho_dm_qobj.dims[0] != [2, 2] or rho_dm_qobj.dims[1] != [2, 2]:
        raise ValueError("Expected a two-qubit system with dims=[[2,2],[2,2]].")

    # Get the underlying 4x4 data
    data = rho_dm_qobj.full()

    # Reshape to a 4-rank tensor (q1_row, q2_row, q1_col, q2_col)
    tensor_data = data.reshape(2, 2, 2, 2)

    # Perform partial transpose for the specified subsystem
    if subsystem_index == 0: # Transpose first qubit
        # Swap q1_row (axis 0) and q1_col (axis 2)
        pt_data = np.transpose(tensor_data, (2, 1, 0, 3))
    elif subsystem_index == 1: # Transpose second qubit
        # Swap q2_row (axis 1) and q2_col (axis 3)
        pt_data = np.transpose(tensor_data, (0, 3, 2, 1))
    else:
        raise ValueError("Subsystem index must be 0 or 1 for a two-qubit system.")

    # Reshape back to a 4x4 matrix
    pt_data_matrix = pt_data.reshape(4, 4)

    # Create a new Qobj with the same dimensions as original
    return qt.Qobj(pt_data_matrix, dims=rho_dm_qobj.dims)

# Custom logarithmic negativity function for discrete systems
def calculate_log_negativity_discrete(rho_dm):
    # Ensure rho_dm is a proper Qobj with its full structure re-evaluated
    rho_dm = qt.Qobj(rho_dm.full(), dims=[[2, 2], [2, 2]])

    # Use the manual partial transpose function
    rho_pt = manual_partial_transpose(rho_dm, 1)
    evals = rho_pt.eigenenergies()
    # Logarithmic negativity is defined as log2(sum of absolute eigenvalues of partially transposed matrix)
    return np.log2(np.sum([np.abs(x) for x in evals]))

# ========================= SIMULAZIONE =========================
print("Avvio simulazione...")

times = np.linspace(0, 80, 800)
# Corrected to be a two-qubit initial state
rho0 = (qt.tensor(qt.basis(2,0), qt.basis(2,0)) + qt.tensor(qt.basis(2,1), qt.basis(2,1))).unit().proj()

entropies, log_neg, concurrences, torques = [], [], [], []

rho = rho0
dt = times[1] - times[0]

for t in times:
    result = qt.mesolve(H_t, rho, [t, t+dt], c_ops_std)
    rho = result.states[-1]

    # Explicitly set dimensions after mesolve for consistency with two-qubit ops
    rho = qt.Qobj(rho.full(), dims=[[2, 2], [2, 2]])

    rho = rho + dt * L_retro(rho)
    rho = rho + dt * L_braid(rho)
    rho = rho / rho.tr()

    entropies.append(qt.entropy_vn(rho))
    log_neg.append(calculate_log_negativity_discrete(rho)) # Using custom function
    concurrences.append(qt.concurrence(rho))
    torques.append(vacuum_torque(rho))

print("Simulazione completata.")

# ========================= PLOT FINALE RAFFINATO =========================
fig, axs = plt.subplots(2, 2, figsize=(15, 11))

# (a) Entropia von Neumann
axs[0,0].plot(times, entropies, 'r-', lw=2.8, label='Entropia di von Neumann')
axs[0,0].set_ylabel('S(ρ)')
axs[0,0].set_title('(a) Entropia di von Neumann')
axs[0,0].legend()

# (b) Logarithmic Negativity
axs[0,1].plot(times, log_neg, 'c-', lw=2.8, label=r'Logarithmic Negativity $E_{\mathcal{N}}(\rho)$')
axs[0,1].set_ylabel(r'$E_{\mathcal{N}}(\rho)$')
axs[0,1].set_title('(b) Logarithmic Negativity')
axs[0,1].legend()

# (c) Concurrence
axs[1,0].plot(times, concurrences, 'orange', lw=2.8, label='Concurrence')
axs[1,0].set_xlabel('Tempo')
axs[1,0].set_ylabel('Concurrence')
axs[1,0].set_title('(c) Concurrence')
axs[1,0].legend()

# (d) Vacuum Torque
axs[1,1].plot(times, torques, 'magenta', lw=2.8, label=r'Vacuum Torque $\tau_{\rm vac}$')
axs[1,1].set_xlabel('Tempo')
axs[1,1].set_ylabel(r'$\tau_{\rm vac}$')
axs[1,1].set_title('(d) Vacuum Torque')
axs[1,1].legend()

plt.suptitle('Evoluzione della Master Equation Estesa TET-CVTL\n'
             'Retrocausalità-Negentropica + Braiding Topologico + Vacuum Torque',
             fontsize=16, y=0.96)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('master_equation_extended.png', dpi=400, bbox_inches='tight', facecolor='#05080f')
plt.show()

print("Figura salvata come 'master_equation_extended_final.png'")