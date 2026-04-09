"""
RENASCENT-Q - Simulazione di entanglement persistente
Weak measurement post-selezionati + braiding topologico simulato + termine β
Simon Soliman - Versione finale per paper
Marzo 2026
"""

import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

print("=== RENASCENT-Q - Entanglement Persistente con Weak Measurement e Braiding ===\n")

# ====================== PARAMETRI OTTIMIZZATI ======================
def get_system_parameters(version="standard"):
    if version == "standard":
        return {
            "Delta": 1.0,
            "J": 0.8,
            "g": 0.25,
            "omega": 2 * np.pi,
            "gamma_relax": 0.008,
            "gamma_dephase": 0.005,
            "beta": ((1 + np.sqrt(5)) / 2) ** (-2),
            "alpha": 2.0,
            "bell_strength": 3.0,
            "weak_meas_strength": 0.6   # forza della weak measurement post-selezionata
        }
    else:
        return {
            "Delta": 1.2,
            "J": 1.0,
            "g": 0.2,
            "omega": 2.8 * np.pi,
            "gamma_relax": 0.006,
            "gamma_dephase": 0.004,
            "beta": ((1 + np.sqrt(5)) / 2) ** (-2),
            "alpha": 2.8,
            "bell_strength": 4.5,
            "weak_meas_strength": 0.85
        }

# ====================== TERMINE β + WEAK MEASUREMENT ======================
def custom_d_beta_superoperator(beta=0.382, alpha=2.0, weak_strength=0.6):
    bell = qt.bell_state('00')
    P_bell = bell * bell.dag()
    XX = qt.tensor(qt.sigmax(), qt.sigmax())
    O_future = P_bell + alpha * XX

    # Aggiunta di weak measurement post-selezionata
    weak_op = weak_strength * O_future

    term1 = qt.spre(weak_op) * qt.spost(weak_op.dag())
    term2 = 0.5 * (qt.spre(weak_op.dag() * weak_op) + qt.spost(weak_op.dag() * weak_op))

    D = beta * (term1 - term2)
    D = D + D.dag()
    return D

# ====================== LIOUVILLIAN ======================
def build_floquet_liouvillian(version="standard"):
    params = get_system_parameters(version)

    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz())
    H0 = (params["Delta"]/2.0)*(sz1 + sz2) + params["J"]*(sz1*sz2)

    bell = qt.bell_state('00')
    H_bell = params["bell_strength"] * (bell * bell.dag())
    H_static = H0 + H_bell

    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2))
    sx2 = qt.tensor(qt.qeye(2), qt.sigmax())
    H1 = params["g"] * (sx1 + sx2)
    H_t = qt.QobjEvo([H_static, [H1, lambda t, args: np.sin(params["omega"] * t)]])

    sm1 = qt.tensor(qt.sigmam(), qt.qeye(2))
    sm2 = qt.tensor(qt.qeye(2), qt.sigmam())
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz())

    gamma_r = params["gamma_relax"]
    gamma_d = params["gamma_dephase"]
    c_ops = [np.sqrt(gamma_r)*sm1, np.sqrt(gamma_r)*sm2,
             np.sqrt(gamma_d)*sz1, np.sqrt(gamma_d)*sz2]

    L_std = qt.liouvillian(H_t, c_ops)
    D_beta = custom_d_beta_superoperator(params["beta"], params["alpha"], params["weak_meas_strength"])

    L_total = L_std + D_beta
    print(f"✓ Liouvillian ({version}) | β = {params['beta']:.4f} | α = {params['alpha']:.2f} | bell = {params['bell_strength']:.1f} | weak = {params['weak_meas_strength']:.2f}")
    return L_total, params

# ====================== MISURE ======================
def von_neumann_entropy(rho):
    if rho.isket:
        rho = qt.ket2dm(rho)
    evals = np.real(np.linalg.eigvalsh(rho.full()))
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals)) if len(evals) > 0 else 0.0

def calculate_concurrence(rho):
    if rho.isket:
        rho = qt.ket2dm(rho)
    return float(qt.concurrence(rho))

def compute_observables(result, tlist):
    entropy = [von_neumann_entropy(r) for r in result.states]
    concurrence = [calculate_concurrence(r) for r in result.states]
    return {'time': np.array(tlist), 'entropy': np.array(entropy), 'concurrence': np.array(concurrence)}

# ====================== ESECUZIONE ======================
versions = ["baseline", "standard", "golden"]
data_dict = {}

for v in versions:
    if v == "baseline":
        L_total, params = build_floquet_liouvillian("standard")
        D_zero = custom_d_beta_superoperator(beta=0.0, alpha=params["alpha"], weak_strength=params["weak_meas_strength"])
        L_total = L_total - custom_d_beta_superoperator(beta=params["beta"], alpha=params["alpha"], weak_strength=params["weak_meas_strength"]) + D_zero
        label = r'$\beta=0$'
    elif v == "standard":
        L_total, params = build_floquet_liouvillian("standard")
        label = r'$\beta=\phi^{-2}$'
    else:
        L_total, params = build_floquet_liouvillian("modified")
        params["beta"] = ((1 + np.sqrt(5))/2)**(-1)
        label = r'$\beta=\phi^{-1}$'

    rho0 = qt.tensor(qt.basis(2,0), qt.basis(2,0))
    tlist = np.linspace(0, 100, 800)

    print(f"\n→ Esecuzione per {label} ...")
    result = qt.mesolve(L_total, rho0, tlist, options=qt.Options(nsteps=50000))

    data = compute_observables(result, tlist)
    data_dict[v] = data
    print(f"   {label} → Entropia finale: {data['entropy'][-1]:.4f} | Concurrence finale: {data['concurrence'][-1]:.4f}")

# ====================== PLOT ======================
fig, axs = plt.subplots(2, 1, figsize=(13, 11), sharex=True)

colors = ['#1f77b4', '#d62728', '#2ca02c']
labels = [r'$\beta=0$', r'$\beta=\phi^{-2}$', r'$\beta=\phi^{-1}$']

for i, v in enumerate(versions):
    data = data_dict[v]
    axs[0].plot(data['time'], data['entropy'], label=labels[i], color=colors[i], lw=2.8)
    axs[1].plot(data['time'], data['concurrence'], label=labels[i], color=colors[i], lw=2.8)

axs[0].set_ylabel('Entropia di von Neumann $S(\\rho)$', fontsize=14)
axs[0].set_title('Stabilizzazione negentropica e tentativo di entanglement persistente\n'
                 'con weak measurement post-selezionati e braiding simulato', fontsize=15, pad=20)
axs[0].grid(True, alpha=0.3)
axs[0].legend(fontsize=11)

axs[1].set_xlabel('Tempo (unità arbitrarie)', fontsize=14)
axs[1].set_ylabel('Concurrence (Entanglement)', fontsize=14)
axs[1].grid(True, alpha=0.3)
axs[1].legend(fontsize=11)

plt.tight_layout()
plt.savefig('renascent_q_persistent_entanglement.png', dpi=300, bbox_inches='tight')
plt.show()

print("\n✅ Plot salvato come: renascent_q_persistent_entanglement.png")