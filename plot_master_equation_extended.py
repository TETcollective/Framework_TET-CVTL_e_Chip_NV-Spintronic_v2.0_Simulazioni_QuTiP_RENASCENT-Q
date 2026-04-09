import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import matplotlib.patheffects as path_effects

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'figure.figsize': (16, 11),
    'figure.dpi': 400,
})

fig, ax = plt.subplots()
ax.set_xlim(0, 16)
ax.set_ylim(0, 11.5)
ax.axis('off')

# Sfondo
ax.add_patch(Rectangle((0, 0), 16, 11.5, facecolor='#05080f', alpha=0.98))

# Titolo principale
title = ax.text(8, 10.85, 'Master Equation Estesa',
                ha='center', va='center', fontsize=19.5, fontweight='bold', color='white')
title.set_path_effects([path_effects.withStroke(linewidth=4, foreground='#00bfff')])

# Sottotitolo
ax.text(8, 10.25, 'con Protezione Topologica e Termine Retrocausale-Negentropico',
        ha='center', va='center', fontsize=13, color='#88ccff')

# Logo TET-CVTL
ax.text(15, 10.9, 'TET-CVTL', fontsize=15, fontweight='bold', color='#00ffff',
        ha='right', va='center')

# Box equazione principale
eq_box = FancyBboxPatch((2.5, 6.8), 11, 2.9, boxstyle="round,pad=1.0",
                        facecolor='#0f1a2e', edgecolor='#00ccff', linewidth=3)
ax.add_patch(eq_box)

ax.text(8, 8.35, r'$\frac{d\rho}{dt} = -i [H(t),\rho] + \mathcal{L}_{\rm std} + \mathcal{L}_{\rm retro} + \mathcal{L}_{\rm braid}$',
        ha='center', va='center', fontsize=17, color='white')

# === I TRE TERMI ===

# 1. L_std
ax.add_patch(Rectangle((2.8, 3.8), 3.1, 2.7, facecolor='#ff444420', edgecolor='#ff6666', linewidth=2))
ax.text(4.35, 5.75, r'$\mathcal{L}_{\rm std}(\rho)$', ha='center', va='center', fontsize=14, fontweight='bold', color='#ff6666')
ax.text(4.35, 4.6, 'Decoerenza termica\n+ SAW modulation', ha='center', va='center', fontsize=10.2, color='#ffaaaa', linespacing=1.5)

# 2. L_retro
ax.add_patch(Rectangle((6.4, 3.8), 3.1, 2.7, facecolor='#44ff4420', edgecolor='#44ff88', linewidth=2))
ax.text(7.95, 5.75, r'$\mathcal{L}_{\rm retro}(\rho)$', ha='center', va='center', fontsize=14, fontweight='bold', color='#44ff88')
ax.text(7.95, 4.7, r'$\mathcal{D}_\beta[\rho]\quad \beta = \phi^{-2} \approx 0.382$',
        ha='center', va='center', fontsize=12, color='#88ffbb', fontweight='bold')
ax.text(7.95, 4.1, '(Two-State Vector Formalism)', ha='center', va='center', fontsize=9.8, color='#88ffbb')

# 3. L_braid
ax.add_patch(Rectangle((10.0, 3.8), 3.1, 2.7, facecolor='#ffaa4420', edgecolor='#ffcc66', linewidth=2))
ax.text(11.55, 5.75, r'$\mathcal{L}_{\rm braid}(\rho)$', ha='center', va='center', fontsize=14, fontweight='bold', color='#ffcc66')
ax.text(11.55, 4.7, r'$\gamma_{\rm braid} \left( U_{\rm braid}\,\rho\,U^\dagger - \rho \right)$',
        ha='center', va='center', fontsize=11, color='#ffdd99')
ax.text(11.55, 4.1, r'$U_{\rm braid} = \sigma_1^{n_1} \sigma_2^{n_2} F^{m}$',
        ha='center', va='center', fontsize=10.8, color='#ffdd99', fontweight='bold')

# Freccia centrale
ax.annotate('', xy=(13.8, 7.6), xytext=(2.7, 7.6),
            arrowprops=dict(arrowstyle='->', lw=4, color='#00ddff', alpha=0.9))

# === LEGENDA REGIMI ===
reg_box = FancyBboxPatch((2.0, 0.9), 11.5, 2.4, boxstyle="round,pad=0.8",
                         facecolor='#1a2538', edgecolor='#6688cc', linewidth=2)
ax.add_patch(reg_box)

ax.text(7.75, 3.0, r'Regimi di Applicazione di $\gamma_{\rm braid}$',
        ha='center', va='center', fontsize=14.5, fontweight='bold', color='#aaccff')

regimi = [
    (3.6, 2.05, 'Adiabatico lento', r'$\gamma_{\rm braid} \ll \omega_{\rm SAW}$',
     'Massima fedeltà topologica\nMinimo leakage', '#88ffaa'),
    (7.75, 2.05, 'Floquet risonante', r'$\gamma_{\rm braid} \approx k \cdot \omega_{\rm SAW}$',
     'Effective Hamiltonians topologici\nIdeale per simulazioni QuTiP', '#ffcc88'),
    (11.9, 2.05, 'Impulsivo (bang-bang)', r'$\gamma_{\rm braid} \gg \omega_{\rm SAW}$',
     'Correzione rapida errori\nEfficace contro decoerenza termica', '#ff8888')
]

for x, y, title, cond, desc, color in regimi:
    ax.text(x, y+0.32, title, ha='center', va='center', fontsize=12.2, fontweight='bold', color=color)
    ax.text(x, y-0.1, cond, ha='center', va='center', fontsize=10.2, color='#dddddd')
    ax.text(x, y-0.75, desc, ha='center', va='center', fontsize=9.2, color='#bbbbbb', linespacing=1.35)

# Footer
plt.figtext(0.5, 0.07,
            "Framework TET-CVTL  •  Protocollo RENASCENT-Q  •  Master Equation Estesa con Braiding Topologico",
            ha='center', fontsize=11, color='#777777')

plt.tight_layout(rect=[0, 0.12, 1, 0.98])
plt.savefig('master_equation_extended.png', dpi=400, bbox_inches='tight', facecolor='#05080f')
plt.show()