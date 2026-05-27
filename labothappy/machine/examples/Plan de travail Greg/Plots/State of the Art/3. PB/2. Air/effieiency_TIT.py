# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:24:59 2026

@author: gregoire.hendrix

Efficiency vs TIT (publication style)
"""

import matplotlib.pyplot as plt
import numpy as np

width = 0.5
factor = 0.85/width

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'        : 'serif',
    'font.size'          : 28*factor,
    'axes.labelsize'     : 28*factor,
    'xtick.labelsize'    : 24*factor,
    'ytick.labelsize'    : 24*factor,
    'legend.fontsize'    : 22*factor,
    'figure.dpi'         : 150,
    'axes.grid'          : True,
    'grid.alpha'         : 0.25,
    'grid.linestyle'     : '--',
    'text.usetex'        : True,
})

# ── Données (linéaire comme sur le schéma) ────────────────────────────────
TIT = np.linspace(800, 1400, 100)          # °C
efficiency = 30 + 0.025 * (TIT - 800)      # % (croissance simple)

# ── Figure ───────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(TIT, efficiency, color='cornflowerblue', linewidth=4.5)

# Labels (identiques au schéma)
ax.set_xlabel(r'TIT ($^\circ$C)')
ax.set_ylabel(r'Cycle efficiency (\%)')

# Limites pour cadrage visuel similaire
ax.set_xlim(750, 1450)
ax.set_ylim(28, 46)

# Pas de ticks visibles (style schématique)
ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()
plt.savefig("efficiency_vs_TIT.pdf", dpi=150, bbox_inches="tight")
plt.show()
