# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:58:13 2026

@author: gregoire.hendrix

T-s diagram of ideal single-stage Brayton cycle

"""

import matplotlib.pyplot as plt
import numpy as np

width = 0.5
factor = 0.85/width

# ── Style (identique) ─────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'        : 'serif',
    'font.size'          : 24*factor,
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

# ── Données (schéma stylisé) ──────────────────────────────────────────────
# Entropie

s = np.linspace(1, 4.8, 200)

# Points imposés
s1, T1 = 1.0, 1.1
s2, T2 = 1.0, 2.2
s3, T3 = 4.8, 3.6
s4, T4 = 4.8, 2.1

# Normalisation
xi = (s - s1) / (s3 - s1)

# Courbure (ajuste si tu veux plus "bombé")
n = 1.8

# Courbes CONSTRAINTES
T_low  = T1 + (T4 - T1) * xi**n   # passe exactement par (1,1.1) et (4.8,2.1)
T_high = T2 + (T3 - T2) * xi**n   # passe exactement par (1,2.2) et (4.8,3.6)


# ── Figure ────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

# Courbes principales
ax.plot(s, T_low, color='cornflowerblue', linewidth=3)
ax.plot(s, T_high, color='cornflowerblue', linewidth=3)

# Lignes verticales (compression / détente)
ax.plot([s1, s2], [T1, T2], color='cornflowerblue', linewidth=3)
ax.plot([s3, s4], [T3, T4], color='cornflowerblue', linewidth=3)

# Points
ax.scatter([s1, s2, s3, s4], [T1, T2, T3, T4],
           color='cornflowerblue', s=120, zorder=5)

# Annotations
ax.text(s1-0.2, T1-0.2, '1')
ax.text(s2-0.2, T2+0.1, '2')
ax.text(s3+0.1, T3, '3')
ax.text(s4+0.1, T4, '4')

ax.text(2.0, 2.9, r'P2 = 10*P1')
ax.text(2.8, 1.5, r'P1')

# Axes
ax.set_xlabel(r'Entropy')
ax.set_ylabel(r'Temperature')

# Nettoyage des ticks pour effet "schéma"
ax.set_xticks([])
ax.set_yticks([])

ax.set_xlim(0, 6)
ax.set_ylim(0.8, 4)

plt.tight_layout()
plt.savefig("brayton_ts.pdf", dpi=150, bbox_inches="tight")
plt.show()