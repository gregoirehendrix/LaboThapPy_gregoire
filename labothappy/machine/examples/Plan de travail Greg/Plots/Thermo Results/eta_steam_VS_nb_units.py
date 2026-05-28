# -*- coding: utf-8 -*-
"""
Rankine efficiency vs number of CSP units (combined piping + cycle)
Author: Grégoire Hendrix
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import savgol_filter
import os

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 28, 'axes.labelsize': 28,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 20,
    'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25, 'grid.linestyle': '--',
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\3. Thermodynamic analysis\2. Power Blocks")
os.makedirs(SAVE_DIR, exist_ok=True)

nb_list = list(range(5, 151))
eta_list = [
    45.70, 45.69, 45.69, 45.66, 45.66, 45.67, 45.67, 45.67, 45.67, 45.67,
    45.64, 45.64, 45.66, 45.66, 45.66, 45.66, 45.65, 45.65, 45.65, 45.62,
    45.62, 45.63, 45.63, 45.63, 45.63, 45.63, 45.63, 45.63, 45.63, 45.62,
    45.62, 45.62, 45.62, 45.63, 45.63, 45.64, 45.63, 45.63, 45.63, 45.63,
    45.63, 45.63, 45.62, 45.62, 45.62, 45.62, 45.62, 45.62, 45.62, 45.62,
    45.62, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61,
    45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61,
    45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.61, 45.60, 45.60,
    45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60,
    45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60,
    45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60, 45.60,
    45.59, 45.59, 45.59, 45.59, 45.59, 45.59, 45.59, 45.59, 45.58, 45.58,
    45.58, 45.58, 45.59, 45.59, 45.59, 45.59, 45.59, 45.59, 45.59, 45.59,
    45.59, 45.59, 45.59, 45.59, 45.59, 45.59, 45.58, 45.58, 45.57, 45.57,
    45.59, 45.59, 45.59, 45.59, 45.59, 45.59,
]

N_REF   = 82
eta_ref = eta_list[N_REF - 5]

eta_smooth = savgol_filter(eta_list, window_length=21, polyorder=2)

C_MAIN = '#1f77b4'
C_REF  = '#ff7f0e'

fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(nb_list, eta_smooth, color=C_MAIN, lw=2.5, ms=8,
        markevery=25, label='Rankine efficiency')
ax.scatter([N_REF], 45.605, marker='o', color=C_REF, s=150,
           edgecolors=C_REF, linewidths=2.5, zorder=6,
           label=f'Reference point ($N$ = {N_REF}, $\\eta$ = {eta_ref:.2f}%)')
ax.set_xlabel(r'Number of solar units $N$ [-]')
ax.set_ylabel(r'Rankine efficiency $\eta$ [%]')
ax.set_xlim(0, 155)
ax.set_ylim(45.5, 45.8)
ax.legend(loc='upper right')
ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:.2f}'))

fig.tight_layout()
path = os.path.join(SAVE_DIR, 'fig_steam_vs_N.pdf')
fig.savefig(path, bbox_inches='tight')
print(f'Saved → {path}')
plt.show()