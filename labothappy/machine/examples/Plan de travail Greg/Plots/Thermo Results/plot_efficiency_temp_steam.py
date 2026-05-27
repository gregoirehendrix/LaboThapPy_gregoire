# -*- coding: utf-8 -*-
"""
Steam Rankine cycle — thermal efficiency vs salt inlet temperature
Author: Grégoire Hendrix
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import numpy as np
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

T_ok = [375.62, 377.80, 390.79, 407.42, 408.86, 411.92, 412.72, 424.78,
        442.21, 443.85, 447.76, 458.85, 475.65, 478.98, 480.36, 482.89,
        493.02, 509.15, 512.63, 516.95, 518.09, 527.24, 542.72, 546.33,
        553.34, 553.59, 561.51, 576.33, 580.07, 588.64, 590.28,
        623.97, 627.01, 663.77]
eta_ok = [38.06, 38.24, 39.11, 39.89, 39.95, 40.06, 40.09, 40.47,
          41.05, 41.13, 41.34, 42.06, 43.44, 43.71, 43.81, 43.98,
          44.76, 45.25, 45.33, 45.42, 45.44, 45.63, 45.95, 46.03,
          46.17, 46.17, 46.33, 46.62, 46.69, 46.86, 46.89,
          47.52, 47.57, 48.23]
T_fail = [309.56, 312.18, 323.22, 342.51, 343.00, 345.28,
          356.93, 371.11, 378.54, 445.41, 595.83]

T_ref = 557.7
eta_ref = np.interp(T_ref, T_ok, eta_ok)

C_MAIN = '#1f77b4'
C_FAIL = '#d62728'
C_REF  = '#2ca02c'

fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(T_ok, eta_ok, '-o', color=C_MAIN, markersize=6,
        linewidth=2.5, label='Converged solutions')
ax.scatter(T_fail, [None] * len(T_fail), marker='x', color=C_FAIL,
           s=120, linewidths=2.5, label='Failed to converge', zorder=5)

ax.axvline(x=375, color=C_FAIL, linestyle='--', linewidth=1.5, alpha=0.7)
ax.text(378, 34.6, r'Convergence threshold $\approx 375\,°C$',
        color=C_FAIL, fontsize=18)

ax.scatter([T_ref], [eta_ref], marker='*', color=C_REF,
           s=350, zorder=6, label=f'Reference point (565°C → {eta_ref:.2f}\%)')
#ax.axvline(x=T_ref, color=C_REF, linestyle=':', linewidth=1.5, alpha=0.7)

ax.set_xlabel(r'Salt inlet temperature $T_{\mathrm{salt,in}}$ [°C]')
ax.set_ylabel(r'Rankine efficiency $\eta$ [%]')
ax.set_xlim(290, 690)
ax.set_ylim(34, 50)
ax.legend(loc='upper left')
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:.0f}%'))

fig.tight_layout()
path = os.path.join(SAVE_DIR, 'fig_steam_TIT.pdf')
fig.savefig(path, bbox_inches='tight')
print(f'Saved → {path}')
plt.show()