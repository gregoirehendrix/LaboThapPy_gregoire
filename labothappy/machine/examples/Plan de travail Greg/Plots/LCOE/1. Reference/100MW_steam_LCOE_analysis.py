# -*- coding: utf-8 -*-
"""
LCOE parametric analysis — Steam Rankine reference
SM x Storage sweep — line plot + heatmap (spline interpolation)
Author: Grégoire Hendrix
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import CubicSpline, RectBivariateSpline
import os
from scipy.interpolate import PchipInterpolator



plt.rcParams.update({
    'font.family': 'serif', 'font.size': 28, 'axes.labelsize': 28,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25,
    'grid.linestyle': '--',
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\4. LCOE")
os.makedirs(SAVE_DIR, exist_ok=True)

# ── Data ─────────────────────────────────────────────────────────────────────
SM_vals = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
ST_vals = np.arange(16, 50, 2)   # 16, 18, 20, … 48  (17 values)

LCOE = np.array([
    [71.73, 70.30, 69.13, 68.30, 67.59, 67.10, 66.82, 66.59,
     66.51, 66.62, 66.86, 67.28, 67.69, 68.10, 68.50, 68.90, 69.29],
    [71.20, 69.63, 68.28, 67.13, 66.12, 65.23, 64.37, 63.54,
     62.75, 62.12, 61.65, 61.24, 60.91, 60.68, 60.80, 61.08, 61.48],
    [72.28, 70.64, 69.18, 68.02, 66.94, 65.91, 64.92, 63.96,
     63.11, 62.46, 62.14, 62.07, 62.25, 62.56, 62.97, 63.47, 63.97],
    [72.73, 71.10, 69.51, 68.23, 66.96, 65.85, 64.98, 64.53,
     64.49, 64.70, 65.16, 65.80, 66.43, 67.06, 67.67, 68.28, 68.88],
    [73.36, 71.56, 69.81, 68.37, 67.28, 66.70, 66.66, 67.04,
     67.67, 68.43, 69.19, 69.93, 70.66, 71.39, 72.10, 72.81, 73.51],
    [74.34, 72.32, 70.49, 69.44, 69.20, 69.53, 70.26, 71.15,
     72.04, 72.91, 73.78, 74.63, 75.48, 76.31, 77.13, 77.94, 78.75],
    [75.48, 73.39, 71.99, 71.82, 72.37, 73.31, 74.33, 75.33,
     76.33, 77.31, 78.28, 79.24, 80.18, 81.12, 82.04, 82.95, 83.85],
    [75.82, 74.22, 73.81, 74.49, 75.64, 76.79, 77.92, 79.04,
     80.15, 81.25, 82.33, 83.40, 84.45, 85.50, 86.53, 87.55, 88.55],
])

# ── Global minimum (on raw data) ─────────────────────────────────────────────
min_idx  = np.unravel_index(np.argmin(LCOE), LCOE.shape)
SM_opt   = SM_vals[min_idx[0]]
ST_opt   = ST_vals[min_idx[1]]
LCOE_opt = LCOE[min_idx]
print(f"Global minimum: {LCOE_opt:.2f} EUR/MWh  |  "
      f"SM = {SM_opt}  |  Storage = {ST_opt} MWh/unit")

# ── Figure 1 — Line plot LCOE vs SM (spline) ─────────────────────────────────
storage_curves = [16, 24, 32, 42, 48]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
SM_fine = np.linspace(SM_vals[0], SM_vals[-1], 400)

fig, ax = plt.subplots(figsize=(14, 10))

for st, col in zip(storage_curves, colors):
    idx = np.where(ST_vals == st)[0][0]
    y   = LCOE[:, idx]
    # Remplace CubicSpline par :
    cs = PchipInterpolator(SM_vals, y)
    ax.plot(SM_fine, cs(SM_fine), color=col, lw=2.5,
            label=f'{st} MWh unit$^{{-1}}$')
    ax.scatter(SM_vals, y, color=col, s=60, zorder=5)

ax.set_xlabel('Solar multiple  $SM$  [--]')
ax.set_ylabel(r'LCOE  [€/MWh]')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.legend(title='Storage capacity', loc='best', fontsize=18)
ax.set_xlim(1.5, 5.0)
fig.tight_layout()
path = os.path.join(SAVE_DIR, "fig_LCOE_lineplot_steam.pdf")
fig.savefig(path, bbox_inches='tight')
print(f"Saved -> {path}")
plt.show()

# ── Figure 2 — Heatmap LCOE (spline on fine grid) ────────────────────────────
interp    = RectBivariateSpline(SM_vals, ST_vals, LCOE)
SM_fine2  = np.linspace(SM_vals[0], SM_vals[-1], 300)
ST_fine2  = np.linspace(ST_vals[0], ST_vals[-1], 300)
LCOE_fine = interp(SM_fine2, ST_fine2)

ST_grid_f, SM_grid_f = np.meshgrid(ST_fine2, SM_fine2)

fig, ax = plt.subplots(figsize=(14, 10))
ax.grid(False)

pcm = ax.pcolormesh(ST_grid_f, SM_grid_f, LCOE_fine,
                    cmap='RdYlGn_r', shading='auto',
                    vmin=LCOE.min(), vmax=LCOE.max())

cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
cbar.set_label(r'LCOE  [€/MWh]', fontsize=24)
cbar.ax.tick_params(labelsize=20)

# Contour lines on fine grid
contours = ax.contour(ST_grid_f, SM_grid_f, LCOE_fine,
                      levels=12, colors='black', linewidths=0.8, alpha=0.5)
ax.clabel(contours, inline=True, fontsize=14, fmt='%.0f')

# Optimum marker
ax.scatter(ST_opt, SM_opt, color='white', s=400, zorder=7,
           marker='*', edgecolors='black', linewidths=1.0,
           label=f'Optimum: {LCOE_opt:.1f} €/MWh')

ax.set_xlabel(r'Storage capacity  [MWh/unit]')
ax.set_ylabel('Solar multiple  $SM$  [--]')
ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.legend(loc='best', fontsize=20)
fig.tight_layout()
path = os.path.join(SAVE_DIR, "fig_LCOE_heatmap_steam.pdf")
fig.savefig(path, bbox_inches='tight')
print(f"Saved -> {path}")
plt.show()