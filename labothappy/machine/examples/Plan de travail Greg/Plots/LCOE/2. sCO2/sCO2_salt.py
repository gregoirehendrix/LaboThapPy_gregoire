# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:34:09 2026
@author: gregoire.hendrix

sCO2 recompression at 565 C - Solar Salt
Figures for section 4.3:
  1. CAPEX breakdown - Steam Rankine vs sCO2 recompression
  2. LCOE heatmap   - SM x Storage sweep
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from scipy.interpolate import RectBivariateSpline
import os

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 28,
    'axes.labelsize': 28,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 20,
    'figure.dpi': 150,
    'axes.grid': True,
    'grid.alpha': 0.25,
    'grid.linestyle': '--',
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\4. LCOE")
os.makedirs(SAVE_DIR, exist_ok=True)

# ==============================================================================
# Figure 1 - CAPEX breakdown
# ==============================================================================
labels  = ['Steam Rankine\n100 MW$_{el}$', 'sCO$_2$ Recompression (salt)\n20 x 5 MW$_{el}$']
solar   = [122.504, 151.078]
storage = [ 75.535, 117.272]
pb      = [ 26.561, 531.230]
hx      = [ 26.339,  56.278]
piping  = [ 18.237,  12.180]
ec      = [  4.685,   7.034]

x      = np.arange(len(labels))
width  = 0.45
colors = ['#e9c46a', '#457b9d', '#e63946', '#adb5bd', '#2a9d8f', '#6d6875']

fig, ax = plt.subplots(figsize=(14, 10))
bottoms = np.zeros(len(labels))

for vals, color, name in zip(
    [solar, storage, pb, hx, piping, ec],
    colors,
    ['Solar field', 'Storage', 'Power block', 'Heat exchanger', 'Piping', 'Electrical & Civil']
):
    vals = np.array(vals)
    ax.bar(x, vals, width, bottom=bottoms, color=color, label=name,
           edgecolor='white', linewidth=0.8)
    bottoms += vals

totals = [sum(v) for v in zip(solar, storage, pb, hx, piping, ec)]
for i, total in enumerate(totals):
    ax.text(i, total + 8, f'{total:.0f} M€', ha='center', va='bottom',
            fontsize=22, fontweight='bold')

lcoe = [60.68, 156.31]
for i, (l, t) in enumerate(zip(lcoe, totals)):
    ax.text(i, -23, f'LCOE = {l:.1f} €/MWh', ha='center', va='top',
            fontsize=20, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                      edgecolor='gray', alpha=0.8))

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('CAPEX [M€]')

ax.set_ylim(-80, 980)
ax.legend(loc='upper left', framealpha=0.9)
ax.grid(axis='x', alpha=0)
plt.tight_layout()

fname = "capex_breakdown_steam_vs_sCO2_salt.pdf"
try:
    plt.savefig(os.path.join(SAVE_DIR, fname), bbox_inches='tight')
    print(f"Saved -> {os.path.join(SAVE_DIR, fname)}")
except Exception:
    plt.savefig(fname, bbox_inches='tight')
    print(f"Saved locally as {fname}")
plt.show()

# ==============================================================================
# Figure 2 - LCOE heatmap (SM x Storage)
# ==============================================================================
SM_vals = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
ST_vals = np.arange(16, 56, 2)   # 16 to 54, 20 values

LCOE = np.array([
    [205.00, 200.23, 196.30, 193.10, 190.32, 188.31, 186.75, 185.47,
     184.55, 184.27, 184.54, 184.88, 185.41, 185.93, 186.45, 186.96, 187.47, 187.97, 188.46, 188.95],
    [198.02, 192.83, 188.40, 184.38, 180.69, 177.44, 174.34, 171.34,
     168.47, 166.24, 164.28, 162.67, 161.21, 160.10, 159.75, 159.94, 160.48, 161.01, 161.53, 162.04],
    [191.78, 186.47, 181.88, 177.82, 174.07, 170.46, 167.00, 163.73,
     160.77, 158.54, 157.23, 156.52, 156.31, 156.56, 157.09, 157.73, 158.36, 158.99, 159.60, 160.20],
    [188.61, 183.45, 178.52, 174.17, 169.97, 166.12, 163.01, 161.08,
     160.16, 159.83, 160.18, 160.98, 161.77, 162.55, 163.32, 164.08, 164.83, 165.56, 166.29, 167.02],
    [185.19, 179.53, 174.26, 169.54, 165.81, 163.70, 162.76, 162.87,
     163.62, 164.58, 165.53, 166.47, 167.39, 168.31, 169.22, 170.11, 170.99, 171.86, 172.72, 173.57],
    [182.59, 176.58, 171.13, 167.52, 165.96, 165.70, 166.47, 167.59,
     168.70, 169.80, 170.88, 171.95, 173.01, 174.06, 175.09, 176.10, 177.11, 178.10, 179.08, 180.05],
    [179.67, 173.60, 169.65, 168.35, 168.68, 169.94, 171.22, 172.49,
     173.74, 174.98, 176.20, 177.41, 178.61, 179.79, 180.95, 182.10, 183.24, 184.36, 185.47, 186.56],
    [177.74, 173.03, 171.21, 171.69, 173.14, 174.58, 176.01, 177.41,
     178.80, 180.18, 181.54, 182.88, 184.20, 185.51, 186.80, 188.08, 189.34, 190.58, 191.81, 193.02],
])

min_idx  = np.unravel_index(np.argmin(LCOE), LCOE.shape)
SM_opt   = SM_vals[min_idx[0]]
ST_opt   = ST_vals[min_idx[1]]
LCOE_opt = LCOE[min_idx]
print(f"Optimum: {LCOE_opt:.2f} EUR/MWh | SM = {SM_opt} | Storage = {ST_opt} MWh/unit")

interp    = RectBivariateSpline(SM_vals, ST_vals, LCOE)
SM_fine   = np.linspace(SM_vals[0], SM_vals[-1], 300)
ST_fine   = np.linspace(ST_vals[0], ST_vals[-1], 300)
LCOE_fine = interp(SM_fine, ST_fine)

ST_grid, SM_grid = np.meshgrid(ST_fine, SM_fine)

fig, ax = plt.subplots(figsize=(14, 10))
ax.grid(False)

pcm = ax.pcolormesh(ST_grid, SM_grid, LCOE_fine,
                    cmap='RdYlGn_r', shading='auto',
                    vmin=LCOE.min(), vmax=LCOE.max())

cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
cbar.set_label(r'LCOE  [€/MWh]', fontsize=24)
cbar.ax.tick_params(labelsize=20)

contours = ax.contour(ST_grid, SM_grid, LCOE_fine,
                      levels=12, colors='black', linewidths=0.8, alpha=0.5)
ax.clabel(contours, inline=True, fontsize=14, fmt='%.0f')

ax.scatter(ST_opt, SM_opt, color='white', s=400, zorder=7,
           marker='*', edgecolors='black', linewidths=1.0,
           label=f'Optimum: {LCOE_opt:.1f} €/MWh\n'
                 f'(SM = {SM_opt}, {ST_opt} MWh/unit)')

ax.set_xlabel(r'Storage capacity  [MWh/unit]')
ax.set_ylabel('Solar multiple  $SM$  [--]')
ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.legend(loc='upper left', fontsize=20)
fig.tight_layout()

fname = "fig_LCOE_heatmap_sCO2_565.pdf"
try:
    fig.savefig(os.path.join(SAVE_DIR, fname), bbox_inches='tight')
    print(f"Saved -> {os.path.join(SAVE_DIR, fname)}")
except Exception:
    fig.savefig(fname, bbox_inches='tight')
    print(f"Saved locally as {fname}")
plt.show()