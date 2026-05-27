# -*- coding: utf-8 -*-
"""
Created on Sat May 16 11:23:46 2026

@author: gregoire.hendrix
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
import os

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 28, 'axes.labelsize': 28,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25,
    'grid.linestyle': '--',
})

SAVE_DIR = r"C:\Users\gregoire.hendrix@johncockerill.com\OneDrive - John Cockerill\Documents\Cockerill\Images\4. LCOE"

# ── Data ───────────────────────────────────────────────────────────────────
TIT = np.array([600, 650, 700, 750, 800, 850, 900, 950, 1000])
LCOE_ref = 60.68  # Steam Rankine reference

powers = [0.2, 0.5, 1, 2, 5]  # MW

LCOE_data = {
    0.2: np.array([128.77, 124.18, 120.54, 117.91, 115.67, 114.03, 113.10, 112.31, 111.71]),
    0.5: np.array([117.36, 112.12, 108.27, 105.35, 103.05, 101.14,  99.74,  98.73,  98.10]),
    1:   np.array([112.57, 107.08, 102.96,  99.84,  97.36,  95.53,  94.05,  92.96,  92.25]),
    2:   np.array([107.89, 102.02,  97.80,  94.51,  91.94,  90.01,  88.57,  87.49,  86.76]),
    5:   np.array([101.91,  96.28,  92.20,  89.05,  86.58,  84.71,  83.24,  82.12,  81.37]),
}

# ── Colours (project palette) ─────────────────────────────────────────────
colors = ['#e9c46a', '#457b9d', '#e63946', '#adb5bd', '#2a9d8f']

# ── Figure ─────────────────────────────────────────────────────────────────
TIT_fine = np.linspace(600, 1000, 500)

fig, ax = plt.subplots(figsize=(14, 10))

for (power, lcoe), color in zip(LCOE_data.items(), colors):
    interp = PchipInterpolator(TIT, lcoe)
    label = f'{power} MW$_{{el}}$'
    ax.plot(TIT_fine, interp(TIT_fine), color=color, linewidth=2.5, label=label)
    ax.scatter(TIT, lcoe, color=color, zorder=5, s=60)

ax.axhline(LCOE_ref, color='#6d6875', linewidth=2, linestyle='--',
           label=f'Steam Rankine reference ({LCOE_ref:.2f} €/MWh)')

ax.set_xlabel('Turbine inlet temperature (°C)')
ax.set_ylabel('Optimal LCOE (€/MWh)')
ax.set_xlim(600, 1000)
ax.set_xticks(TIT)
ax.legend(loc='lower left', title='Power block size')

plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, 'brayton_air_LCOE_TIT.pdf'), format='pdf')
plt.show()

#Figure 2; capex breakdown
labels  = ['Steam Rankine\n100 MW$_{el}$', 'Air (Recup + IC + RH)\n20 x 5 MW$_{el}$']
solar   = [122.504, 145.938]
storage = [ 75.535, 101.552]
pb      = [ 26.561, 125.402]
hx      = [ 26.339,  20.707]
piping  = [ 18.237,  0]
ec      = [  4.685,   4.926]

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

lcoe = [60.68, 81.83]
for i, (l, t) in enumerate(zip(lcoe, totals)):
    ax.text(i, -10, f'LCOE = {l:.1f} €/MWh', ha='center', va='top',
            fontsize=20, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                      edgecolor='gray', alpha=0.8))

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('CAPEX [M€]')

ax.set_ylim(-40, 500)
ax.legend(loc='upper left', framealpha=0.9)
ax.grid(axis='x', alpha=0)
plt.tight_layout()

fname = "capex_breakdown_steam_vs_air.pdf"
try:
    plt.savefig(os.path.join(SAVE_DIR, fname), bbox_inches='tight')
    print(f"Saved -> {os.path.join(SAVE_DIR, fname)}")
except Exception:
    plt.savefig(fname, bbox_inches='tight')
    print(f"Saved locally as {fname}")
plt.show()