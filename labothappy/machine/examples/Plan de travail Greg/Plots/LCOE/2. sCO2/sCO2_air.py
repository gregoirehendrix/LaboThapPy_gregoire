# -*- coding: utf-8 -*-
"""
Created on Sat May 16 11:15:42 2026

@author: gregoire.hendrix
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
import os

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 28,
    'axes.labelsize': 28,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 18,
    'figure.dpi': 150,
    'axes.grid': True,
    'grid.alpha': 0.25,
    'grid.linestyle': '--',
})

SAVE_DIR = r"C:\Users\gregoire.hendrix@johncockerill.com\OneDrive - John Cockerill\Documents\Cockerill\Images\4. LCOE"

# ── Data ───────────────────────────────────────────────────────────────────
TIT = np.array([600, 650, 700, 750, 800, 850, 900, 950, 1000])
LCOE = np.array([157.56, 154.79, 152.68, 150.96, 149.68, 148.66, 147.97, 147.43, 147.16])

LCOE_ref = 60.68  # Steam Rankine reference

# ── Interpolation ──────────────────────────────────────────────────────────
TIT_fine = np.linspace(600, 1000, 500)
interp = PchipInterpolator(TIT, LCOE)
LCOE_fine = interp(TIT_fine)

# ── Figure ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(TIT_fine, LCOE_fine, color='steelblue', linewidth=2.5,
        label='sCO$_2$ recompression — air HTF')
ax.scatter(TIT, LCOE, color='steelblue', zorder=5, s=80)

ax.axhline(LCOE_ref, color='firebrick', linewidth=2, linestyle='--',
           label=f'Steam Rankine reference ({LCOE_ref:.2f} €/MWh)')

ax.set_xlabel('Turbine inlet temperature (°C)')
ax.set_ylabel('Optimal LCOE (€/MWh)')
ax.set_xlim(600, 1000)
ax.set_xticks(TIT)
ax.legend(loc='center right')

plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, 'sCO2_air_LCOE_TIT.pdf'), format='pdf')
plt.show()


#Figure 2; capex breakdown
labels  = ['Steam Rankine\n100 MW$_{el}$', 'sCO$_2$ Recompression (air)\n20 x 5 MW$_{el}$']
solar   = [122.504, 143.551]
storage = [ 75.535, 109.278]
pb      = [ 26.561, 531.230]
hx      = [ 26.339,  36.813]
piping  = [ 18.237,  9.157]
ec      = [  4.685,   5.413]

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

lcoe = [60.68, 148.59]
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

fname = "capex_breakdown_steam_vs_sCO2_air.pdf"
try:
    plt.savefig(os.path.join(SAVE_DIR, fname), bbox_inches='tight')
    print(f"Saved -> {os.path.join(SAVE_DIR, fname)}")
except Exception:
    plt.savefig(fname, bbox_inches='tight')
    print(f"Saved locally as {fname}")
plt.show()