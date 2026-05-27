# -*- coding: utf-8 -*-
"""
Created on Fri May 22 13:32:22 2026

@author: gregoire.hendrix
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
import matplotlib.ticker as ticker
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


labels  = ['Steam Rankine\n100 MW$_{el}$', 'Salt Engine\n 1 unit/tower']
solar   = [122.504, 135.750]
storage = [ 75.535, 87.212]
pb      = [ 26.561, 76.064]
hx      = [ 26.339,  0]
piping  = [ 18.237,  0]
ec      = [  4.685,   5.349]

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

lcoe = [60.68, 68.12]
for i, (l, t) in enumerate(zip(lcoe, totals)):
    ax.text(i, -23, f'LCOE = {l:.1f} €/MWh', ha='center', va='top',
            fontsize=20, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                      edgecolor='gray', alpha=0.8))

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('CAPEX [M€]')

ax.set_ylim(-80, 350)
ax.legend(loc='best', framealpha=0.9)
ax.grid(axis='x', alpha=0)
plt.tight_layout()

fname = "capex_breakdown_steam_vs_salt_engine.pdf"
try:
    plt.savefig(os.path.join(SAVE_DIR, fname), bbox_inches='tight')
    print(f"Saved -> {os.path.join(SAVE_DIR, fname)}")
except Exception:
    plt.savefig(fname, bbox_inches='tight')
    print(f"Saved locally as {fname}")
plt.show()