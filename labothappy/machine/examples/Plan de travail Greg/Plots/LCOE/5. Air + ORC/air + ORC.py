# -*- coding: utf-8 -*-
"""
Created on Wed May 20 09:14:09 2026
CAPEX breakdown + LCOE vs turbine size — Combined cycle air + steam

@author: gregoire.hendrix
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

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

# ── Data (ordered small → large) ─────────────────────────────────────────────
sizes     = ['1 MW', '2 MW', '5 MW']
LCOE_vals = np.array([86.45, 82.84, 81.88])

# CAPEX in M€
solar_field  = np.array([148.110, 141.876, 141.876])
storage      = np.array([99.609,  94.278,  94.278])
power_block  = np.array([136.271,  134.772,  129.003])
piping       = np.array([ 6.110,  5.599,  5.559])
hx           = np.array([ 37.563,  35.354,  35.354])
elec_civil   = np.array([  4.891,   4.646,   4.646])

total_capex = solar_field + storage + power_block + piping + hx + elec_civil

categories = {
    'Solar field':        (solar_field,  '#e9c46a'),
    'Storage':            (storage,      '#457b9d'),
    'Power block':        (power_block,  '#e63946'),
    'Piping':             (piping,       '#adb5bd'),
    'Heat exchangers':    (hx,           '#2a9d8f'),
    'Electrical & civil': (elec_civil,   '#6d6875'),
}

# ── Figure ────────────────────────────────────────────────────────────────────
x     = np.arange(len(sizes))
width = 0.5

fig, ax1 = plt.subplots(figsize=(14, 10))
ax1.grid(axis='x', visible=False)

# Stacked bars
bottoms = np.zeros(len(sizes))
for label, (data, color) in categories.items():
    ax1.bar(x, data, width, bottom=bottoms,
            color=color, label=label, zorder=3)
    bottoms += data

# Total CAPEX annotation on top of each bar
for i, tot in enumerate(total_capex):
    ax1.text(x[i], tot + 6, f'{tot:.0f} M€',
             ha='center', va='bottom', fontsize=18, fontweight='bold')

# LCOE on secondary axis
ax2 = ax1.twinx()
ax2.plot(x, LCOE_vals, color='black', linestyle='--',
         linewidth=2.5, marker='o', markersize=10, zorder=5, label='LCOE')

# Annotate LCOE values above markers
for i, lcoe in enumerate(LCOE_vals):
    ax2.text(x[i] + 0.07, lcoe + 1.0, f'{lcoe:.2f}',
             ha='left', va='bottom', fontsize=18, color='black')

ax1.set_ylabel('CAPEX [M€]')
ax1.set_ylim(0, max(total_capex) * 1.18)
ax1.set_xticks(x)
ax1.set_xticklabels(sizes)

ax2.set_ylabel('LCOE [€/MWh]')
ax2.set_ylim(50, 115)
ax2.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax2.tick_params(labelsize=24)

# Combined legend
bars_h, bars_l = ax1.get_legend_handles_labels()
line_h, line_l = ax2.get_legend_handles_labels()
ax1.legend(bars_h + line_h, bars_l + line_l,
           loc='best', fontsize=16, ncol=1)

fig.tight_layout()
path = os.path.join(SAVE_DIR, "fig_LCOE_turbine_size_air+ORC.pdf")
fig.savefig(path, bbox_inches='tight')
print(f"Saved -> {path}")
plt.show()