# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 09:25:19 2026

@author: gregoire.hendrix
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator as PI
import os

factor = 1.5
# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'        : 'serif',
    'font.size'          : 28*factor,
    'axes.labelsize'     : 26*factor,
    'xtick.labelsize'    : 24*factor,
    'ytick.labelsize'    : 24*factor,
    'legend.fontsize'    : 22*factor,
    'figure.dpi'         : 150,
    'axes.grid'          : True,
    'grid.alpha'         : 0.25,
    'grid.linestyle'     : '--',
    'text.usetex'        : True,
})


SAVE_DIR = r"C:\Users\gregoire.hendrix@johncockerill.com\OneDrive - John Cockerill\Documents\Cockerill\Images\1. State of the Art\3. Power blocks\Reference solution"

# ── Data ───────────────────────────────────────────────────────────────────
efficiency = [37.63, 44.09, 45.61, 46.28]
elec_power = [4, 25, 50, 100]

interp = PI(elec_power, efficiency)

x = np.linspace(min(elec_power), max(elec_power), 300)
y = interp(x)

# ── Figure ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(x, y, color='cornflowerblue', label='Interpolated efficiency', linewidth=3)
ax.scatter(elec_power, efficiency, s=100, label='Data points', zorder=5)

ax.set_xlabel('Electrical power output [MW$_{\mathrm{el}}$]')
ax.set_ylabel(r'Turbine efficiency [\%]')
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, 'steam_interpolated_efficiency.pdf'),
            dpi=150, bbox_inches='tight')
plt.show()

# ── Query ──────────────────────────────────────────────────────────────────
x_query = 18.5
y_query = interp(x_query)
print(f"Interpolated efficiency at {x_query} MW = {y_query:.3f} \%")