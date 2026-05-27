# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:00:18 2026
@author: gregoire.hendrix
"""
# sweep_1D.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from full_piping_losses import (
    CSPOptimizerRows, compute_pipe_charges_and_diameters,
    compute_cascaded_temperatures,
    T_HOT, T_AMB, M_DOT_UNIT, DIAM_TABLE, DIST_CC
)

import os
SAVE_DIR = r"C:\Users\gregoire.hendrix@johncockerill.com\OneDrive - John Cockerill\Documents\Cockerill\Images\3. Thermodynamic analysis\1. Piping"

plt.close("all")
plt.rcParams.update({
    'font.size'      : 28,
    'axes.labelsize' : 28,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 24,
})

# =============================================================================
# PARAMETERS
# =============================================================================
CAP      = 5
DIST     = DIST_CC / 2
nb_range = range(5, 151, 1)
T_C      = T_HOT - 273.15
Cp_salt  = 1443.0 + 0.172 * T_C

# =============================================================================
# SWEEP
# =============================================================================
nb_list, length_list, Q_list, dT_list = [], [], [], []

for nb in nb_range:
    print(f"Running nb = {nb} ...")
    opt         = CSPOptimizerRows(DIST, CAP)
    blue_pos, M = opt.generate_blue_towers(nb)
    sol         = opt.find_best_pb(blue_pos, M)

    pipe_info, charge, children = compute_pipe_charges_and_diameters(
        sol["parent"], nb, DIAM_TABLE)

    node_T_out, pipe_results, T_PB = compute_cascaded_temperatures(
        parent     = sol["parent"],
        coords     = sol["coords"],
        pipe_info  = pipe_info,
        children   = children,
        T_hot      = T_HOT,
        T_amb      = T_AMB,
        m_dot_unit = M_DOT_UNIT,
        Cp_salt    = Cp_salt,
    )

    total_Q = sum(info["Q_loss"] for info in pipe_results.values())
    nb_list.append(nb)
    length_list.append(sol["cost"])
    Q_list.append(total_Q / 1e6)
    dT_list.append(T_HOT - T_PB)

# =============================================================================
# SMOOTHING
# =============================================================================
length_smooth = savgol_filter(length_list, window_length=21, polyorder=2)
Q_smooth      = savgol_filter(Q_list,      window_length=21, polyorder=2)
dT_smooth     = savgol_filter(dT_list,     window_length=30, polyorder=2)

# =============================================================================
# FIGURES
# =============================================================================
def save_fig(filename):
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, filename), bbox_inches="tight")
    plt.show()

# --- Figure 1 : Piping length ---
fig, ax = plt.subplots(figsize=(14, 10))
#ax.plot(nb_list, length_list,  color='steelblue', linewidth=1,   alpha=0.2)
ax.plot(nb_list, length_smooth, color='steelblue', linewidth=2.5)
ax.set_xlabel("Number of solar units [-]")
ax.set_ylabel("Total piping length [m]")
ax.grid(True, alpha=0.3)
save_fig("sweep_1D_length.pdf")

# --- Figure 2 : Thermal losses ---
fig, ax = plt.subplots(figsize=(14, 10))
#ax.plot(nb_list, Q_list,   color='firebrick', linewidth=1,   alpha=0.2)
ax.plot(nb_list, Q_smooth, color='firebrick', linewidth=2.5)
ax.set_xlabel("Number of solar units [-]")
ax.set_ylabel("Total heat losses [MW]")
ax.grid(True, alpha=0.3)
save_fig("sweep_1D_Qloss.pdf")

# --- Figure 3 : Temperature drop ---
fig, ax = plt.subplots(figsize=(14, 10))
#ax.plot(nb_list, dT_list,  color='darkorange', linewidth=1,   alpha=0.2)
ax.plot(nb_list, dT_smooth, color='darkorange', linewidth=2.5)
ax.set_xlabel("Number of solar units [-]")
ax.set_ylabel(r"Temperature drop at PB inlet $\Delta T$ [K]")
ax.grid(True, alpha=0.3)
save_fig("sweep_1D_dT.pdf")