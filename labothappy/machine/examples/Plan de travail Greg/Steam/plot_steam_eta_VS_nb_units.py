# -*- coding: utf-8 -*-
"""
Combined sweep : piping losses → Rankine efficiency vs number of CSP units
Author: Grégoire Hendrix
08/05/2026
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from labothappy.connector.mass_connector import MassConnector
from full_piping_losses import (
    CSPOptimizerRows, compute_pipe_charges_and_diameters,
    compute_cascaded_temperatures,
    T_HOT, T_AMB, M_DOT_UNIT, DIAM_TABLE, DIST_CC
)
from steam_LPfwh_1open_4closed_HPfwh_2closed import solve_cycle, compute_fwh
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

# --- Piping parameters ---
CAP  = 5
DIST = DIST_CC / 2
T_C      = T_HOT - 273.15
Cp_salt  = 1443.0 + 0.172 * T_C

# --- Rankine parameters ---
eta_g      = 0.986
eta_pump   = 0.8
eta_turb   = 0.93
eta_hx     = 0.95
pinch_hx   = 5.0
pinch_cond = 5.0
pinch_cfwh = 5.0
P_high   = 160.00e5
P_rh     =  38.41e5
P_low    =   0.095e5
P_hp_ext =  40.50e5
P_fwh5   =  19.11e5
P_dea    =  11.64e5
P_lp4    =   6.08e5
P_lp3    =   3.17e5
P_lp2    =   1.38e5
P_lp1    =   0.44e5
f_rh           = 0.20
m_dot_st       = 1.0
P_salt         = 2e5
m_dot_salt_tot = 9.2

CSource = MassConnector()
CSource.set_properties(fluid='Water', T=20+273.15, P=3e5, m_dot=100.0)

N_REF = 82

# --- Combined sweep ---
nb_range = range(5, 151, 1)
nb_list, T_PB_list, eta_list = [], [], []
nb_fail = []

for nb in nb_range:
    print(f"N = {nb} ...", end=' ', flush=True)

    # Step 1 : piping → T_PB
    opt         = CSPOptimizerRows(DIST, CAP)
    blue_pos, M = opt.generate_blue_towers(nb)
    sol         = opt.find_best_pb(blue_pos, M)
    pipe_info, charge, children = compute_pipe_charges_and_diameters(
        sol["parent"], nb, DIAM_TABLE)
    _, pipe_results, T_PB = compute_cascaded_temperatures(
        parent=sol["parent"], coords=sol["coords"], pipe_info=pipe_info,
        children=children, T_hot=T_HOT, T_amb=T_AMB,
        m_dot_unit=M_DOT_UNIT, Cp_salt=Cp_salt)

    T_PB_C = T_PB - 273.15
    print(f"T_PB = {T_PB_C:.1f}°C ...", end=' ', flush=True)

    # Step 2 : Rankine cycle at T_PB
    try:
        pump, eco, eva, sh, rh, turb_hp, turb_ic, cond = \
            solve_cycle(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                        T_PB, P_salt, m_dot_salt_tot, f_rh,
                        CSource, P_low, P_high, P_rh, m_dot_st)

        perf = compute_fwh(pump, eco, eva, sh, rh, turb_hp, turb_ic, cond,
                           eta_pump, eta_turb, P_high, P_low,
                           P_hp_ext, P_fwh5, P_dea, P_lp4, P_lp3, P_lp2, P_lp1,
                           eta_g, pinch_cfwh=pinch_cfwh)

        if perf['eta_el'] is not None and 0.2 < perf['eta_el'] < 0.65:
            nb_list.append(nb)
            T_PB_list.append(T_PB_C)
            eta_list.append(perf['eta_el'] * 100)
            print(f"eta = {perf['eta_el']*100:.2f}%")
        else:
            nb_fail.append(nb)
            print("rejected")
    except Exception as e:
        nb_fail.append(nb)
        print(f"failed ({e})")

# --- Reference point ---
eta_ref = np.interp(N_REF, nb_list, eta_list) if N_REF in nb_list or (nb_list[0] < N_REF < nb_list[-1]) else None

# --- Plot ---
C_MAIN = '#1f77b4'
C_FAIL = '#d62728'
C_REF  = '#2ca02c'

fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(nb_list, eta_list, '-o', color=C_MAIN, markersize=5,
        linewidth=2.5, label='Converged solutions')

if nb_fail:
    ax.scatter(nb_fail, [None] * len(nb_fail), marker='x', color=C_FAIL,
               s=120, linewidths=2.5, label='Failed to converge', zorder=5)

if eta_ref is not None:
    ax.scatter([N_REF], [eta_ref], marker='*', color=C_REF, s=350, zorder=6,
               label=f'Reference point (N = {N_REF} → {eta_ref:.1f}%)')
    ax.axvline(x=N_REF, color=C_REF, linestyle=':', linewidth=1.5, alpha=0.7)

ax.set_xlabel(r'Number of solar units $N$ [-]')
ax.set_ylabel(r'Rankine efficiency $\eta$ [%]')
ax.set_xlim(0, 155)
ax.legend(loc='upper right')
ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:.0f}%'))

fig.tight_layout()
path = os.path.join(SAVE_DIR, 'fig_steam_vs_N.pdf')
fig.savefig(path, bbox_inches='tight')
print(f'\nSaved → {path}')
plt.show()
