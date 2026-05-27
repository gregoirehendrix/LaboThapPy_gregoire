# -*- coding: utf-8 -*-
"""
Steam Rankine — sweep T_salt from 350 to 650°C
Author: Grégoire Hendrix
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from CoolProp.CoolProp import PropsSI
from labothappy.connector.mass_connector import MassConnector
import os

# paste your solve_cycle and compute_fwh functions here
from steam_LPfwh_1open_4closed_HPfwh_2closed import solve_cycle, compute_fwh   # adjust import to your file name

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 28, 'axes.labelsize': 28,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 20,
    'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25, 'grid.linestyle': '--',
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\3. Thermodynamic analysis\2. Power Blocks")
os.makedirs(SAVE_DIR, exist_ok=True)

# --- Fixed parameters ---
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

T_ref = 557.7   # reference point after piping losses

# --- Sweep ---
T_sweep   = np.arange(350, 655, 5)
T_ok, eta_ok   = [], []
T_fail         = []

for T_C in T_sweep:
    print(f"T_salt = {T_C:.0f}°C ... ", end='', flush=True)
    try:
        pump, eco, eva, sh, rh, turb_hp, turb_ic, cond = \
            solve_cycle(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                        T_C + 273.15, P_salt, m_dot_salt_tot, f_rh,
                        CSource, P_low, P_high, P_rh, m_dot_st)

        perf = compute_fwh(pump, eco, eva, sh, rh, turb_hp, turb_ic, cond,
                           eta_pump, eta_turb, P_high, P_low,
                           P_hp_ext, P_fwh5, P_dea, P_lp4, P_lp3, P_lp2, P_lp1,
                           eta_g, pinch_cfwh=pinch_cfwh)

        if perf['eta_el'] is not None and 0.2 < perf['eta_el'] < 0.65:
            T_ok.append(T_C)
            eta_ok.append(perf['eta_el'] * 100)
            print(f"eta = {perf['eta_el']*100:.2f}%")
        else:
            T_fail.append(T_C)
            print("rejected (out of range)")
    except Exception as e:
        T_fail.append(T_C)
        print(f"failed ({e})")

# --- Reference point ---
eta_ref = np.interp(T_ref, T_ok, eta_ok) if T_ok else None

# --- Plot ---
C_MAIN = '#1f77b4'
C_FAIL = '#d62728'
C_REF  = '#2ca02c'

fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(T_ok, eta_ok, '-o', color=C_MAIN, markersize=6,
        linewidth=2.5, label='Converged solutions')

if T_fail:
    ax.scatter(T_fail, [None] * len(T_fail), marker='x', color=C_FAIL,
               s=120, linewidths=2.5, label='Failed to converge', zorder=5)

if T_ok:
    T_thresh = T_ok[0]
    ax.axvline(x=T_thresh, color=C_FAIL, linestyle='--', linewidth=1.5, alpha=0.7)
    ax.text(T_thresh + 3, min(eta_ok) - 0.3,
            fr'Convergence threshold $\approx {T_thresh:.0f}\,°C$',
            color=C_FAIL, fontsize=18)
"""
if eta_ref is not None:
    ax.scatter([T_ref], [eta_ref], marker='*', color=C_REF,
               s=350, zorder=6,
               label=f'Reference point ({T_ref}°C → {eta_ref:.1f}\%)')
    #ax.axvline(x=T_ref, color=C_REF, linestyle=':', linewidth=1.5, alpha=0.7)
"""
ax.set_xlabel(r'Salt inlet temperature $T_{\mathrm{salt,in}}$ [°C]')
ax.set_ylabel(r'Rankine efficiency $\eta$ [%]')
ax.set_xlim(340, 660)
ax.legend(loc='upper left')
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:.0f}%'))

fig.tight_layout()
path = os.path.join(SAVE_DIR, 'fig_steam_TIT.pdf')
fig.savefig(path, bbox_inches='tight')
print(f'\nSaved → {path}')
plt.show()