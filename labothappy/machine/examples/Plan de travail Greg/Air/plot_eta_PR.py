# -*- coding: utf-8 -*-
"""
Air Brayton (Recuperated + IC + RH) — efficiency vs T_hot for 5 PR values
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

from labothappy.connector.mass_connector       import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector

from air_Brayton import (brayton_recuperated_rh_ic, compute_cycle_performance)

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 28, 'axes.labelsize': 28,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 20,
    'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25, 'grid.linestyle': '--',
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\3. Thermodynamic analysis\2. Power Blocks")
os.makedirs(SAVE_DIR, exist_ok=True)

# ── PR values and corresponding power outputs ─────────────────────────────────
PR_configs = [
    (3.0000, 0.2),
    (3.7567, 0.5),
    (4.4373, 1.0),
    (4.8853, 2.0),
    (5.0000, 5.0),
]

# ── Fixed parameters ──────────────────────────────────────────────────────────
T_amb        = 13.2 + 273.15
P_low        = 1.01325e5
P_hot        = 1e5
m_dot_air    = 1.0
T_salt_limit = 565 + 273.15

eta_cp     = 0.75
eta_tb     = 0.85
eta_heater = 0.90
eta_recup  = 0.85
eta_cooler = 0.95

T_hot_values = np.arange(600, 1001, 25)  # °C

# ── Colormap: qualitative palette + distinct markers ─────────────────────────
COLORS  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
MARKERS = ["v", "D", "^", "s", "o"]

# ── Sweep ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

for (PR, W_MW), color, mk in zip(PR_configs, COLORS, MARKERS):
    P_high = P_low * PR
    P_mid  = np.sqrt(P_low * P_high)
    eta_list = []

    for T_hot_C in T_hot_values:
        T_hot_su = T_hot_C + 273.15
        hot_fluid = 'SolarSalt' if T_hot_su <= T_salt_limit else 'Air'

        def make_hot():
            if hot_fluid == 'SolarSalt':
                h = SolarSaltConnector(); h.set_properties(T=T_hot_su, p=P_hot, m_dot=1.0)
            else:
                h = MassConnector(); h.set_properties(fluid='Air', T=T_hot_su, p=P_hot, m_dot=1.0)
            return h

        AirInlet = MassConnector()
        AirInlet.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=m_dot_air)
        CSource = MassConnector()
        CSource.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=10.0)

        try:
            cycle, c1, c2, t1, t2, heater, reheater, recup, ic = brayton_recuperated_rh_ic(
                eta_cp, eta_tb, eta_heater, eta_recup, eta_cooler,
                make_hot(), make_hot(), AirInlet, CSource,
                P_low, P_high, P_mid, T_tb_su_guess=T_hot_su)
            perf = compute_cycle_performance(
                c1, t1, heater, hot_fluid, T_hot_su, P_hot,
                recuperator=recup, comp2=c2, turb2=t2,
                reheater=reheater, intercooler=ic)
            eta_list.append(perf['eta'] * 100)
        except Exception as e:
            print(f"PR={PR} | T={T_hot_C}°C | FAILED: {e}")
            eta_list.append(np.nan)

    ax.plot(T_hot_values, eta_list, color=color, lw=2.5, marker=mk, ms=8,
            label=f'{W_MW} MW  (PR = {PR:.2f})', zorder=3)

ax.set_xlabel(r'Hot source temperature $T_{\mathrm{hot,in}}$ [°C]')
ax.set_ylabel(r'Cycle efficiency $\eta$ [%]')
ax.legend(loc='best')
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:.0f}'))

fig.tight_layout()
path = os.path.join(SAVE_DIR, 'fig_air_brayton_PR.pdf')
fig.savefig(path, bbox_inches='tight')
print(f'Saved → {path}')
plt.show()