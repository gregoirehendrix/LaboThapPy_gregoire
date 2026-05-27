# -*- coding: utf-8 -*-
"""
Air Brayton efficiency vs T_hot — comparison of Simple / Recuperated / Recuperated+IC+RH
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

from labothappy.connector.mass_connector       import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector

from air_Brayton import (brayton_simple, brayton_recuperated,
                                brayton_recuperated_rh_ic, compute_cycle_performance)

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 28, 'axes.labelsize': 28,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 20,
    'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25, 'grid.linestyle': '--',
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\3. Thermodynamic analysis\2. Power Blocks")
os.makedirs(SAVE_DIR, exist_ok=True)

# ── Fixed parameters ──────────────────────────────────────────────────────────
T_amb     = 13.2 + 273.15
P_low     = 1.01325e5
PR        = 5
P_high    = P_low * PR
P_mid     = np.sqrt(P_low * P_high)
P_hot     = 1e5
m_dot_air = 1.0

eta_cp     = 0.75
eta_tb     = 0.85
eta_heater = 0.90
eta_recup  = 0.85
eta_cooler = 0.95

T_salt_limit = 565 + 273.15

T_hot_values = np.arange(600, 1001, 25)   # °C

# ── Sweep ─────────────────────────────────────────────────────────────────────
eta_simple, eta_recup_list, eta_full = [], [], []

for T_hot_C in T_hot_values:
    T_hot_su = T_hot_C + 273.15

    if T_hot_su <= T_salt_limit:
        def make_hot():
            h = SolarSaltConnector(); h.set_properties(T=T_hot_su, p=P_hot, m_dot=1.0); return h
        hot_fluid = 'SolarSalt'
    else:
        def make_hot():
            h = MassConnector(); h.set_properties(fluid='Air', T=T_hot_su, p=P_hot, m_dot=1.0); return h
        hot_fluid = 'Air'

    AirInlet = MassConnector()
    AirInlet.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=m_dot_air)
    CSource = MassConnector()
    CSource.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=10.0)

    # Simple
    try:
        cycle, comp, turb, heater = brayton_simple(
            eta_cp, eta_tb, eta_heater, make_hot(), AirInlet, P_low, P_high)
        perf = compute_cycle_performance(comp, turb, heater, hot_fluid, T_hot_su, P_hot)
        eta_simple.append(perf['eta'] * 100)
    except Exception as e:
        print(f"Simple failed at {T_hot_C}°C: {e}")
        eta_simple.append(np.nan)

    # Recuperated
    try:
        AirInlet2 = MassConnector()
        AirInlet2.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=m_dot_air)
        cycle, comp, turb, heater, recup = brayton_recuperated(
            eta_cp, eta_tb, eta_heater, eta_recup,
            make_hot(), AirInlet2, P_low, P_high, T_tb_su_guess=T_hot_su)
        perf = compute_cycle_performance(comp, turb, heater, hot_fluid, T_hot_su, P_hot,
                                         recuperator=recup)
        eta_recup_list.append(perf['eta'] * 100)
    except Exception as e:
        print(f"Recuperated failed at {T_hot_C}°C: {e}")
        eta_recup_list.append(np.nan)

    # Recuperated + IC + RH
    try:
        AirInlet3 = MassConnector()
        AirInlet3.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=m_dot_air)
        cycle, c1, c2, t1, t2, heater, reheater, recup2, ic = brayton_recuperated_rh_ic(
            eta_cp, eta_tb, eta_heater, eta_recup, eta_cooler,
            make_hot(), make_hot(), AirInlet3, CSource,
            P_low, P_high, P_mid, T_tb_su_guess=T_hot_su)
        perf = compute_cycle_performance(c1, t1, heater, hot_fluid, T_hot_su, P_hot,
                                         recuperator=recup2, comp2=c2, turb2=t2,
                                         reheater=reheater, intercooler=ic)
        eta_full.append(perf['eta'] * 100)
    except Exception as e:
        print(f"Recuperated+IC+RH failed at {T_hot_C}°C: {e}")
        eta_full.append(np.nan)

print(f"\n{'TIT':>6} | {'Simple':>10} | {'Recuperated':>12} | {'Rec+IC+RH':>10}")
print("-" * 50)
for i, T in enumerate(T_hot_values):
    s  = f"{eta_simple[i]:.2f}%" if not np.isnan(eta_simple[i]) else "  N/A"
    r  = f"{eta_recup_list[i]:.2f}%" if not np.isnan(eta_recup_list[i]) else "  N/A"
    f  = f"{eta_full[i]:.2f}%" if not np.isnan(eta_full[i]) else "  N/A"
    print(f"{T:>6}°C | {s:>10} | {r:>12} | {f:>10}")
    
# ── Plot ──────────────────────────────────────────────────────────────────────
C_SIMPLE = '#d62728'
C_RECUP  = '#ff7f0e'
C_FULL   = '#1f77b4'

fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(T_hot_values, eta_simple,     color=C_SIMPLE, linewidth=2.5, label='Simple')
ax.plot(T_hot_values, eta_recup_list, color=C_RECUP,  linewidth=2.5, label='Recuperated')
ax.plot(T_hot_values, eta_full,       color=C_FULL,   linewidth=2.5, label='Recuperated + IC + RH')

ax.set_xlabel(r'Hot source temperature $T_{\mathrm{hot,in}}$ [°C]')
ax.set_ylabel(r'Cycle efficiency $\eta$ [%]')
ax.legend(loc='upper left')
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:.0f}%'))

fig.tight_layout()
path = os.path.join(SAVE_DIR, 'fig_air_brayton_TIT.pdf')
fig.savefig(path, bbox_inches='tight')
print(f'Saved → {path}')
plt.show()