# -*- coding: utf-8 -*-
"""
T-s diagram — sCO2 Recompression (Hanwha type)
T_hot = 565°C   —   Author: Grégoire Hendrix

State numbering:
  1  Turbine inlet         (HP, ~563°C)
  2  Turbine outlet        (LP, ~457°C)
  3  RecupHT hot exit      (LP, ~228°C)
  4  Split point           (LP, ~105°C)
  5  MC suction            (LP,  37.1°C)
  6  MC discharge          (HP,  ~87°C)
  7  RC exit               (HP, ~240°C)
  8  Mixing point MC+RC    (HP, ~113°C)   ← formerly unlabelled
  9  LTR cold exit         (HP, ~220°C)   ← formerly state 8
  10 HTR cold exit         (HP, ~420°C)   ← formerly state 9
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from CoolProp.CoolProp import PropsSI

# ── Model import ──────────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sCO2_clean import sco2_cycle

# ── Plot style ────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 28, 'axes.labelsize': 28,
    'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 18,
    'figure.dpi': 150, 'axes.grid': True, 'grid.alpha': 0.25,
    'grid.linestyle': '--',
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\3. Thermodynamic analysis\2. Power Blocks")
os.makedirs(SAVE_DIR, exist_ok=True)

F = 'CO2'
def _s(T_C, P_bar):
    return PropsSI('S', 'T', T_C + 273.15, 'P', P_bar * 1e5, F) / 1e3

def _ib(Ta, Tb, P, n=300):
    T = np.linspace(Ta, Tb, n)
    return np.array([PropsSI('S', 'T', t + 273.15, 'P', P * 1e5, F) / 1e3
                     for t in T]), T

# ── State points from model ───────────────────────────────────────────────────
r   = sco2_cycle(T_hot=565, m_dot=1.0)
P_LP = 89.60
P_HP = 215.60

ST = {
    1:  (r['T1'],    P_HP),   # turbine inlet
    2:  (r['T2'],    P_LP),   # turbine outlet
    3:  (r['T4'],    P_LP),   # RecupHT hot exit
    4:  (r['T5'],    P_LP),   # split point
    5:  (37.12,      P_LP),   # MC suction
    6:  (r['T9'],    P_HP),   # MC discharge
    7:  (r['T_rc'],  P_HP),   # RC exit
    8:  (r['T_mix'], P_HP),   # mixing point MC+RC  ← now a labelled state
    9:  (r['T11'],   P_HP),   # LTR cold exit
    10: (r['T15'],   P_HP),   # HTR cold exit / heater inlet
}

S = {k: _s(*v) for k, v in ST.items()}

# ── Console summary ───────────────────────────────────────────────────────────
print(f"\n{'St':<4}{'T[°C]':>9}{'P[bar]':>8}{'s[kJ/kgK]':>12}")
for k in ST:
    print(f"  {k:<3} {ST[k][0]:>9.2f} {ST[k][1]:>8.2f} {S[k]:>12.4f}")

# ── Saturation dome ───────────────────────────────────────────────────────────
T_sat  = np.linspace(220.0, 303.9, 400)
s_liq  = np.array([PropsSI('S', 'T', T, 'Q', 0, F) for T in T_sat]) / 1e3
s_vap  = np.array([PropsSI('S', 'T', T, 'Q', 1, F) for T in T_sat]) / 1e3
s_crit = PropsSI('S', 'T', 304.0, 'Q', 0.5, F) / 1e3

# ── Colours ───────────────────────────────────────────────────────────────────
C_HEAT  = '#d62728'
C_TURB  = '#1f77b4'
C_MC    = '#17becf'
C_COOL  = '#2ca02c'
C_RECUP = '#9467bd'
C_RC    = '#ff7f0e'
C_DOME  = 'k'

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

# Saturation dome
ax.plot(s_liq, T_sat - 273.15, color=C_DOME, lw=2.0, zorder=5)
ax.plot(s_vap, T_sat - 273.15, color=C_DOME, lw=2.0, zorder=5)
ax.scatter(s_crit, 30.85, color=C_DOME, s=60, zorder=6)
ax.annotate('Critical point (30.98 °C)', xy=(s_crit, 30.85),
            xytext=(s_crit - 0.03, -40), fontsize=16,
            arrowprops=dict(arrowstyle='->', lw=1.2))

# ── Cycle processes ───────────────────────────────────────────────────────────

# Heater 10→1
s, T = _ib(ST[10][0], ST[1][0], ST[10][1])
ax.plot(s, T, color=C_HEAT, lw=2.5, zorder=4)

# Turbine 1→2
ax.plot([S[1], S[2]], [ST[1][0], ST[2][0]], color=C_TURB, lw=2.5, zorder=4)

# RecupHT hot 2→3
s, T = _ib(ST[2][0], ST[3][0], ST[2][1])
ax.plot(s, T, color=C_RECUP, lw=2.5, zorder=4)

# RecupLT hot 3→4
s, T = _ib(ST[3][0], ST[4][0], ST[3][1])
ax.plot(s, T, color=C_RECUP, lw=2.5, zorder=4)

# Gas cooler 4→5
s, T = _ib(ST[4][0], ST[5][0], ST[4][1])
ax.plot(s, T, color=C_COOL, lw=2.5, zorder=4)

# Main compressor 5→6
ax.plot([S[5], S[6]], [ST[5][0], ST[6][0]], color=C_MC, lw=2.5, zorder=4)

# Recompressor 4→7
ax.plot([S[4], S[7]], [ST[4][0], ST[7][0]], color=C_RC, lw=2.5, zorder=4)

# Mixing: dotted lines 6→8 and 7→8 (topology indicator, not a T-s process)
ax.plot([S[6], S[8]], [ST[6][0], ST[8][0]], color=C_MC, lw=1.5, ls=':', zorder=3)
ax.plot([S[7], S[8]], [ST[7][0], ST[8][0]], color=C_RC,    lw=1.5, ls=':', zorder=3)

# LTR cold 8→9
s, T = _ib(ST[8][0], ST[9][0], P_HP)
ax.plot(s, T, color=C_RECUP, lw=2.5, ls='--', zorder=4)

# HTR cold 9→10
s, T = _ib(ST[9][0], ST[10][0], ST[9][1])
ax.plot(s, T, color=C_RECUP, lw=2.5, ls='--', zorder=4)

# ── State point markers + labels ─────────────────────────────────────────────
off = {
    1:  ( 0.01,  12),
    2:  ( 0.03,  12),
    3:  ( 0.03, -10),
    4:  ( 0.00, -35),
    5:  (-0.05,   4),
    6:  (-0.03,  12),
    7:  (-0.03,  12),
    8:  (-0.03,  12),   # below, avoids overlap with 6 and dotted lines
    9:  (-0.03,  12),
    10: (-0.06,  12),
}
for k, (T, P) in ST.items():
    ax.scatter(S[k], T, color='k', s=70, zorder=7)
    ds, dT = off[k]
    ax.text(S[k] + ds, T + dT, str(k), fontsize=21, fontweight='bold', zorder=8)

# ── Annotations ───────────────────────────────────────────────────────────────
ax.annotate('RC in\n(split from 4)', xy=(S[4], ST[4][0]),
            xytext=(2.3, 75), fontsize=15, color=C_RC, ha='center',
            arrowprops=dict(arrowstyle='->', color=C_RC, lw=1.2))

# Pressure labels
sall = list(S.values())
ax.text(max(sall) - 0.35, (ST[1][0] + ST[10][0]) / 2, '215.6 bar',
        fontsize=17, color=C_HEAT, style='italic')
ax.text(S[4] - 0.25, 36, '89.6 bar',
        fontsize=17, color=C_COOL, style='italic')

# ── Legend ────────────────────────────────────────────────────────────────────
leg = [
    Line2D([0],[0], color=C_HEAT,  lw=2.5,         label='Heat addition (heater)'),
    Line2D([0],[0], color=C_TURB,  lw=2.5,         label='Turbine'),
    Line2D([0],[0], color=C_MC,    lw=2.5,         label='Main compressor'),
    Line2D([0],[0], color=C_RC,    lw=2.5,         label='Recompressor'),
    Line2D([0],[0], color=C_RECUP, lw=2.5,         label='Recuperators — hot side'),
    Line2D([0],[0], color=C_RECUP, lw=2.5, ls='--',label='Recuperators — cold side'),
    Line2D([0],[0], color=C_COOL,  lw=2.5,         label='Gas cooler'),
    Line2D([0],[0], color=C_DOME,  lw=2.0,         label='Saturation dome (CO₂)'),
]
ax.legend(handles=leg, loc='upper left', fontsize=17)

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xlim(min(sall) - 0.15, max(sall) + 0.25)
ax.set_ylim(-80, 610)
ax.set_xlabel(r'Specific entropy  $s$  [kJ kg$^{-1}$ K$^{-1}$]')
ax.set_ylabel(r'Temperature  $T$  [°C]')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_locator(ticker.MultipleLocator(50))

fig.tight_layout()
path = os.path.join(SAVE_DIR, "fig_sco2_Ts.pdf")
fig.savefig(path, bbox_inches='tight')
print(f"Saved → {path}")
plt.show()