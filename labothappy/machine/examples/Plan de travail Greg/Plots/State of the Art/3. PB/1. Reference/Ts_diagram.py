# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:08:09 2026

@author: gregoire.hendrix
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from CoolProp.CoolProp import PropsSI
import os

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'        : 'serif',
    'font.size'          : 28,
    'axes.labelsize'     : 28,
    'xtick.labelsize'    : 20,
    'ytick.labelsize'    : 24,
    'legend.fontsize'    : 18,
    'figure.dpi'         : 150,
    'axes.grid'          : True,
    'grid.alpha'         : 0.25,
    'grid.linestyle'     : '--',
    'text.usetex'        : True,
    'text.latex.preamble': r'\usepackage{amsmath}',
})

SAVE_DIR   = r"C:\Users\gregoire.hendrix@johncockerill.com\OneDrive - John Cockerill\Documents\Cockerill\Images\1. State of the Art\3. Power blocks\Reference solution"
font_graph = 23

fluid = "Water"

# ============================================================
#  PARAMETRES
# ============================================================
P_haute     = 40e5   # Pression haute (Pa)
P_basse     = 2e5    # Pression basse (Pa)
T_surech    = 400    # Température de surchauffe °C
eta_turbine = 0.80   # Efficacité isentropique turbine

T_surech_K = T_surech + 273.15

# --- Courbe de saturation ---
T_sat_range = np.linspace(273.16, PropsSI("Tcrit", fluid) - 0.05, 500)
s_liq, s_vap = [], []
for T in T_sat_range:
    try:
        s_liq.append(PropsSI("S","T",T,"Q",0,fluid)/1000)
        s_vap.append(PropsSI("S","T",T,"Q",1,fluid)/1000)
    except:
        s_liq.append(np.nan); s_vap.append(np.nan)
s_liq = np.array(s_liq); s_vap = np.array(s_vap)
s_crit = PropsSI("S","P",PropsSI("Pcrit",fluid),"Q",0.5,fluid)/1000
T_crit = PropsSI("Tcrit", fluid) - 273.15

# ============================================================
#  POINTS COMMUNS
# ============================================================
T1 = PropsSI("T","P",P_basse,"Q",0,fluid) - 273.15
s1 = PropsSI("S","P",P_basse,"Q",0,fluid) / 1000
h1 = PropsSI("H","P",P_basse,"Q",0,fluid)

s2 = s1
h2 = PropsSI("H","P",P_haute,"S",s2*1000,fluid)
T2 = PropsSI("T","P",P_haute,"S",s2*1000,fluid) - 273.15

T3 = PropsSI("T","P",P_haute,"Q",1,fluid) - 273.15
s3 = PropsSI("S","P",P_haute,"Q",1,fluid) / 1000
h3 = PropsSI("H","P",P_haute,"Q",1,fluid)

# ============================================================
#  CYCLE SANS SURCHAUFFE — turbine 3→4
# ============================================================
h4s = PropsSI("H","P",P_basse,"S",s3*1000,fluid)
h4  = h3 - eta_turbine * (h3 - h4s)
T4  = PropsSI("T","P",P_basse,"H",h4,fluid) - 273.15
s4  = PropsSI("S","P",P_basse,"H",h4,fluid) / 1000

# ============================================================
#  CYCLE AVEC SURCHAUFFE — turbine 3'→4'
# ============================================================
s3p = PropsSI("S","T",T_surech_K,"P",P_haute,fluid) / 1000
h3p = PropsSI("H","T",T_surech_K,"P",P_haute,fluid)
T3p = T_surech

h4ps = PropsSI("H","P",P_basse,"S",s3p*1000,fluid)
h4p  = h3p - eta_turbine * (h3p - h4ps)
T4p  = PropsSI("T","P",P_basse,"H",h4p,fluid) - 273.15
s4p  = PropsSI("S","P",P_basse,"H",h4p,fluid) / 1000

# ============================================================
#  CHEMINS INTERMEDIAIRES
# ============================================================
T_liqsat = PropsSI("T","P",P_haute,"Q",0,fluid)

T_liq_range = np.linspace(T2+273.15, T_liqsat - 0.01, 100)
s_liq_range = np.array([PropsSI("S","T",T,"P",P_haute,fluid)/1000
                         for T in T_liq_range])

s_diph = np.linspace(PropsSI("S","P",P_haute,"Q",0,fluid)/1000, s3, 80)
T_diph = np.full(80, T3)

T_surech_range = np.linspace(T3+273.15+0.01, T_surech_K, 100)
s_surech_range = np.array([PropsSI("S","T",T,"P",P_haute,fluid)/1000
                            for T in T_surech_range])

T_cond = PropsSI("T","P",P_basse,"Q",0,fluid) - 273.15

s_cond_1 = np.linspace(s4,  s1, 100)
T_cond_1  = np.full(100, T_cond)

s_cond_2 = np.linspace(s4p, s1, 100)
T_cond_2  = np.full(100, T_cond)

# ============================================================
#  FIGURE
# ============================================================
fig, ax = plt.subplots(figsize=(14, 8))

# Courbe de saturation
ax.fill_betweenx(T_sat_range-273.15, s_liq, s_vap, alpha=0.07, color="#60a5fa")
ax.plot(s_liq, T_sat_range-273.15, color="#3b82f6", lw=1.8)
ax.plot(s_vap, T_sat_range-273.15, color="#3b82f6", lw=1.8)
ax.plot(s_crit, T_crit, "ko", ms=6, zorder=5)

# --- CYCLE SANS SURCHAUFFE (gris) ---
cg = "#555555"
ax.plot(s_liq_range, T_liq_range-273.15, color=cg, lw=2)
ax.plot(s_diph,      T_diph,             color=cg, lw=2)
ax.plot([s3,  s4 ], [T3,  T4 ],         color=cg, lw=2)
ax.plot(s_cond_1,    T_cond_1,           color=cg, lw=2)
ax.plot([s1,  s2 ], [T1,  T2 ],         color=cg, lw=2, ls="--")

h4s_T = PropsSI("T","P",P_basse,"H",h4s,fluid) - 273.15
h4s_s = PropsSI("S","P",P_basse,"H",h4s,fluid) / 1000
ax.plot([s3, h4s_s], [T3, h4s_T], color=cg, lw=1, ls=":")

for s, T, lbl, dx, dy in [
    (s1,T1,"1",-12,-12),(s2,T2,"2",-16,4),
    (s3,T3,"3",5,4),(s4,T4,"4",5,-14),
]:
    ax.plot(s, T, "o", color=cg, ms=7, zorder=6)
    ax.annotate(lbl,(s,T),textcoords="offset points",
                xytext=(dx,dy),fontsize=16,color=cg,fontweight="bold")

# --- CYCLE AVEC SURCHAUFFE (rouge) ---
cr = "#dc2626"
ax.plot(s_liq_range,    T_liq_range-273.15,    color=cr, lw=2)
ax.plot(s_diph,         T_diph,                color=cr, lw=2)
ax.plot(s_surech_range, T_surech_range-273.15, color=cr, lw=2)
ax.plot([s3p, s4p],    [T3p, T4p],             color=cr, lw=2)
ax.plot(s_cond_2,       T_cond_2,              color=cr, lw=2)
ax.plot([s1, s2],      [T1, T2],               color=cr, lw=2, ls="--")

h4ps_T = PropsSI("T","P",P_basse,"H",h4ps,fluid) - 273.15
h4ps_s = PropsSI("S","P",P_basse,"H",h4ps,fluid) / 1000
ax.plot([s3p, h4ps_s], [T3p, h4ps_T], color=cr, lw=1, ls=":")

for s, T, lbl, dx, dy in [
    (s3p,T3p,"3'",5,4),(s4p,T4p,"4'",5,-14),
]:
    ax.plot(s, T, "o", color=cr, ms=7, zorder=6)
    ax.annotate(lbl,(s,T),textcoords="offset points",
                xytext=(dx,dy),fontsize=16,color=cr,fontweight="bold")

# Ellipse zone surchauffe
ell = Ellipse(xy=(6.8,250), width=1.8, height=360,
              fill=False, edgecolor="#f97316", lw=2.5, zorder=4)
ax.add_patch(ell)

# Annotations zones
ax.text(1,   200, "Subcooled\nliquid",  fontsize=font_graph, color="#64748b",
        ha="center", style="italic")
ax.text(4.0,  80, "Two-phase region",   fontsize=font_graph, color="#1d4ed8",
        ha="center", style="italic")
ax.text(9,   280, "Superheated\nvapor", fontsize=font_graph, color="#64748b",
        ha="center", style="italic")

ax.set_xlabel("Entropy")
ax.set_ylabel("Temperature")
ax.set_xticks([])
ax.set_yticks([])

for side in ["top","right","left","bottom"]:
    ax.spines[side].set_visible(False)

arrow_style = dict(arrowstyle="->", lw=2, color="black", mutation_scale=30)
ax.annotate("", xy=(1,0), xytext=(0,0),
            xycoords=("axes fraction","data"),
            textcoords=("axes fraction","data"),
            arrowprops=arrow_style)
ax.annotate("", xy=(-0.3,1), xytext=(-0.3,0),
            xycoords=("data","axes fraction"),
            textcoords=("data","axes fraction"),
            arrowprops=arrow_style)

ax.set_xlim(-0.3, 10.5)
ax.set_ylim(0, T_surech + 120)
ax.grid(True, linestyle=":", alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, "rankine_Ts.pdf"), dpi=150, bbox_inches="tight")
plt.show()