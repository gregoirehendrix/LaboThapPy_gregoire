

# -*- coding: utf-8 -*-
"""
extract_TQ_data.py
==================
Extrait les températures et chaleurs nécessaires pour le diagramme T-Q
du cycle combiné Air Brayton Simple + Steam Rankine Siemens.

Utilisation :
    Lancer ce script depuis le répertoire où se trouvent
    combined_progressive.py, csp_piping.py, steam_rankine.py.
    Il affiche un bloc de données Python prêt à copier-coller
    dans plot_TQ_combined.py.

    python extract_TQ_data.py
"""

import sys
import numpy as np
from CoolProp.CoolProp import PropsSI

# ── Import du code principal ──────────────────────────────────────────────────
from combined_cycles_T_PR_sweep1 import (
    brayton_simple, compute_brayton_perf, _solve_siemens,
    MassConnector,
)

# ════════════════════════════════════════════════════════════════════════════
#   POINTS DE FONCTIONNEMENT À EXTRAIRE
#   Ajouter ou supprimer des entrées selon les cas souhaités dans le plot
# ════════════════════════════════════════════════════════════════════════════

OPERATING_POINTS = [
    dict(T_hot_C=700,  PR=3.7567, label="700 °C"),
    dict(T_hot_C=800,  PR=3.7567, label="800 °C"),
    dict(T_hot_C=1000, PR=3.7567, label="1000 °C"),
]

# ════════════════════════════════════════════════════════════════════════════
#   PARAMÈTRES FIXES (mêmes que dans combined_progressive.py)
# ════════════════════════════════════════════════════════════════════════════

T_AMB   = 13.2 + 273.15   # K
P_LOW_B = 1.01325e5        # Pa  — pression basse Brayton
ETA_CP  = 0.75
ETA_TB  = 0.85
ETA_HX  = 0.90

ETA_PUMP_ST = 0.80
ETA_TURB_ST = 0.93
ETA_HX_ST   = 0.95
PINCH_HX    = 5.0          # K
PINCH_COND  = 5.0          # K
P_HIGH_ST   = 160e5        # Pa
P_RH_ST     = 38.41e5      # Pa
P_LOW_ST    = 0.095e5      # Pa
M_DOT_AIR   = 20.0         # kg/s  (interne au solveur Siemens)
M_DOT_ST    = 1.0          # kg/s  (base de normalisation)


# ════════════════════════════════════════════════════════════════════════════
#   FONCTION D'EXTRACTION
# ════════════════════════════════════════════════════════════════════════════

def extract_one(T_hot_C, PR, label):
    """
    Résout le top cycle + bottom cycle à un point de fonctionnement donné
    et retourne un dictionnaire de toutes les températures et chaleurs.
    """
    P_high_b = P_LOW_B * PR
    T_hot_K  = T_hot_C + 273.15

    # ── Top cycle ────────────────────────────────────────────────────────────
    AirInlet = MassConnector()
    AirInlet.set_properties(fluid='Air', T=T_AMB, P=P_LOW_B, m_dot=1.0)
    HSource = MassConnector()
    HSource.set_properties(fluid='Air', T=T_hot_K, p=1e5, m_dot=1.0)

    _, comp, turb, heater = brayton_simple(
        ETA_CP, ETA_TB, ETA_HX, HSource, AirInlet, P_LOW_B, P_high_b
    )
    br = compute_brayton_perf(comp, turb, heater, 'Air', T_hot_K, 1e5)
    T_exhaust = br['T_exhaust']

    # ── Bottom cycle ─────────────────────────────────────────────────────────
    pump, eco, eva, sh, rh, turb_hp, turb_ic, cond = _solve_siemens(
        eta_pump   = ETA_PUMP_ST,
        eta_turb   = ETA_TURB_ST,
        eta_hx     = ETA_HX_ST,
        pinch_hx   = PINCH_HX,
        pinch_cond = PINCH_COND,
        T_hot_K    = T_exhaust,
        T_amb_K    = T_AMB,
        m_dot_air  = M_DOT_AIR,
        P_hot      = P_LOW_B,
        P_low      = P_LOW_ST,
        P_high     = P_HIGH_ST,
        P_rh       = P_RH_ST,
        m_dot_st   = M_DOT_ST,
    )

    # ── Chaleurs [kW / (kg_steam/s)] ─────────────────────────────────────────
    q_eco = M_DOT_ST * (eco.ex_C.h - eco.su_C.h) / 1000
    q_eva = M_DOT_ST * (eva.ex_C.h - eva.su_C.h) / 1000
    q_sh  = M_DOT_ST * (sh.ex_C.h  - sh.su_C.h)  / 1000
    q_rh  = M_DOT_ST * (rh.ex_C.h  - rh.su_C.h)  / 1000

    return {
        'label':    label,
        'T_hot_C':  T_hot_C,
        'PR':       PR,
        # Brayton
        'T_exhaust': T_exhaust - 273.15,
        'eta_br':    br['eta'] * 100,
        # Chaleurs
        'q_eco': q_eco,
        'q_eva': q_eva,
        'q_sh':  q_sh,
        'q_rh':  q_rh,
        # Températures vapeur [°C] — chemin principal ECO → EVA → SH
        'T_pump_ex':  pump.ex.T  - 273.15,
        'T_eco_ex_C': eco.ex_C.T - 273.15,
        'T_sat':      eva.su_C.T - 273.15,
        'T_sh_ex_C':  sh.ex_C.T  - 273.15,   # = TIT
        # Températures air [°C] — contre-courant
        'T_stack':    eco.ex_H.T - 273.15,
        'T_eco_su_H': eco.su_H.T - 273.15,
        'T_eva_ex_H': eva.ex_H.T - 273.15,
        'T_eva_su_H': eva.su_H.T - 273.15,
        'T_sh_ex_H':  sh.ex_H.T  - 273.15,
        'T_sh_su_H':  sh.su_H.T  - 273.15,   # = T_exhaust côté SG
        # Températures RH [°C]
        'T_rh_su_C':  rh.su_C.T - 273.15,    # sortie HP turbine
        'T_rh_ex_C':  rh.ex_C.T - 273.15,    # entrée LP turbine
        'T_rh_ex_H':  rh.ex_H.T - 273.15,    # air sortant côté RH
        'T_rh_su_H':  rh.su_H.T - 273.15,    # air entrant côté RH
    }


# ════════════════════════════════════════════════════════════════════════════
#   AFFICHAGE — FORMAT PRÊT À COPIER DANS plot_TQ_combined.py
# ════════════════════════════════════════════════════════════════════════════

def print_data_block(data_list):
    print("\n" + "=" * 70)
    print("  COPIER-COLLER dans plot_TQ_combined.py  →  section DATA")
    print("=" * 70)
    print("\nDATA = [")
    for d in data_list:
        print(f"    dict(")
        print(f"        label        = '{d['label']}',")
        print(f"        T_hot_C      = {d['T_hot_C']},")
        print(f"        PR           = {d['PR']},")
        print(f"        T_exhaust    = {d['T_exhaust']:.2f},   # °C — sortie turbine Brayton")
        print(f"        eta_br       = {d['eta_br']:.2f},     # %")
        print(f"        # Chaleurs [kW / kg_steam/s]")
        print(f"        q_eco        = {d['q_eco']:.2f},")
        print(f"        q_eva        = {d['q_eva']:.2f},")
        print(f"        q_sh         = {d['q_sh']:.2f},")
        print(f"        q_rh         = {d['q_rh']:.2f},")
        print(f"        # Températures vapeur [°C]")
        print(f"        T_pump_ex    = {d['T_pump_ex']:.2f},   # sortie pompe HP")
        print(f"        T_eco_ex_C   = {d['T_eco_ex_C']:.2f},   # sortie éco ≈ liq sat")
        print(f"        T_sat        = {d['T_sat']:.2f},   # T_sat(P_high)")
        print(f"        T_TIT        = {d['T_sh_ex_C']:.2f},   # TIT = sortie SH")
        print(f"        # Températures air [°C]")
        print(f"        T_stack      = {d['T_stack']:.2f},   # sortie éco côté air")
        print(f"        T_eco_su_H   = {d['T_eco_su_H']:.2f},   # entrée éco = sortie EVA côté air")
        print(f"        T_eva_ex_H   = {d['T_eva_ex_H']:.2f},   # sortie EVA côté air")
        print(f"        T_eva_su_H   = {d['T_eva_su_H']:.2f},   # entrée EVA = sortie SH côté air")
        print(f"        T_sh_ex_H    = {d['T_sh_ex_H']:.2f},   # sortie SH côté air")
        print(f"        T_sh_su_H    = {d['T_sh_su_H']:.2f},   # entrée SH = T_exhaust")
        print(f"        # RH [°C]")
        print(f"        T_rh_su_C    = {d['T_rh_su_C']:.2f},   # sortie HP turbine")
        print(f"        T_rh_ex_C    = {d['T_rh_ex_C']:.2f},   # entrée LP turbine (après RH)")
        print(f"        T_rh_ex_H    = {d['T_rh_ex_H']:.2f},   # air sortant RH")
        print(f"        T_rh_su_H    = {d['T_rh_su_H']:.2f},   # air entrant RH = T_exhaust")
        print(f"    ),")
    print("]\n")


# ════════════════════════════════════════════════════════════════════════════
#   MAIN
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    data_list = []
    for op in OPERATING_POINTS:
        print(f"→ Extraction : T={op['T_hot_C']}°C, PR={op['PR']} ...")
        try:
            d = extract_one(**op)
            data_list.append(d)
            print(f"   T_exhaust = {d['T_exhaust']:.1f}°C  |  "
                  f"T_TIT = {d['T_TIT']:.1f}°C  |  "
                  f"T_stack = {d['T_stack']:.1f}°C  |  "
                  f"q_eco+eva+sh = {d['q_eco']+d['q_eva']+d['q_sh']:.0f} kW")
        except Exception as e:
            print(f"   ⚠ FAILED : {e}")

    print_data_block(data_list)