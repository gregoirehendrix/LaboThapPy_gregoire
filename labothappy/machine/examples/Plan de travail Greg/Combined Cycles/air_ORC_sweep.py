# -*- coding: utf-8 -*-
"""
combined_progressive.py
=======================
Étude progressive des cycles combinés Air Brayton + bottom cycle.

Logique centrale :
  brayton_simple           → T_exhaust élevée (~500 °C) → sCO2 (Hanwha type)
  brayton_recuperated      → T_exhaust moyenne (~300 °C) → Rankine vapeur simple
  brayton_recuperated_rh_ic→ T_exhaust basse  (~200 °C) → ORC R245fa

Sélection automatique du bottom cycle :
  T_exhaust > T_THRESH_STEAM (300 °C) → Rankine vapeur
  T_exhaust ≤ T_THRESH_STEAM          → ORC R245fa

Transport thermique :
  Intégré via heat_transport.py (NTU cylindrique rockwool, diamètre selon N).

Bottom cycles analytiques (scipy/CoolProp) — robustes et rapides,
indépendants du solveur LaboThapPy.

Author : Grégoire Hendrix
"""

#%%
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_it   import IterativeCircuit
from labothappy.machine.circuit_rec  import RecursiveCircuit
from labothappy.connector.mass_connector       import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.compressor.compressor_csteff          import CompressorCstEff
from labothappy.component.expander.expander_csteff              import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff             import HexCstEff
from labothappy.component.heat_exchanger.hex_csteff_disc        import HexCstEffDisc
from labothappy.component.heat_exchanger.hex_cstpinch           import HexCstPinch
from labothappy.component.pump.pump_csteff                      import PumpCstEff

from csp_piping import network_heat_loss, select_fluid, fluid_label


# ════════════════════════════════════════════════════════════════════════════
#   SEUILS DE SÉLECTION DU BOTTOM CYCLE
# ════════════════════════════════════════════════════════════════════════════

T_THRESH_STEAM = 300.0 + 273.15   # K  — T_exhaust > seuil → Steam, sinon ORC

BOTTOM_STEAM = "Steam"
BOTTOM_ORC   = "ORC_R245fa"


def select_bottom_cycle(T_exhaust_K):
    """Sélection automatique du bottom cycle selon T_exhaust [K]."""
    if T_exhaust_K > T_THRESH_STEAM:
        return BOTTOM_STEAM
    else:
        return BOTTOM_ORC


def bottom_label(bc, bp=None):
    if bc == BOTTOM_STEAM:
        return "Rankine Vapeur (Siemens)"
    else:
        fluid_lbl = bp.get('fluid_lbl', 'R245fa') if bp else 'ORC'
        return f"ORC {fluid_lbl}"


# ════════════════════════════════════════════════════════════════════════════
#   TOP CYCLES — BRAYTON AIR  (versions validées)
# ════════════════════════════════════════════════════════════════════════════

def brayton_simple(eta_cp, eta_tb, eta_heater,
                   HSource, AirInlet, P_low, P_high):
    """Brayton simple — IterativeCircuit, pas de récupérateur."""
    cycle      = IterativeCircuit(fluid='Air')
    compressor = CompressorCstEff()
    heater     = HexCstEff()
    turbine    = ExpanderCstEff()

    compressor.set_parameters(eta_is=eta_cp)
    turbine.set_parameters(eta_is=eta_tb)
    heater.set_parameters(eta=eta_heater)

    cycle.add_component(compressor, "Compressor")
    cycle.add_component(heater,     "Heater")
    cycle.add_component(turbine,    "Turbine")

    cycle.link_components("Compressor", "m-ex",   "Heater",  "m-su_C")
    cycle.link_components("Heater",     "m-ex_C", "Turbine", "m-su")

    cycle.add_source("AirInlet",  AirInlet, cycle.components["Compressor"], "m-su")
    cycle.add_source("HotSource", HSource,  cycle.components["Heater"],     "m-su_H")

    cycle.set_cycle_input(target="Compressor:ex", p=P_high)
    cycle.set_cycle_input(target="Turbine:ex",    p=P_low)

    cycle._build_solve_order()
    for name in cycle.solve_start_components:
        cycle.components[name].solve()

    return cycle, compressor, turbine, heater


def brayton_recuperated(eta_cp, eta_tb, eta_heater, eta_recup,
                        HSource, AirInlet, P_low, P_high, T_tb_su_guess):
    """Brayton récupéré — RecursiveCircuit."""
    cycle       = RecursiveCircuit('Air')
    compressor  = CompressorCstEff()
    recuperator = HexCstEff()
    heater      = HexCstEff()
    turbine     = ExpanderCstEff()

    compressor.set_parameters(eta_is=eta_cp)
    turbine.set_parameters(eta_is=eta_tb)
    heater.set_parameters(eta=eta_heater)
    recuperator.set_parameters(eta=eta_recup)

    cycle.add_component(compressor,  "Compressor")
    cycle.add_component(recuperator, "Recuperator")
    cycle.add_component(heater,      "Heater")
    cycle.add_component(turbine,     "Turbine")

    cycle.link_components("Compressor",  "m-ex",   "Recuperator", "m-su_C")
    cycle.link_components("Recuperator", "m-ex_C", "Heater",      "m-su_C")
    cycle.link_components("Heater",      "m-ex_C", "Turbine",     "m-su")
    cycle.link_components("Turbine",     "m-ex",   "Recuperator", "m-su_H")

    cycle.add_source("AirInlet",  AirInlet, cycle.components["Compressor"], "m-su")
    cycle.add_source("HotSource", HSource,  cycle.components["Heater"],     "m-su_H")

    cycle.set_fixed_properties(target="Compressor:ex", p=P_high)
    cycle.set_fixed_properties(target="Turbine:ex",    p=P_low)

    T_amb = AirInlet.T
    m_dot = AirInlet.m_dot

    h_cs    = PropsSI('H', 'T', T_amb,         'P', P_low,  'Air')
    s_cs    = PropsSI('S', 'T', T_amb,         'P', P_low,  'Air')
    h_cx_is = PropsSI('H', 'P', P_high,        'S', s_cs,   'Air')
    h_cx    = h_cs + (h_cx_is - h_cs) / eta_cp
    T_cx    = PropsSI('T', 'H', h_cx,          'P', P_high, 'Air')

    h_ts    = PropsSI('H', 'T', T_tb_su_guess, 'P', P_high, 'Air')
    s_ts    = PropsSI('S', 'T', T_tb_su_guess, 'P', P_high, 'Air')
    h_tx_is = PropsSI('H', 'P', P_low,         'S', s_ts,   'Air')
    h_tx    = h_ts - eta_tb * (h_ts - h_tx_is)
    T_tx    = PropsSI('T', 'H', h_tx,          'P', P_low,  'Air')

    T_rec_xC = T_cx + eta_recup * (T_tx - T_cx)
    T_rec_xH = T_tx - eta_recup * (T_tx - T_cx)

    cycle.set_cycle_guess(target='Recuperator:su_C', T=T_cx,     p=P_high, m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:su_H', T=T_tx,     p=P_low,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:ex_C', T=T_rec_xC, p=P_high)
    cycle.set_cycle_guess(target='Recuperator:ex_H', T=T_rec_xH, p=P_low)

    cycle.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance=1e-3)

    cycle.mute_print()
    cycle.solve(max_iter=50)

    return cycle, compressor, turbine, heater, recuperator


def brayton_recuperated_rh_ic(eta_cp, eta_tb, eta_heater, eta_recup, eta_cooler,
                               HSource_heater, HSource_reheater, AirInlet, CSource,
                               P_low, P_high, P_mid, T_tb_su_guess):
    """Brayton récupéré + réchauffage intermédiaire + refroidissement intermédiaire."""
    cycle       = RecursiveCircuit('Air')
    comp1       = CompressorCstEff()
    comp2       = CompressorCstEff()
    intercooler = HexCstEff()
    recuperator = HexCstEff()
    heater      = HexCstEff()
    turb1       = ExpanderCstEff()
    reheater    = HexCstEff()
    turb2       = ExpanderCstEff()

    comp1.set_parameters(eta_is=eta_cp)
    comp2.set_parameters(eta_is=eta_cp)
    turb1.set_parameters(eta_is=eta_tb)
    turb2.set_parameters(eta_is=eta_tb)
    heater.set_parameters(eta=eta_heater)
    reheater.set_parameters(eta=eta_heater)
    recuperator.set_parameters(eta=eta_recup)
    intercooler.set_parameters(eta=eta_cooler)

    cycle.add_component(comp1,       "Comp1")
    cycle.add_component(intercooler, "Intercooler")
    cycle.add_component(comp2,       "Comp2")
    cycle.add_component(recuperator, "Recuperator")
    cycle.add_component(heater,      "Heater")
    cycle.add_component(turb1,       "Turb1")
    cycle.add_component(reheater,    "Reheater")
    cycle.add_component(turb2,       "Turb2")

    cycle.link_components("Comp1",       "m-ex",   "Intercooler", "m-su_H")
    cycle.link_components("Intercooler", "m-ex_H", "Comp2",       "m-su")
    cycle.link_components("Comp2",       "m-ex",   "Recuperator", "m-su_C")
    cycle.link_components("Recuperator", "m-ex_C", "Heater",      "m-su_C")
    cycle.link_components("Heater",      "m-ex_C", "Turb1",       "m-su")
    cycle.link_components("Turb1",       "m-ex",   "Reheater",    "m-su_C")
    cycle.link_components("Reheater",    "m-ex_C", "Turb2",       "m-su")
    cycle.link_components("Turb2",       "m-ex",   "Recuperator", "m-su_H")

    cycle.add_source("AirInlet",   AirInlet,         cycle.components["Comp1"],       "m-su")
    cycle.add_source("HotSource",  HSource_heater,   cycle.components["Heater"],      "m-su_H")
    cycle.add_source("ReheatSrc",  HSource_reheater, cycle.components["Reheater"],    "m-su_H")
    cycle.add_source("ColdSource", CSource,          cycle.components["Intercooler"], "m-su_C")

    cycle.set_fixed_properties(target="Comp1:ex", p=P_mid)
    cycle.set_fixed_properties(target="Comp2:ex", p=P_high)
    cycle.set_fixed_properties(target="Turb1:ex", p=P_mid)
    cycle.set_fixed_properties(target="Turb2:ex", p=P_low)

    T_amb = AirInlet.T
    m_dot = AirInlet.m_dot

    h_c1s    = PropsSI('H', 'T', T_amb,         'P', P_low,  'Air')
    s_c1s    = PropsSI('S', 'T', T_amb,         'P', P_low,  'Air')
    h_c1x_is = PropsSI('H', 'P', P_mid,         'S', s_c1s,  'Air')
    h_c1x    = h_c1s + (h_c1x_is - h_c1s) / eta_cp
    T_c1x    = PropsSI('T', 'H', h_c1x,         'P', P_mid,  'Air')
    T_ic_ex  = T_amb + (1 - eta_cooler) * (T_c1x - T_amb)

    h_c2s    = PropsSI('H', 'T', T_ic_ex,       'P', P_mid,  'Air')
    s_c2s    = PropsSI('S', 'T', T_ic_ex,       'P', P_mid,  'Air')
    h_c2x_is = PropsSI('H', 'P', P_high,        'S', s_c2s,  'Air')
    h_c2x    = h_c2s + (h_c2x_is - h_c2s) / eta_cp
    T_c2x    = PropsSI('T', 'H', h_c2x,         'P', P_high, 'Air')

    h_t1s    = PropsSI('H', 'T', T_tb_su_guess, 'P', P_high, 'Air')
    s_t1s    = PropsSI('S', 'T', T_tb_su_guess, 'P', P_high, 'Air')
    h_t1x_is = PropsSI('H', 'P', P_mid,         'S', s_t1s,  'Air')
    h_t1x    = h_t1s - eta_tb * (h_t1s - h_t1x_is)
    T_t1x    = PropsSI('T', 'H', h_t1x,         'P', P_mid,  'Air')

    h_t2s    = PropsSI('H', 'T', T_tb_su_guess, 'P', P_mid,  'Air')
    s_t2s    = PropsSI('S', 'T', T_tb_su_guess, 'P', P_mid,  'Air')
    h_t2x_is = PropsSI('H', 'P', P_low,         'S', s_t2s,  'Air')
    h_t2x    = h_t2s - eta_tb * (h_t2s - h_t2x_is)
    T_t2x    = PropsSI('T', 'H', h_t2x,         'P', P_low,  'Air')

    T_rec_xC = T_c2x + eta_recup * (T_t2x - T_c2x)
    T_rec_xH = T_t2x - eta_recup * (T_t2x - T_c2x)

    cycle.set_cycle_guess(target='Comp1:su',         T=T_amb,         p=P_low,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Comp1:ex',         T=T_c1x,         p=P_mid)
    cycle.set_cycle_guess(target='Comp2:su',         T=T_ic_ex,       p=P_mid,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Comp2:ex',         T=T_c2x,         p=P_high)
    cycle.set_cycle_guess(target='Recuperator:su_C', T=T_c2x,         p=P_high, m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:ex_C', T=T_rec_xC,      p=P_high)
    cycle.set_cycle_guess(target='Recuperator:su_H', T=T_t2x,         p=P_low,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:ex_H', T=T_rec_xH,      p=P_low)
    cycle.set_cycle_guess(target='Turb1:su',         T=T_tb_su_guess, p=P_high, m_dot=m_dot)
    cycle.set_cycle_guess(target='Turb1:ex',         T=T_t1x,         p=P_mid)
    cycle.set_cycle_guess(target='Turb2:su',         T=T_tb_su_guess, p=P_mid,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Turb2:ex',         T=T_t2x,         p=P_low)

    cycle.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Turb1:ex',         variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Turb1:ex',         variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Comp2:su',         variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Comp2:su',         variable='p', tolerance=1e-3)

    cycle.mute_print()
    cycle.solve(max_iter=50)

    return cycle, comp1, comp2, turb1, turb2, heater, reheater, recuperator, intercooler


# ════════════════════════════════════════════════════════════════════════════
#   PERFORMANCE BRAYTON
# ════════════════════════════════════════════════════════════════════════════

def compute_brayton_perf(compressor, turbine, heater,
                         hot_fluid, T_hot_su, P_hot,
                         recuperator=None,
                         comp2=None, turb2=None,
                         reheater=None, intercooler=None):
    """Bilan Brayton — ref à 1 kg_air/s. Retourne dict avec clés courtes."""
    m = compressor.su.m_dot

    W_comp  = m * (compressor.ex.h - compressor.su.h) / 1000
    W_turb  = m * (turbine.su.h    - turbine.ex.h)    / 1000
    W_comp2 = m * (comp2.ex.h - comp2.su.h) / 1000 if comp2  is not None else 0.0
    W_turb2 = m * (turb2.su.h - turb2.ex.h) / 1000 if turb2  is not None else 0.0

    W_net   = (W_turb + W_turb2) - (W_comp + W_comp2)
    Q_in    = m * (heater.ex_C.h - heater.su_C.h) / 1000
    Q_rh    = m * (reheater.ex_C.h  - reheater.su_C.h)  / 1000 if reheater    is not None else 0.0
    Q_ic    = m * (intercooler.su_H.h - intercooler.ex_H.h) / 1000 if intercooler is not None else 0.0
    Q_total = Q_in + Q_rh

    # Température et enthalpie d'échappement (après récupérateur si présent)
    if recuperator is not None:
        T_exhaust = recuperator.ex_H.T
        h_exhaust = recuperator.ex_H.h
    elif turb2 is not None:
        T_exhaust = turb2.ex.T
        h_exhaust = turb2.ex.h
    else:
        T_exhaust = turbine.ex.T
        h_exhaust = turbine.ex.h

    return {
        'W_comp': W_comp, 'W_comp2': W_comp2,
        'W_turb': W_turb, 'W_turb2': W_turb2,
        'W_net':  W_net,
        'Q_in': Q_in, 'Q_rh': Q_rh, 'Q_ic': Q_ic, 'Q_total': Q_total,
        'Q_cooler': Q_total - W_net,
        'eta': W_net / Q_total,
        'T_exhaust': T_exhaust,
        'h_exhaust': h_exhaust,
    }


# ════════════════════════════════════════════════════════════════════════════
#   MODÈLES ORC PAR TAILLE — 1 / 2 / 5 MWel
#   Sources : Quoilin et al. (2013), Macchi & Astolfi (2017),
#             données Turboden / Ormat
# ════════════════════════════════════════════════════════════════════════════

ORC_MODELS = {
    # Clé : puissance nominale [MWel]
    # eta_exp   : rendement isentropique expandeur (turbine axiale)
    # eta_pump  : rendement isentropique pompe
    # eta_recup : efficacité récupérateur interne
    # dT_pinch  : pinch point échangeur évaporateur [K]
    # dT_sh     : surchauffe au-dessus de la saturation [K]
    1.0: dict(eta_exp=0.62, eta_pump=0.62, eta_recup=0.6, dT_pinch=15, dT_sh=5.0, label="1 MWel"),
    2.0: dict(eta_exp=0.67, eta_pump=0.65, eta_recup=0.65, dT_pinch=12, dT_sh=5.0, label="2 MWel"),
    5.0: dict(eta_exp=0.72, eta_pump=0.69, eta_recup=0.70, dT_pinch=10, dT_sh=5.0, label="5 MWel"),
}


def get_orc_params(P_nom_MW):
    """
    Retourne les paramètres ORC pour la taille nominale demandée.
    Valeurs disponibles : 1, 2, 5 MWel.
    Si P_nom_MW ne correspond pas exactement, retourne la clé la plus proche.

    Parameters
    ----------
    P_nom_MW : float — puissance nominale souhaitée [MWel]

    Returns
    -------
    dict avec eta_exp, eta_pump, eta_recup, dT_pinch, dT_sh, label
    """
    keys = sorted(ORC_MODELS.keys())   # [1.0, 2.0, 5.0]
    closest = min(keys, key=lambda k: abs(k - P_nom_MW))
    return dict(ORC_MODELS[closest])


# ════════════════════════════════════════════════════════════════════════════
#   BOTTOM CYCLES — analytiques (CoolProp + scipy)
# ════════════════════════════════════════════════════════════════════════════

# ── Helpers ─────────────────────────────────────────────────────────────────

def _h(fluid, T_K, P_Pa):
    return PropsSI('H', 'T', T_K, 'P', P_Pa, fluid)

def _s(fluid, T_K, P_Pa):
    return PropsSI('S', 'T', T_K, 'P', P_Pa, fluid)

def _hPS(fluid, P_Pa, s_J):
    return PropsSI('H', 'P', P_Pa, 'S', s_J, fluid)

def _expand(h_su, P_su, P_ex, eta, fluid):
    s_su    = PropsSI('S', 'H', h_su, 'P', P_su, fluid)
    h_ex_is = PropsSI('H', 'P', P_ex, 'S', s_su, fluid)
    return h_su - eta * (h_su - h_ex_is)

def _pump(h_su, P_su, P_ex, eta, fluid):
    s_su    = PropsSI('S', 'H', h_su, 'P', P_su, fluid)
    h_ex_is = PropsSI('H', 'P', P_ex, 'S', s_su, fluid)
    return h_su + (h_ex_is - h_su) / eta


# ── sCO2 Recompression (Hanwha-type, scipy) ──────────────────────────────────

# Constantes cycle Hanwha (back-calculées depuis PFD)
_SCO2 = dict(
    P_TB_SU=215.6e5, P_TB_EX=89.6e5,
    P_MC_SU=85.7e5,  P_MC_EX=224.6e5,
    P_RC_SU=86.7e5,  P_RC_EX=223.0e5,
    P_RHT_HI=220.0e5, P_RHT_LO=89.0e5,
    P_RLT_HI=222.0e5, P_RLT_LO=88.0e5,
    ETA_TB=0.9102, ETA_MC=0.8824, ETA_RC=0.7544,
    ETA_RHT=0.9385, ETA_RLT=0.8987,
    T_MC_SU=37.1, PINCH_HTR=1.2,
    SPLIT=47.68/68.11,
)


def bottom_sco2(T_hot_K, Q_available_kW, verbose=False):
    """
    Bottom cycle sCO2 recompression (Hanwha-type).

    T_hot_K         : température du fluide chaud arrivant au heater [K]
    Q_available_kW  : chaleur disponible côté air [kW / kg_air/s]

    Returns
    -------
    dict : W_net [kW], Q_heater [kW], eta [-], m_dot_co2 [kg/s / kg_air/s]
    """
    p = _SCO2
    alpha = p['SPLIT']
    beta  = 1.0 - alpha
    T1    = T_hot_K - p['PINCH_HTR']

    # Turbine
    h1  = _h('CO2', T1, p['P_TB_SU'])
    h2  = _expand(h1, p['P_TB_SU'], p['P_TB_EX'], p['ETA_TB'], 'CO2')

    # Main compressor
    h7  = _h('CO2', p['T_MC_SU'] + 273.15, p['P_MC_SU'])
    h9  = _pump(h7, p['P_MC_SU'], p['P_MC_EX'], p['ETA_MC'], 'CO2')

    def residuals(x):
        h4, h_mix, h_rc = x
        T4   = PropsSI('T', 'H', h4,    'P', p['P_RHT_LO'], 'CO2')
        T_mx = PropsSI('T', 'H', h_mix, 'P', p['P_RLT_HI'], 'CO2')

        Qmax_lt = min(h4 - _h('CO2', T_mx, p['P_RLT_LO']),
                      alpha * (_h('CO2', T4, p['P_RLT_HI']) - h_mix))
        Q_lt    = p['ETA_RLT'] * Qmax_lt
        h5      = h4    - Q_lt
        h11     = h_mix + Q_lt / alpha

        T11     = PropsSI('T', 'H', h11, 'P', p['P_RLT_HI'], 'CO2')
        Qmax_ht = min(h2 - _h('CO2', T11, p['P_RHT_LO']),
                      _h('CO2', PropsSI('T','H',h2,'P',p['P_RHT_HI'],'CO2'),
                         p['P_RHT_HI']) - h11)
        Q_ht    = p['ETA_RHT'] * Qmax_ht
        h4_c    = h2 - Q_ht

        h_rc_c  = _pump(h5, p['P_RC_SU'], p['P_RC_EX'], p['ETA_RC'], 'CO2')
        h_mix_c = alpha * h9 + beta * h_rc

        return [h4 - h4_c, h_mix - h_mix_c, h_rc - h_rc_c]

    x0  = [_h('CO2', 202.4+273.15, p['P_RHT_LO']),
            _h('CO2', 103.1+273.15, p['P_RLT_HI']),
            _h('CO2', 184.2+273.15, p['P_RC_EX'])]
    sol, _, ier, _ = fsolve(residuals, x0, full_output=True)
    if ier != 1:
        return None

    h4, h_mix, h_rc_ex = sol
    T4   = PropsSI('T', 'H', h4,    'P', p['P_RHT_LO'], 'CO2')
    T_mx = PropsSI('T', 'H', h_mix, 'P', p['P_RLT_HI'], 'CO2')

    Qmax_lt = min(h4 - _h('CO2', T_mx, p['P_RLT_LO']),
                  alpha * (_h('CO2', T4, p['P_RLT_HI']) - h_mix))
    Q_lt    = p['ETA_RLT'] * Qmax_lt
    h5      = h4 - Q_lt
    h11     = h_mix + Q_lt / alpha
    T11     = PropsSI('T', 'H', h11, 'P', p['P_RLT_HI'], 'CO2')

    Qmax_ht = min(h2 - _h('CO2', T11, p['P_RHT_LO']),
                  _h('CO2', PropsSI('T','H',h2,'P',p['P_RHT_HI'],'CO2'),
                     p['P_RHT_HI']) - h11)
    Q_ht    = p['ETA_RHT'] * Qmax_ht
    h15     = h11 + Q_ht

    Q_htr_sp = _h('CO2', T1, p['P_TB_SU']) - h15          # J/kg_CO2
    W_tb_sp  = h1 - h2
    W_mc_sp  = alpha * (h9 - h7)
    W_rc_sp  = beta  * (h_rc_ex - h5)
    W_net_sp = W_tb_sp - W_mc_sp - W_rc_sp
    eta_sp   = W_net_sp / Q_htr_sp

    m_dot_co2 = Q_available_kW * 1000 / Q_htr_sp
    W_net     = W_net_sp  * m_dot_co2 / 1000
    Q_heater  = Q_htr_sp  * m_dot_co2 / 1000

    return {
        'type':       'sCO2',
        'eta':        eta_sp,
        'W_net':      W_net,
        'Q_heater':   Q_heater,
        'm_dot':      m_dot_co2,
        'fluid':      'CO2',
        'T_hot_in':   T1,
    }


# ── Rankine Vapeur avec resurchauffe + dégazeur (regeneration) ──────────────

def bottom_steam(T_hot_K, Q_available_kW,
                 eta_pump=0.85, eta_turb=0.85,
                 P_high=160e5, P_rh=40e5, P_dea=5e5, P_low=0.10e5,
                 dT_pinch=10.0):
    """
    Bottom cycle Rankine vapeur avec resurchauffe + dégazeur.
    """
    fluid  = 'Water'
    T_crit = 647.096
    P_crit = 220.64e5

    T_su = T_hot_K - dT_pinch

    T_sat_max = min(T_su - 5.0, T_crit - 5.0)
    try:
        P_high_max = PropsSI('P', 'T', T_sat_max, 'Q', 1, fluid)
    except Exception:
        P_high_max = P_crit * 0.95
    P_high = min(P_high, P_high_max)
    P_high = max(P_high, P_rh * 1.5)

    P_rh  = min(P_rh,  P_high * 0.6)
    P_rh  = max(P_rh,  P_dea  * 2.0)
    P_dea = min(P_dea, P_rh   * 0.4)
    P_dea = max(P_dea, P_low  * 5.0)

    try:
        h_sat_liq_cond = PropsSI('H', 'P', P_low,  'Q', 0, fluid)
        h_sat_vap_cond = PropsSI('H', 'P', P_low,  'Q', 1, fluid)
        h_fw_dea       = PropsSI('H', 'P', P_dea,  'Q', 0, fluid)

        h_plp_ex = _pump(h_sat_liq_cond, P_low,  P_dea,  eta_pump, fluid)
        h_php_ex = _pump(h_fw_dea,       P_dea,  P_high, eta_pump, fluid)

        h_hp_su  = _h(fluid, T_su, P_high)
        h_hp_ex  = _expand(h_hp_su, P_high, P_rh,  eta_turb, fluid)
        h_lp_su  = _h(fluid, T_su, P_rh)
        h_ext    = _expand(h_lp_su, P_rh,   P_dea, eta_turb, fluid)
        h_lp_ex  = _expand(h_ext,   P_dea,  P_low, eta_turb, fluid)

        y = (h_fw_dea - h_plp_ex) / (h_ext - h_plp_ex)
        y = max(0.0, min(y, 0.5))

        W_turb_hp = 1.0     * (h_hp_su - h_hp_ex)
        W_turb_lp = 1.0     * (h_lp_su - h_ext) + (1.0-y) * (h_ext - h_lp_ex)
        W_pump_hp = 1.0     * (h_php_ex - h_fw_dea)
        W_pump_lp = (1.0-y) * (h_plp_ex - h_sat_liq_cond)
        W_net_sp  = W_turb_hp + W_turb_lp - W_pump_hp - W_pump_lp

        Q_boiler  = h_hp_su - h_php_ex
        Q_rh_sp   = h_lp_su - h_hp_ex
        Q_sg_sp   = Q_boiler + Q_rh_sp
        eta_sp    = W_net_sp / Q_sg_sp

        if h_lp_ex >= h_sat_vap_cond:
            x_ex = 1.0
        elif h_lp_ex <= h_sat_liq_cond:
            x_ex = 0.0
        else:
            x_ex = (h_lp_ex - h_sat_liq_cond) / (h_sat_vap_cond - h_sat_liq_cond)

    except Exception as e:
        print(f"    ⚠ bottom_steam exception : {e}")
        return None

    if W_net_sp <= 0:
        print(f"    ⚠ bottom_steam : W_net_sp={W_net_sp/1000:.2f} kJ/kg ≤ 0")
        return None
    if x_ex < 0.80:
        print(f"    ⚠ bottom_steam : x_ex={x_ex:.4f} < 0.80")
        return None

    m_dot_st = Q_available_kW * 1000 / Q_sg_sp
    W_net    = W_net_sp * m_dot_st / 1000
    Q_heater = Q_sg_sp  * m_dot_st / 1000

    try:
        T_sat_rh = PropsSI('T', 'P', P_rh, 'Q', 0.5, fluid)
    except Exception:
        T_sat_rh = 0.0

    return {
        'type':      BOTTOM_STEAM,
        'eta':       eta_sp,
        'W_net':     W_net,
        'Q_heater':  Q_heater,
        'm_dot':     m_dot_st,
        'fluid':     fluid,
        'P_high':    P_high,
        'P_rh':      P_rh,
        'P_dea':     P_dea,
        'T_su':      T_su,
        'T_sat_rh':  T_sat_rh,
        'x_ex':      x_ex,
        'y_dea':     y,
        'T_hot_in':  T_hot_K,
    }


# ── ORC — sélection automatique du fluide selon T_source ────────────────────

_ORC_FLUIDS = [
    (60,  120, 'R245fa',      'R245fa'),
    (120, 180, 'R1233zd(E)',  'R1233zd'),
    (180, 250, 'Cyclopentane','Cyclopentane'),
    (250, 350, 'Toluene',     'Toluène'),
    (350, 400, 'Water',       'Vapeur (basse P)'),
]

def select_orc_fluid(T_hot_K):
    """Retourne (fluid_CoolProp, label) optimal pour T_source donnée."""
    T_C = T_hot_K - 273.15
    for T_min, T_max, fluid, label in _ORC_FLUIDS:
        if T_min <= T_C < T_max:
            return fluid, label
    return 'Water', 'Vapeur (basse P)'


def bottom_orc(T_hot_K, Q_available_kW,
               P_nom_MW=None,
               eta_pump=0.71, eta_exp=0.75, eta_recup_orc=0.70,
               P_low_default=None,
               dT_pinch=10.0, dT_superheat=5.0,
               eta_gen=0.97):
    """
    ORC avec récupérateur interne + sélection automatique du fluide.

    Si P_nom_MW est fourni (1, 2 ou 5 MWel), les paramètres de performance
    sont chargés depuis ORC_MODELS via get_orc_params() — les kwargs
    eta_pump / eta_exp / eta_recup_orc / dT_pinch / dT_superheat
    sont alors ignorés.

    Layout :
      Pump → Recup(cold) → Evaporateur → Expandeur → Recup(hot) → Condenseur

    Parameters
    ----------
    T_hot_K         : température source chaude [K]
    Q_available_kW  : chaleur disponible [kW / kg_air/s]
    P_nom_MW        : taille nominale ORC [MWel] — 1, 2 ou 5
    eta_pump        : rendement pompe (ignoré si P_nom_MW fourni)
    eta_exp         : rendement expandeur (ignoré si P_nom_MW fourni)
    eta_recup_orc   : efficacité récupérateur (ignoré si P_nom_MW fourni)
    dT_pinch        : pinch point [K] (ignoré si P_nom_MW fourni)
    dT_superheat    : surchauffe [K]  (ignoré si P_nom_MW fourni)
    eta_gen         : rendement générateur [-] (défaut 0.97)
    """
    # ── Sélection des paramètres selon taille ─────────────────────────
    if P_nom_MW is not None:
        orc_params    = get_orc_params(P_nom_MW)
        eta_exp       = orc_params['eta_exp']
        eta_pump      = orc_params['eta_pump']
        eta_recup_orc = orc_params['eta_recup']
        dT_pinch      = orc_params['dT_pinch']
        dT_superheat  = orc_params['dT_sh']
        _size_label   = orc_params['label']
    else:
        _size_label   = "taille non spécifiée (paramètres manuels)"

    fluid, fluid_lbl = 'Toluene', 'Toluène'

    _P_low_defaults = {
        'R245fa':       1.8e5,
        'R1233zd(E)':   1.0e5,
        'Cyclopentane': 0.20e5,
        'Toluene':      0.04e5,
        'Water':        0.10e5,
    }
    P_low = P_low_default or _P_low_defaults.get(fluid, 1.0e5)

    # ── Pression haute ────────────────────────────────────────────────
    T_sat_max = T_hot_K - dT_pinch - dT_superheat
    try:
        T_crit_f   = PropsSI('Tcrit', '', 0, '', 0, fluid)
        P_crit_f   = PropsSI('Pcrit', '', 0, '', 0, fluid)
        T_sat_max  = min(T_sat_max, T_crit_f - 5.0)
        P_high_max = PropsSI('P', 'T', T_sat_max, 'Q', 1, fluid)
        P_high     = min(P_high_max, P_crit_f * 0.95)
    except Exception as e:
        print(f"    ⚠ bottom_orc P_high calc : {e}")
        return None
    P_high = max(P_high, P_low * 2)

    try:
        T_sat    = PropsSI('T', 'P', P_high, 'Q', 0.5, fluid)
        T_cond   = PropsSI('T', 'P', P_low,  'Q', 0,   fluid)

        h_pu_su  = PropsSI('H', 'P', P_low,  'Q', 0, fluid)
        h_pu_ex  = _pump(h_pu_su, P_low, P_high, eta_pump, fluid)

        h_exp_su = _h(fluid, T_sat + dT_superheat, P_high)
        h_exp_ex = _expand(h_exp_su, P_high, P_low, eta_exp, fluid)

        T_exp_ex = PropsSI('T', 'H', h_exp_ex, 'P', P_low, fluid)

        # ── Récupérateur interne ──────────────────────────────────────
        h_exp_ex_sat = PropsSI('H', 'P', P_low, 'Q', 1, fluid)
        Q_max_hot    = h_exp_ex - h_pu_su
        T_cold_max   = T_exp_ex - 5.0
        h_cold_max   = PropsSI('H', 'T', T_cold_max, 'P', P_high, fluid)
        Q_max_cold   = h_cold_max - h_pu_ex

        Q_recup  = eta_recup_orc * min(Q_max_hot, Q_max_cold)
        Q_recup  = max(0.0, Q_recup)

        h_recup_cold_ex = h_pu_ex  + Q_recup
        h_recup_hot_ex  = h_exp_ex - Q_recup

        # ── Bilan énergétique ─────────────────────────────────────────
        Q_sg_sp  = h_exp_su - h_recup_cold_ex
        W_net_sp = (h_exp_su - h_exp_ex) - (h_pu_ex - h_pu_su)

        # ── Rendement générateur ──────────────────────────────────────
        W_net_sp = W_net_sp * eta_gen
        eta_sp   = W_net_sp / Q_sg_sp

        if Q_sg_sp <= 0 or W_net_sp <= 0:
            print(f"    ⚠ bottom_orc ({fluid_lbl}) : Q_sg ou W_net ≤ 0")
            return None

        T_recup_cold_ex = PropsSI('T', 'H', h_recup_cold_ex, 'P', P_high, fluid)
        T_recup_hot_ex  = PropsSI('T', 'H', h_recup_hot_ex,  'P', P_low,  fluid)

    except Exception as e:
        print(f"    ⚠ bottom_orc ({fluid_lbl}) exception : {e}")
        return None

    m_dot_orc = Q_available_kW * 1000 / Q_sg_sp
    W_net     = W_net_sp * m_dot_orc / 1000
    Q_heater  = Q_sg_sp  * m_dot_orc / 1000

    return {
        'type':             BOTTOM_ORC,
        'eta':              eta_sp,
        'W_net':            W_net,
        'Q_heater':         Q_heater,
        'm_dot':            m_dot_orc,
        'fluid':            fluid,
        'fluid_lbl':        fluid_lbl,
        'P_high':           P_high,
        'P_low':            P_low,
        'T_sat':            T_sat,
        'T_exp_ex':         T_exp_ex,
        'T_recup_cold_ex':  T_recup_cold_ex,
        'T_recup_hot_ex':   T_recup_hot_ex,
        'Q_recup_sp':       Q_recup / 1000,
        'T_hot_in':         T_hot_K,
        'P_nom_MW':         P_nom_MW,
        'size_label':       _size_label,
        'eta_gen':          eta_gen,
    }


# ════════════════════════════════════════════════════════════════════════════
#   BOTTOM STEAM — ARCHITECTURE SIEMENS
# ════════════════════════════════════════════════════════════════════════════

def _ew(h_su, p_su, p_ex, eta):
    """Expansion isentropique Water."""
    s_su    = PropsSI('S', 'H', h_su, 'P', p_su, 'Water')
    h_ex_is = PropsSI('H', 'P', p_ex, 'S', s_su, 'Water')
    return h_su - eta * (h_su - h_ex_is)

def _pw(h_su, p_su, p_ex, eta):
    """Pompage isentropique Water."""
    s_su    = PropsSI('S', 'H', h_su, 'P', p_su, 'Water')
    h_ex_is = PropsSI('H', 'P', p_ex, 'S', s_su, 'Water')
    return h_su + (h_ex_is - h_su) / eta


def _build_siemens_circuit(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond_v,
                            T_sh_K, T_eva_K, T_eco_K, T_rh_K,
                            T_amb_K, m_dot_air, P_hot,
                            P_low, P_high, P_rh, m_dot_st, n_disc=10):
    """Construit le circuit LaboThapPy Siemens avec source Air."""
    T_sat_hi = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')

    def _disc():
        hx = HexCstEffDisc()
        hx.set_parameters(eta_max=eta_hx, Pinch_min=pinch_hx, n_disc=n_disc)
        return hx

    def _cond():
        hx = HexCstPinch()
        hx.set_parameters(Pinch=pinch_cond_v, Delta_T_sh_sc=2.0, HX_type='condenser')
        return hx

    def _air(T_K):
        c = MassConnector()
        c.set_properties(fluid='Air', T=T_K, p=P_hot, m_dot=m_dot_air)
        return c

    cycle   = RecursiveCircuit('Water')
    pump    = PumpCstEff();  pump.set_parameters(eta_is=eta_pump)
    eco     = _disc()
    eva     = _disc()
    sh      = _disc()
    rh      = _disc()
    turb_hp = ExpanderCstEff(); turb_hp.set_parameters(eta_is=eta_turb)
    turb_ic = ExpanderCstEff(); turb_ic.set_parameters(eta_is=eta_turb)
    cond    = _cond()

    for nm, cp in [("Pump", pump), ("Eco", eco), ("Eva", eva), ("SH", sh),
                   ("Turb_HP", turb_hp), ("RH", rh), ("Turb_IC", turb_ic), ("Cond", cond)]:
        cycle.add_component(cp, nm)

    cycle.link_components("Pump",    "m-ex",   "Eco",     "m-su_C")
    cycle.link_components("Eco",     "m-ex_C", "Eva",     "m-su_C")
    cycle.link_components("Eva",     "m-ex_C", "SH",      "m-su_C")
    cycle.link_components("SH",      "m-ex_C", "Turb_HP", "m-su")
    cycle.link_components("Turb_HP", "m-ex",   "RH",      "m-su_C")
    cycle.link_components("RH",      "m-ex_C", "Turb_IC", "m-su")
    cycle.link_components("Turb_IC", "m-ex",   "Cond",    "m-su_H")
    cycle.link_components("Cond",    "m-ex_H", "Pump",    "m-su")

    CSource_cond = MassConnector()
    CSource_cond.set_properties(fluid='Water', T=T_amb_K + 15, P=3e5, m_dot=100.0)

    cycle.add_source("Air_SH",  _air(T_sh_K),  cycle.components["SH"],   "m-su_H")
    cycle.add_source("Air_EVA", _air(T_eva_K), cycle.components["Eva"],  "m-su_H")
    cycle.add_source("Air_ECO", _air(T_eco_K), cycle.components["Eco"],  "m-su_H")
    cycle.add_source("Air_RH",  _air(T_rh_K),  cycle.components["RH"],   "m-su_H")
    cycle.add_source("ColdSrc", CSource_cond,  cycle.components["Cond"], "m-su_C")

    cycle.set_fixed_properties(target="Pump:ex",    p=P_high)
    cycle.set_fixed_properties(target="Turb_HP:ex", p=P_rh)
    cycle.set_fixed_properties(target="Turb_IC:ex", p=P_low)
    cycle.set_fixed_properties(target="Cond:ex_H",  p=P_low)

    T_su_g  = T_sh_K - pinch_hx
    h_sl    = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    h_pu_ex = _pw(h_sl, P_low, P_high, eta_pump)
    T_pu_ex = PropsSI('T', 'H', h_pu_ex, 'P', P_high, 'Water')
    h_HP_su = PropsSI('H', 'T', T_su_g,  'P', P_high, 'Water')
    h_HP_ex = _ew(h_HP_su, P_high, P_rh,  eta_turb)
    T_HP_ex = PropsSI('T', 'H', h_HP_ex, 'P', P_rh,   'Water')
    h_IC_su = PropsSI('H', 'T', T_su_g,  'P', P_rh,   'Water')
    h_IC_ex = _ew(h_IC_su, P_rh,   P_low, eta_turb)
    T_IC_ex = PropsSI('T', 'H', h_IC_ex, 'P', P_low,  'Water')

    cycle.set_cycle_guess(target='Pump:su',    p=P_low,  SC=3,      m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',    p=P_high, T=T_pu_ex)
    cycle.set_cycle_guess(target='Eco:ex_C',   p=P_high, T=T_sat_hi - 5)
    cycle.set_cycle_guess(target='Eva:ex_C',   p=P_high, T=T_sat_hi, x=1)
    cycle.set_cycle_guess(target='SH:ex_C',    p=P_high, T=T_su_g)
    cycle.set_cycle_guess(target='Turb_HP:su', p=P_high, T=T_su_g,  m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_HP:ex', p=P_rh,   T=T_HP_ex)
    cycle.set_cycle_guess(target='RH:ex_C',    p=P_rh,   T=T_su_g)
    cycle.set_cycle_guess(target='Turb_IC:su', p=P_rh,   T=T_su_g,  m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_IC:ex', p=P_low,  T=T_IC_ex)
    cycle.set_cycle_guess(target='Cond:ex_H',  p=P_low,  T=T_sat_lo - 2)

    for tgt, var in [
        ('Turb_HP:ex','h'), ('Turb_HP:ex','p'),
        ('Turb_IC:ex','h'), ('Turb_IC:ex','p'),
        ('Pump:ex',   'h'), ('Pump:ex',   'p'),
        ('SH:ex_C',   'h'), ('SH:ex_C',   'p'),
        ('RH:ex_C',   'h'), ('RH:ex_C',   'p'),
        ('Eva:ex_C',  'h'), ('Eva:ex_C',  'p'),
        ('Eco:ex_C',  'h'), ('Eco:ex_C',  'p'),
        ('Cond:ex_H', 'h'), ('Cond:ex_H', 'p'),
    ]:
        cycle.set_residual_variable(target=tgt, variable=var, tolerance=1e-6)

    return cycle, pump, eco, eva, sh, rh, turb_hp, turb_ic, cond


def _solve_siemens(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                   T_hot_K, T_amb_K, m_dot_air, P_hot,
                   P_low, P_high, P_rh, m_dot_st=1.0,
                   tol_outer=0.5, max_outer=15, n_disc=10):
    """Résolution itérative Siemens avec source Air."""
    T_sat_hi = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sh_K   = T_hot_K
    T_rh_K   = T_hot_K
    T_eva_K  = T_hot_K - 10.0
    T_eco_K  = T_sat_hi + 5.0

    for _ in range(max_outer):
        cycle, pump, eco, eva, sh, rh, turb_hp, turb_ic, cond = \
            _build_siemens_circuit(
                eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                T_sh_K, T_eva_K, T_eco_K, T_rh_K,
                T_amb_K, m_dot_air, P_hot,
                P_low, P_high, P_rh, m_dot_st, n_disc
            )
        cycle.solve()

        T_eva_new = sh.ex_H.T
        T_eco_new = eva.ex_H.T
        err = max(abs(T_eva_new - T_eva_K), abs(T_eco_new - T_eco_K))
        T_eva_K   = T_eva_new
        T_eco_K   = T_eco_new
        if err < tol_outer:
            break

    return pump, eco, eva, sh, rh, turb_hp, turb_ic, cond


def bottom_steam_siemens(T_hot_K, Q_available_kW, T_amb_K,
                          eta_pump=0.80, eta_turb=0.93, eta_hx=0.95,
                          pinch_hx=5.0, pinch_cond=5.0, pinch_cfwh=5.0,
                          P_high=160e5, P_rh=38.41e5, P_low=0.095e5,
                          P_hp_ext=40.50e5, P_fwh5=19.11e5, P_dea=11.64e5,
                          P_lp4=6.08e5, P_lp3=3.17e5, P_lp2=1.38e5, P_lp1=0.44e5):
    """Bottom cycle Rankine vapeur style Siemens."""
    T_crit = 647.096
    T_su   = T_hot_K - pinch_hx
    if T_su <= 273.15 + 200:
        return None

    T_sat_max = min(T_su - 5.0, T_crit - 5.0)
    try:
        P_high_max = PropsSI('P', 'T', T_sat_max, 'Q', 1, 'Water')
    except Exception:
        P_high_max = 220e5
    P_high = min(P_high, P_high_max)
    P_high = max(P_high, P_rh * 1.5)

    P_hp_ext = min(P_hp_ext, P_high  * 0.28)
    P_fwh5   = min(P_fwh5,   P_rh    * 0.50)
    P_dea    = min(P_dea,    P_fwh5  * 0.60)
    P_lp4    = min(P_lp4,    P_dea   * 0.52)
    P_lp3    = min(P_lp3,    P_lp4   * 0.52)
    P_lp2    = min(P_lp2,    P_lp3   * 0.43)
    P_lp1    = min(P_lp1,    P_lp2   * 0.32)
    P_lp1    = max(P_lp1,    P_low   * 3.0)

    P_hot     = 1.01325e5
    m_dot_air = 20.0
    m_dot_st  = 1.0

    try:
        pump, eco, eva, sh, rh, turb_hp, turb_ic, cond = \
            _solve_siemens(
                eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                T_hot_K, T_amb_K, m_dot_air, P_hot,
                P_low, P_high, P_rh, m_dot_st
            )
    except Exception as e:
        print(f"    ⚠ bottom_steam_siemens Pass1 : {e}")
        return None

    try:
        from steam_rankine import compute_fwh as _cfwh
        perf = _cfwh(pump, eco, eva, sh, rh, turb_hp, turb_ic, cond,
                     eta_pump, eta_turb, P_high, P_low,
                     P_hp_ext, P_fwh5, P_dea, P_lp4, P_lp3, P_lp2, P_lp1,
                     eta_g=1.0, pinch_cfwh=pinch_cfwh)
    except Exception as e:
        print(f"    ⚠ bottom_steam_siemens Pass2 : {e}")
        return None

    eta_th    = perf['eta_th']
    Q_boil_sp = perf['Q_boil']
    if Q_boil_sp <= 0 or eta_th <= 0:
        return None

    m_dot_st_real = Q_available_kW / Q_boil_sp
    W_net         = eta_th * Q_available_kW

    return {
        'type':        BOTTOM_STEAM,
        'eta':         eta_th,
        'W_net':       W_net,
        'Q_heater':    Q_available_kW,
        'm_dot':       m_dot_st_real,
        'fluid':       'Water',
        'P_high':      P_high,
        'P_rh':        P_rh,
        'T_su':        turb_hp.su.T,
        'x_ex':        1.0,
        'eta_el':      perf['eta_el'],
        'discrepancy': perf['discrepancy'],
    }


def available_heat(h_exhaust, T_return_min_K, P_air=1.01325e5):
    """Chaleur récupérable côté air [kW / kg_air/s]."""
    h_return = PropsSI('H', 'T', T_return_min_K, 'P', P_air, 'Air')
    return max(0.0, (h_exhaust - h_return) / 1000)


def T_return_for_bottom(bc_type, P_high_bottom, dT_pinch=10.0, fluid='Water'):
    """Température minimale de retour de l'air [K] selon le bottom cycle."""
    if bc_type == BOTTOM_STEAM:
        T_sat = PropsSI('T', 'P', P_high_bottom, 'Q', 0, 'Water')
        return T_sat - dT_pinch
    else:
        T_sat = PropsSI('T', 'P', 0.04e5, 'Q', 0, 'Toluene')
        return T_sat - dT_pinch


# ════════════════════════════════════════════════════════════════════════════
#   RUNNER PRINCIPAL
# ════════════════════════════════════════════════════════════════════════════

def run_combined(brayton_case, params, N_towers=20, verbose=True):
    """
    Lance le top cycle Brayton, auto-sélectionne le bottom cycle,
    intègre le transport thermique et retourne le bilan combiné.
    """
    T_amb  = params['T_amb']
    P_low  = params['P_low']
    P_high = params['P_high']
    P_mid  = params['P_mid']
    T_hot  = params['T_hot_su']
    P_hot  = params['P_hot']
    hot_fluid = params['hot_fluid']

    def _make_source(T, P, m=1.0):
        if hot_fluid == 'SolarSalt':
            c = SolarSaltConnector(); c.set_properties(T=T, p=P, m_dot=m)
        else:
            c = MassConnector(); c.set_properties(fluid='Air', T=T, p=P, m_dot=m)
        return c

    AirInlet = MassConnector()
    AirInlet.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=1.0)

    CSource = MassConnector()
    CSource.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=10.0)

    if brayton_case == "Simple":
        HSource = _make_source(T_hot, P_hot)
        cycle, comp, turb, heater = brayton_simple(
            params['eta_cp'], params['eta_tb'], params['eta_heater'],
            HSource, AirInlet, P_low, P_high
        )
        br = compute_brayton_perf(comp, turb, heater, hot_fluid, T_hot, P_hot)

    elif brayton_case == "Recuperated":
        HSource = _make_source(T_hot, P_hot)
        cycle, comp, turb, heater, recup = brayton_recuperated(
            params['eta_cp'], params['eta_tb'], params['eta_heater'],
            params['eta_recup'], HSource, AirInlet, P_low, P_high,
            T_tb_su_guess=T_hot
        )
        br = compute_brayton_perf(comp, turb, heater, hot_fluid, T_hot, P_hot,
                                  recuperator=recup)

    elif brayton_case == "Recuperated_RH_IC":
        HSource_h = _make_source(T_hot, P_hot)
        HSource_r = _make_source(T_hot, P_hot)
        cycle, comp1, comp2, turb1, turb2, heater, reheater, recup, ic = \
            brayton_recuperated_rh_ic(
                params['eta_cp'], params['eta_tb'], params['eta_heater'],
                params['eta_recup'], params['eta_cooler'],
                HSource_h, HSource_r, AirInlet, CSource,
                P_low, P_high, P_mid, T_tb_su_guess=T_hot
            )
        br = compute_brayton_perf(comp1, turb1, heater, hot_fluid, T_hot, P_hot,
                                  recuperator=recup, comp2=comp2, turb2=turb2,
                                  reheater=reheater, intercooler=ic)
    else:
        raise ValueError(f"brayton_case inconnu : {brayton_case}")

    tr = network_heat_loss(
        T_hot_K              = br['T_exhaust'],
        T_amb_K              = T_amb,
        m_dot_per_tower_kg_s = params['m_dot_transport_per_tower'],
        N                    = N_towers,
        verbose              = False,
    )
    T_bottom_in = tr['T_pb_in']

    bc_type = select_bottom_cycle(T_bottom_in)

    P_high_bottom_default = {'Steam': 40e5, 'ORC_R245fa': 20e5}.get(bc_type, 40e5)
    T_ret_min = T_return_for_bottom(bc_type, P_high_bottom_default)

    Q_avail = available_heat(br['h_exhaust'], T_ret_min)
    Q_avail = max(0.0, Q_avail)

    if bc_type == BOTTOM_STEAM:
        bp = bottom_steam_siemens(T_bottom_in, Q_avail, T_amb)
        if bp is None:
            print("    → Fallback sur bottom_steam simplifié")
            bp = bottom_steam(T_bottom_in, Q_avail)
    else:
        bp = bottom_orc(T_bottom_in, Q_avail,
                        P_nom_MW=params.get('orc_size_MW', None))

    if bp is not None:
        W_combined = br['W_net'] + bp['W_net']
        eta_comb   = W_combined / br['Q_total']
        gain_pp    = (eta_comb - br['eta']) * 100
    else:
        W_combined = br['W_net']
        eta_comb   = br['eta']
        gain_pp    = 0.0

    result = {
        'brayton_case': brayton_case,
        'br':           br,
        'tr':           tr,
        'bc_type':      bc_type,
        'bp':           bp,
        'T_bottom_in':  T_bottom_in,
        'Q_avail':      Q_avail,
        'W_combined':   W_combined,
        'eta_combined': eta_comb,
        'eta_gain_pp':  gain_pp,
    }

    if verbose:
        _print_combined(result, params)

    return result


# ════════════════════════════════════════════════════════════════════════════
#   AFFICHAGE
# ════════════════════════════════════════════════════════════════════════════

def _print_combined(r, params):
    br = r['br']
    tr = r['tr']
    bp = r['bp']

    print()
    print("╔" + "═"*68 + "╗")
    print(f"║  TOP CYCLE : Brayton {r['brayton_case']:<45}║")
    print(f"║  BOTTOM    : {bottom_label(r['bc_type'], bp):<55}║")
    print("╚" + "═"*68 + "╝")

    print("\n── Brayton (ref 1 kg_air/s) ────────────────────────────────────────")
    print(f"  η Brayton    : {br['eta']*100:.2f} %")
    print(f"  W_net        : {br['W_net']:.2f} kW")
    print(f"  Q_total      : {br['Q_total']:.2f} kW")
    print(f"  T_exhaust    : {br['T_exhaust']-273.15:.1f} °C")

    print("\n── Transport thermique ─────────────────────────────────────────────")
    print(f"  Fluide       : {tr['fluid_label']}")
    print(f"  L_réseau     : {tr['L_network']:.1f} m  |  L_moy = {tr['L_avg']:.1f} m  |  D : {tr['D_pipe']*1e3:.0f} mm")
    print(f"  T_PB_in      : {tr['T_pb_in']-273.15:.2f} °C  (ΔT = {tr['dT_loss']:.3f} K)")
    print(f"  Q_perte      : {tr['Q_loss_kW']:.3f} kW")

    print(f"\n── Sélection bottom cycle ──────────────────────────────────────────")
    print(f"  T_bottom_in  : {r['T_bottom_in']-273.15:.1f} °C")
    print(f"  Seuils       : Steam si T > {T_THRESH_STEAM-273.15:.0f} °C  |  "
          f"ORC si T ≤ {T_THRESH_STEAM-273.15:.0f} °C")
    print(f"  → Sélectionné : {bottom_label(r['bc_type'])}")
    print(f"  Q_disponible : {r['Q_avail']:.2f} kW / kg_air/s")

    print(f"\n── Bottom cycle ────────────────────────────────────────────────────")
    if bp is not None:
        if 'size_label' in bp:
            print(f"  Modèle taille: {bp['size_label']}")
        print(f"  η bottom     : {bp['eta']*100:.2f} %")
        print(f"  W_net        : {bp['W_net']:.2f} kW / kg_air/s")
        print(f"  m_dot_bottom : {bp['m_dot']:.4f} kg_{bp['fluid']}/s / kg_air/s")
        if 'P_high' in bp:
            p_low_bar = bp.get('P_low', 0) / 1e5
            print(f"  P_high       : {bp['P_high']/1e5:.1f} bar  |  P_low : {p_low_bar:.2f} bar")
        if 'T_sat' in bp:
            print(f"  T_sat_evap   : {bp['T_sat']-273.15:.1f} °C")
        if 'T_recup_cold_ex' in bp:
            print(f"  Récupérateur : T_cold_ex={bp['T_recup_cold_ex']-273.15:.1f}°C  "
                  f"T_hot_ex={bp['T_recup_hot_ex']-273.15:.1f}°C  "
                  f"Q={bp['Q_recup_sp']:.1f} kJ/kg")
        if 'P_rh' in bp:
            print(f"  P_reheat     : {bp['P_rh']/1e5:.0f} bar  "
                  f"(T_sat = {bp.get('T_sat_rh',273.15)-273.15:.1f} °C)")
        if 'P_dea' in bp:
            print(f"  P_DEA        : {bp['P_dea']/1e5:.1f} bar  "
                  f"(y_dea = {bp.get('y_dea',0):.4f})")
        if 'T_su' in bp:
            print(f"  T_turb_su    : {bp['T_su']-273.15:.1f} °C")
        if 'x_ex' in bp:
            print(f"  x_exp_ex     : {bp['x_ex']:.3f}")
    else:
        print("  ⚠ Bottom cycle : pas de solution convergente.")

    print(f"\n── Bilan combiné ───────────────────────────────────────────────────")
    print(f"  W_Brayton    : {br['W_net']:8.2f} kW / kg_air/s")
    print(f"  W_bottom     : {bp['W_net'] if bp else 0:8.2f} kW / kg_air/s")
    print(f"  W_combiné    : {r['W_combined']:8.2f} kW / kg_air/s")
    print(f"  η combiné    : {r['eta_combined']*100:.2f} %   "
          f"(η Brayton seul = {br['eta']*100:.2f} %,  gain = +{r['eta_gain_pp']:.2f} pp)")
    print("─"*70)


def print_comparison_table(results):
    """Tableau comparatif des configurations."""
    print()
    print("╔" + "═"*100 + "╗")
    print("║  TABLEAU COMPARATIF — Brayton progressif + bottom cycle auto-sélectionné"
          + " "*27 + "║")
    print("╠" + "═"*100 + "╣")
    hdr = (f"  {'Brayton':25s}  {'T_exh[°C]':10s}  {'η_br[%]':8s}  "
           f"{'Bottom':22s}  {'η_bot[%]':8s}  {'η_comb[%]':10s}  {'Gain[pp]':8s}")
    print("║" + hdr + " "*(100 - len(hdr)) + "║")
    print("╠" + "═"*100 + "╣")
    for r in results:
        br = r['br']
        bp = r['bp']
        eta_bot = bp['eta']*100 if bp else float('nan')
        row = (f"  {r['brayton_case']:25s}  {br['T_exhaust']-273.15:10.1f}  "
               f"{br['eta']*100:8.2f}  {bottom_label(r['bc_type'], bp):22s}  "
               f"{eta_bot:8.2f}  {r['eta_combined']*100:10.2f}  "
               f"{r['eta_gain_pp']:8.2f}")
        print("║" + row + " "*(100 - len(row)) + "║")
    print("╚" + "═"*100 + "╝")


# ════════════════════════════════════════════════════════════════════════════
#   MAIN
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    # ── Sweep variables ──────────────────────────────────────────────
    PR_values      = [3, 3.7567, 4.4373, 4.8853, 5]
    T_hot_values_C = np.arange(600, 1001, 50)
    ORC_SIZES_MW   = [1.0, 2.0, 5.0]   # tailles ORC à comparer

    # ── Paramètres fixes ─────────────────────────────────────────────
    T_amb        = 13.2 + 273.15
    P_low        = 1.01325e5
    T_salt_limit = 565 + 273.15
    P_hot        = 1e5
    N_towers     = 50

    base_params = dict(
        T_amb      = T_amb,
        P_low      = P_low,
        P_hot      = P_hot,
        eta_cp     = 0.75,
        eta_tb     = 0.85,
        eta_heater = 0.90,
        eta_recup  = 0.85,
        eta_cooler = 0.95,
        m_dot_transport_per_tower = 5.99,
    )

    # ── Sweep ────────────────────────────────────────────────────────
    results = []

    for PR in PR_values:
        P_high = P_low * PR
        P_mid  = np.sqrt(P_low * P_high)

        for T_hot_C in T_hot_values_C:
            T_hot_su  = T_hot_C + 273.15
            hot_fluid = 'SolarSalt' if T_hot_su <= T_salt_limit else 'Air'

            for orc_size in ORC_SIZES_MW:
                params = {
                    **base_params,
                    'T_hot_su':    T_hot_su,
                    'hot_fluid':   hot_fluid,
                    'P_high':      P_high,
                    'P_mid':       P_mid,
                    'orc_size_MW': orc_size,
                }

                try:
                    r  = run_combined("Recuperated", params,
                                      N_towers=N_towers, verbose=False)
                    br = r['br']
                    
                    # Passe 1 : obtenir T_recup_cold_ex
                    bp_pass1 = bottom_orc(r['T_bottom_in'], 999.0, P_nom_MW=orc_size)
                    
                    if bp_pass1 is not None:
                        T_ret_real   = bp_pass1['T_recup_cold_ex'] + 10.0
                        Q_avail_real = available_heat(br['h_exhaust'], T_ret_real)
                    else:
                        Q_avail_real = 0.0
                    
                    # Passe 2 : résoudre avec le bon Q_avail
                    bp = bottom_orc(r['T_bottom_in'], Q_avail_real, P_nom_MW=orc_size)

                    if bp is not None:
                        W_combined = br['W_net'] + bp['W_net']
                        eta_comb   = W_combined / br['Q_total']
                        gain_pp    = (eta_comb - br['eta']) * 100
                        bc_type    = BOTTOM_ORC
                        eta_orc    = bp['eta'] * 100
                    else:
                        W_combined = br['W_net']
                        eta_comb   = br['eta']
                        gain_pp    = 0.0
                        bc_type    = 'None'
                        eta_orc    = float('nan')

                    results.append({
                        'PR':           round(PR, 4),
                        'T_hot_C':      T_hot_C,
                        'hot_fluid':    hot_fluid,
                        'orc_size_MW':  orc_size,
                        'bc_type':      bc_type,
                        'orc_fluid':    bp.get('fluid_lbl', bp['fluid']) if bp else None,
                        'eta_br_%':     br['eta'] * 100,
                        'W_br_kW':      br['W_net'],
                        'T_exh_br_C':   br['T_exhaust'] - 273.15,
                        'T_bot_in_C':   r['T_bottom_in'] - 273.15,
                        'Q_avail_kW':   r['Q_avail'],
                        'eta_bot_%':    eta_orc,
                        'W_bot_kW':     bp['W_net'] if bp else np.nan,
                        'eta_comb_%':   eta_comb * 100,
                        'W_comb_kW':    W_combined,
                        'gain_pp':      gain_pp,
                    })

                    orc_lbl = bp.get('fluid_lbl', '?') if bp else 'None'
                    print(f"PR={PR:.4f} | T={T_hot_C}°C | {hot_fluid:9s} | "
                          f"ORC={orc_size:.0f}MW | {orc_lbl:12s} | "
                          f"η_br={br['eta']*100:5.2f}% | "
                          f"η_orc={eta_orc:5.2f}% | "
                          f"η_comb={eta_comb*100:5.2f}% | "
                          f"+{gain_pp:.2f}pp")

                except Exception as e:
                    print(f"PR={PR:.4f} | T={T_hot_C}°C | ORC={orc_size:.0f}MW | FAILED: {e}")
                    results.append({
                        'PR': round(PR, 4), 'T_hot_C': T_hot_C,
                        'hot_fluid': hot_fluid, 'orc_size_MW': orc_size,
                        'bc_type': None, 'orc_fluid': None,
                        'eta_br_%': np.nan, 'W_br_kW': np.nan,
                        'T_exh_br_C': np.nan, 'T_bot_in_C': np.nan,
                        'Q_avail_kW': np.nan, 'eta_bot_%': np.nan,
                        'W_bot_kW': np.nan, 'eta_comb_%': np.nan,
                        'W_comb_kW': np.nan, 'gain_pp': np.nan,
                    })

    # ── DataFrame & pivots ───────────────────────────────────────────
    df = pd.DataFrame(results)

    print("\n=== Résultats complets ===")
    print(df.to_string(index=False))

    # Pivot par taille ORC pour chaque métrique clé
    for col, label in [
        ('eta_br_%',   'η Brayton seul [%]'),
        ('eta_bot_%',  'η ORC [%]'),
        ('eta_comb_%', 'η combiné Brayton+ORC [%]'),
        ('gain_pp',    'Gain [pp]'),
        ('W_comb_kW',  'W_net combiné [kW / kg_air/s]'),
    ]:
        for orc_size in ORC_SIZES_MW:
            df_sub = df[df['orc_size_MW'] == orc_size]
            pivot  = df_sub.pivot(index='T_hot_C', columns='PR', values=col)
            print(f"\n=== {label} — ORC {orc_size:.0f} MW "
                  f"(lignes = T_hot [°C], colonnes = PR) ===")
            print(pivot.round(2).to_string())

    df.to_csv('sweep_recup_orc.csv', index=False)
    print("\nSauvegardé dans sweep_recup_orc.csv")