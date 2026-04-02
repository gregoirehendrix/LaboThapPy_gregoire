# -*- coding: utf-8 -*-
"""
Created on Apr 01 2026
Updated on Apr 02 2026

@author: gregoire.hendrix

Steam Rankine cycle — Solar Salt heat source

Available configurations (choose in __main__):
    study_case = "Simple"       # Pump → Boiler → Turbine → Condenser
    study_case = "Superheated"  # Pump → Eco → Eva → SH → Turbine → Condenser
    study_case = "Reheat"       # Pump → Eco → Eva → SH → Turb_HP → RH → Turb_LP → Condenser

HX components used:
    - HexCstEffDisc  : Boiler, Eco, Eva, SH, RH  (eta_max + Pinch_min enforced internally)
    - HexCstPinch    : Condenser  (HX_type='condenser', pinch-based solver with phase change)
    Both support SolarSalt as hot-side fluid natively.
"""

#%%
import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.expander.expander_csteff       import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc
from labothappy.component.heat_exchanger.hex_cstpinch    import HexCstPinch
from labothappy.component.pump.pump_csteff               import PumpCstEff


#%%
# =============================================================================
# SHARED HELPERS
# =============================================================================

def _salt_connector(T_K, p_Pa, m_dot):
    """Return a SolarSaltConnector fully initialised at (T, p, m_dot)."""
    c = SolarSaltConnector()
    c.set_properties(T=T_K, p=p_Pa, m_dot=m_dot)
    return c


def _make_disc_hx(eta_max=0.99, pinch_min=5.0, n_disc=10):
    """Return a configured HexCstEffDisc instance."""
    hx = HexCstEffDisc()
    hx.set_parameters(eta_max=eta_max, Pinch_min=pinch_min, n_disc=n_disc)
    return hx


def _make_pinch_cond(pinch=5.0, dT_sc=2.0):
    """Return a configured HexCstPinch condenser instance."""
    hx = HexCstPinch()
    hx.set_parameters(Pinch=pinch, Delta_T_sh_sc=dT_sc, HX_type='condenser')
    return hx


#%%
# =============================================================================
# CONFIG 1 — SIMPLE
# Pump → Boiler (HexCstEffDisc) → Turbine → Condenser (HexCstPinch) → Pump
# =============================================================================

def rankine_simple(eta_pump, eta_turb,
                   eta_boiler, pinch_boiler,
                   pinch_cond,
                   HSource, CSource,
                   P_low, P_high, m_dot_st,
                   T_turb_su_guess=None, T_turb_ex_guess=None,
                   n_disc=10):
    """
    Basic steam Rankine cycle.
    Pump → Boiler → Turbine_St → Condenser → Pump

    Parameters
    ----------
    eta_pump, eta_turb  : isentropic efficiencies [-]
    eta_boiler          : HexCstEffDisc eta_max for Boiler [-]
    pinch_boiler        : Pinch_min in Boiler [K]
    pinch_cond          : Pinch in Condenser [K]
    HSource             : SolarSaltConnector — hot source
    CSource             : MassConnector      — cooling water
    P_low, P_high       : cycle pressures [Pa]
    m_dot_st            : steam mass flow rate [kg/s]
    n_disc              : HexCstEffDisc discretisation segments
    """
    T_sat_hi = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')

    if T_turb_su_guess is None:
        T_turb_su_guess = T_sat_hi + 30
    if T_turb_ex_guess is None:
        T_turb_ex_guess = T_sat_lo + 5

    cycle     = RecursiveCircuit('Water')
    pump      = PumpCstEff()
    boiler    = _make_disc_hx(eta_max=eta_boiler, pinch_min=pinch_boiler, n_disc=n_disc)
    turbine   = ExpanderCstEff()
    condenser = _make_pinch_cond(pinch=pinch_cond)

    pump.set_parameters(eta_is=eta_pump)
    turbine.set_parameters(eta_is=eta_turb)

    cycle.add_component(pump,      "Pump")
    cycle.add_component(boiler,    "Boiler")
    cycle.add_component(turbine,   "Turbine_St")
    cycle.add_component(condenser, "Condenser")

    cycle.link_components("Pump",       "m-ex",   "Boiler",     "m-su_C")
    cycle.link_components("Boiler",     "m-ex_C", "Turbine_St", "m-su")
    cycle.link_components("Turbine_St", "m-ex",   "Condenser",  "m-su_H")
    cycle.link_components("Condenser",  "m-ex_H", "Pump",       "m-su")

    cycle.add_source("HotSource",  HSource, cycle.components["Boiler"],    "m-su_H")
    cycle.add_source("ColdSource", CSource, cycle.components["Condenser"], "m-su_C")

    cycle.set_fixed_properties(target="Pump:ex",       p=P_high)
    cycle.set_fixed_properties(target="Turbine_St:ex", p=P_low)

    # Analytical guesses
    h_sat_lo     = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    s_pump_su    = PropsSI('S', 'P', P_low,  'Q', 0, 'Water')
    h_pump_ex_is = PropsSI('H', 'P', P_high, 'S', s_pump_su, 'Water')
    h_pump_ex    = h_sat_lo + (h_pump_ex_is - h_sat_lo) / eta_pump
    T_pump_ex    = PropsSI('T', 'H', h_pump_ex, 'P', P_high, 'Water')

    s_turb_su    = PropsSI('S', 'T', T_turb_su_guess, 'P', P_high, 'Water')
    h_turb_su    = PropsSI('H', 'T', T_turb_su_guess, 'P', P_high, 'Water')
    h_turb_ex_is = PropsSI('H', 'P', P_low, 'S', s_turb_su, 'Water')
    h_turb_ex    = h_turb_su - eta_turb * (h_turb_su - h_turb_ex_is)
    T_turb_ex    = PropsSI('T', 'H', h_turb_ex, 'P', P_low, 'Water')

    cycle.set_cycle_guess(target='Pump:su',       p=P_low,  SC=3,              m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',       p=P_high, T=T_pump_ex)
    cycle.set_cycle_guess(target='Turbine_St:su', p=P_high, T=T_turb_su_guess, m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turbine_St:ex', p=P_low,  T=T_turb_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',p=P_low,  T=T_sat_lo - 2)

    for target, var in [
        ('Turbine_St:ex', 'h'), ('Turbine_St:ex', 'p'),
        ('Pump:ex',       'h'), ('Pump:ex',       'p'),
        ('Boiler:ex_C',   'h'), ('Boiler:ex_C',   'p'),
        ('Condenser:ex_H','h'), ('Condenser:ex_H','p'),
    ]:
        cycle.set_residual_variable(target=target, variable=var, tolerance=1e-6)

    return cycle, pump, boiler, turbine, condenser


def compute_cycle_performance(pump, boiler, turbine, condenser):
    m_dot = turbine.su.m_dot

    W_pump    = m_dot * (pump.ex.h      - pump.su.h)            / 1000
    W_turb    = m_dot * (turbine.su.h   - turbine.ex.h)         / 1000
    W_net     = W_turb - W_pump
    Q_boil    = m_dot * (boiler.ex_C.h  - boiler.su_C.h)        / 1000
    Q_cond    = m_dot * (condenser.su_H.h - condenser.ex_H.h)   / 1000
    eta_th    = W_net / Q_boil
    x_turb_ex = PropsSI('Q', 'H', turbine.ex.h, 'P', turbine.ex.p, 'Water')

    return {
        'W_pump':    W_pump,
        'W_turb':    W_turb,
        'W_net':     W_net,
        'Q_boiler':  Q_boil,
        'Q_cond':    Q_cond,
        'eta_th':    eta_th,
        'x_turb_ex': x_turb_ex,
    }


def print_states(pump, boiler, turbine, condenser):
    def _row(label, conn):
        print(f"  {label:<32} T={conn.T-273.15:7.2f} °C   p={conn.p/1e5:7.3f} bar   h={conn.h:10.2f} J/kg")
    print("=== Steam states ===")
    _row("Pump:su",               pump.su)
    _row("Pump:ex / Boiler:su_C", pump.ex)
    _row("Boiler:ex_C / Turb:su", boiler.ex_C)
    _row("Turbine:ex",            turbine.ex)
    _row("Condenser:ex_H",        condenser.ex_H)
    print("\n=== Salt state (hot side) ===")
    _row("Salt → Boiler:su_H",    boiler.su_H)
    _row("Salt → Boiler:ex_H",    boiler.ex_H)


def print_efficiency(perf):
    print("=== Simple Rankine — Performance ===")
    print(f"  W_pump      : {perf['W_pump']:.2f} kW/(kg/s)")
    print(f"  W_turb      : {perf['W_turb']:.2f} kW/(kg/s)")
    print(f"  W_net       : {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Q_boiler    : {perf['Q_boiler']:.2f} kW/(kg/s)")
    print(f"  Q_cond      : {perf['Q_cond']:.2f} kW/(kg/s)")
    print(f"  eta_th      : {perf['eta_th'] * 100:.2f} %")
    print(f"  x_turb_ex   : {perf['x_turb_ex']:.4f}  [-]")
    print("  --- First law check ---")
    balance = perf['Q_boiler'] - perf['Q_cond']
    print(f"  Q_boil - Q_cond : {balance:.4f} kW")
    print(f"  W_net           : {perf['W_net']:.4f} kW")
    print(f"  Discrepancy     : {abs(perf['W_net'] - balance):.6f} kW")
    print("=====================================")


def scale_to_power(W_net_target_MW, perf, PRINT_SCALE=True):
    scale = (W_net_target_MW * 1000) / perf['W_net']
    if PRINT_SCALE:
        print("=== Scaling to target net power ===")
        print(f"  W_net_target         : {W_net_target_MW:.2f} MW")
        print(f"  W_net_ref (1 kg/s)   : {perf['W_net']:.2f} kW")
        print(f"  Scale factor         : {scale:.4f}")
        print(f"  m_dot_steam required : {scale:.4f} kg/s")
        print(f"  W_pump               : {perf['W_pump']  * scale / 1000:.2f} MW")
        print(f"  W_turb               : {perf['W_turb']  * scale / 1000:.2f} MW")
        print(f"  W_net                : {perf['W_net']   * scale / 1000:.2f} MW")
        print(f"  Q_boiler             : {perf['Q_boiler']* scale / 1000:.2f} MW")
        print(f"  Q_cond               : {perf['Q_cond']  * scale / 1000:.2f} MW")
        print(f"  eta_th               : {perf['eta_th'] * 100:.2f} %")
        print("====================================")
    return scale


#%%
# =============================================================================
# CONFIG 2 — SUPERHEATED
# Pump → Eco → Eva → SH → Turbine → Condenser → Pump
# Salt series: HSource → SH(H) → Eva(H) → Eco(H)
# =============================================================================

def rankine_superheated(eta_pump, eta_turb,
                        eta_eco, eta_eva, eta_sh,
                        pinch_hx, pinch_cond,
                        HSource, CSource,
                        P_low, P_high, m_dot_st,
                        T_turb_ex_guess=None,
                        n_disc=10):
    """
    Superheated steam Rankine cycle.

    Steam side : Pump → Eco → Eva → SH → Turbine → Condenser → Pump
    Salt side  : HSource → SH(H) → Eva(H) → Eco(H)  (series counter-flow)

    HexCstEffDisc enforces Pinch_min internally on each HX — no workaround needed.

    Parameters
    ----------
    eta_pump, eta_turb      : isentropic efficiencies [-]
    eta_eco, eta_eva, eta_sh: eta_max for Eco, Eva, SH [-]
    pinch_hx                : Pinch_min applied to Eco, Eva, SH [K]
    pinch_cond              : Pinch for Condenser [K]
    HSource                 : SolarSaltConnector — hot salt (to SH)
    CSource                 : MassConnector      — cooling water
    P_low, P_high           : cycle pressures [Pa]
    m_dot_st                : steam mass flow rate [kg/s]
    T_turb_ex_guess         : turbine exit T guess [K]
    n_disc                  : discretisation segments
    """
    T_sat_hi   = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo   = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_salt_su  = HSource.T
    m_dot_salt = HSource.m_dot
    p_salt     = HSource.p

    if T_turb_ex_guess is None:
        T_turb_ex_guess = T_sat_lo + 5

    # Intermediate salt connectors (series: SH → Eva → Eco)
    Salt_eva_su = _salt_connector(T_sat_hi + 30, p_salt, m_dot_salt)
    Salt_eco_su = _salt_connector(T_sat_hi +  5, p_salt, m_dot_salt)

    cycle       = RecursiveCircuit('Water')
    pump        = PumpCstEff()
    eco         = _make_disc_hx(eta_max=eta_eco, pinch_min=pinch_hx, n_disc=n_disc)
    eva         = _make_disc_hx(eta_max=eta_eva, pinch_min=pinch_hx, n_disc=n_disc)
    superheater = _make_disc_hx(eta_max=eta_sh,  pinch_min=pinch_hx, n_disc=n_disc)
    turbine     = ExpanderCstEff()
    condenser   = _make_pinch_cond(pinch=pinch_cond)

    pump.set_parameters(eta_is=eta_pump)
    turbine.set_parameters(eta_is=eta_turb)

    cycle.add_component(pump,        "Pump")
    cycle.add_component(eco,         "Economiser")
    cycle.add_component(eva,         "Evaporator")
    cycle.add_component(superheater, "Superheater")
    cycle.add_component(turbine,     "Turbine_St")
    cycle.add_component(condenser,   "Condenser")

    cycle.link_components("Pump",        "m-ex",   "Economiser",  "m-su_C")
    cycle.link_components("Economiser",  "m-ex_C", "Evaporator",  "m-su_C")
    cycle.link_components("Evaporator",  "m-ex_C", "Superheater", "m-su_C")
    cycle.link_components("Superheater", "m-ex_C", "Turbine_St",  "m-su")
    cycle.link_components("Turbine_St",  "m-ex",   "Condenser",   "m-su_H")
    cycle.link_components("Condenser",   "m-ex_H", "Pump",        "m-su")

    cycle.add_source("SaltSource_SH",  HSource,     cycle.components["Superheater"], "m-su_H")
    cycle.add_source("SaltSource_EVA", Salt_eva_su, cycle.components["Evaporator"],  "m-su_H")
    cycle.add_source("SaltSource_ECO", Salt_eco_su, cycle.components["Economiser"],  "m-su_H")
    cycle.add_source("ColdSource",     CSource,     cycle.components["Condenser"],   "m-su_C")

    cycle.set_fixed_properties(target="Pump:ex",       p=P_high)
    cycle.set_fixed_properties(target="Turbine_St:ex", p=P_low)

    # Analytical guesses
    T_turb_su_guess = T_salt_su - pinch_hx
    h_sat_lo     = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    s_pump_su    = PropsSI('S', 'P', P_low,  'Q', 0, 'Water')
    h_pump_ex_is = PropsSI('H', 'P', P_high, 'S', s_pump_su, 'Water')
    h_pump_ex    = h_sat_lo + (h_pump_ex_is - h_sat_lo) / eta_pump
    T_pump_ex    = PropsSI('T', 'H', h_pump_ex, 'P', P_high, 'Water')

    s_turb_su    = PropsSI('S', 'T', T_turb_su_guess, 'P', P_high, 'Water')
    h_turb_su    = PropsSI('H', 'T', T_turb_su_guess, 'P', P_high, 'Water')
    h_turb_ex_is = PropsSI('H', 'P', P_low, 'S', s_turb_su, 'Water')
    h_turb_ex    = h_turb_su - eta_turb * (h_turb_su - h_turb_ex_is)
    T_turb_ex    = PropsSI('T', 'H', h_turb_ex, 'P', P_low, 'Water')

    cycle.set_cycle_guess(target='Pump:su',          p=P_low,  SC=3,              m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',          p=P_high, T=T_pump_ex)
    cycle.set_cycle_guess(target='Economiser:ex_C',  p=P_high, T=T_sat_hi - 5)
    cycle.set_cycle_guess(target='Evaporator:ex_C',  p=P_high, T=T_sat_hi,   x=1)
    cycle.set_cycle_guess(target='Superheater:ex_C', p=P_high, T=T_turb_su_guess)
    cycle.set_cycle_guess(target='Turbine_St:su',    p=P_high, T=T_turb_su_guess, m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turbine_St:ex',    p=P_low,  T=T_turb_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',   p=P_low,  T=T_sat_lo - 2)

    for target, var in [
        ('Turbine_St:ex',    'h'), ('Turbine_St:ex',    'p'),
        ('Pump:ex',          'h'), ('Pump:ex',          'p'),
        ('Superheater:ex_C', 'h'), ('Superheater:ex_C', 'p'),
        ('Evaporator:ex_C',  'h'), ('Evaporator:ex_C',  'p'),
        ('Economiser:ex_C',  'h'), ('Economiser:ex_C',  'p'),
        ('Condenser:ex_H',   'h'), ('Condenser:ex_H',   'p'),
    ]:
        cycle.set_residual_variable(target=target, variable=var, tolerance=1e-6)

    return cycle, pump, eco, eva, superheater, turbine, condenser


def compute_cycle_performance_sh(pump, eco, eva, superheater, turbine, condenser):
    m_dot = turbine.su.m_dot

    W_pump  = m_dot * (pump.ex.h          - pump.su.h)           / 1000
    W_turb  = m_dot * (turbine.su.h       - turbine.ex.h)        / 1000
    W_net   = W_turb - W_pump
    Q_eco   = m_dot * (eco.ex_C.h         - eco.su_C.h)          / 1000
    Q_eva   = m_dot * (eva.ex_C.h         - eva.su_C.h)          / 1000
    Q_sh    = m_dot * (superheater.ex_C.h - superheater.su_C.h)  / 1000
    Q_boil  = Q_eco + Q_eva + Q_sh
    Q_cond  = m_dot * (condenser.su_H.h   - condenser.ex_H.h)    / 1000
    eta_th  = W_net / Q_boil
    x_turb_ex = PropsSI('Q', 'H', turbine.ex.h, 'P', turbine.ex.p, 'Water')

    return {
        'W_pump':    W_pump,
        'W_turb':    W_turb,
        'W_net':     W_net,
        'Q_eco':     Q_eco,
        'Q_eva':     Q_eva,
        'Q_sh':      Q_sh,
        'Q_boiler':  Q_boil,
        'Q_cond':    Q_cond,
        'eta_th':    eta_th,
        'x_turb_ex': x_turb_ex,
    }


def print_states_sh(pump, eco, eva, superheater, turbine, condenser):
    def _row(label, conn):
        print(f"  {label:<32} T={conn.T-273.15:7.2f} °C   p={conn.p/1e5:7.3f} bar   h={conn.h:10.2f} J/kg")
    print("=== Steam states ===")
    _row("Pump:su",              pump.su)
    _row("Pump:ex / Eco:su_C",   pump.ex)
    _row("Eco:ex_C / Eva:su_C",  eco.ex_C)
    _row("Eva:ex_C / SH:su_C",   eva.ex_C)
    _row("SH:ex_C / Turb:su",    superheater.ex_C)
    _row("Turbine:ex",           turbine.ex)
    _row("Condenser:ex_H",       condenser.ex_H)
    print("\n=== Salt states (hot side) ===")
    _row("Salt → SH:su_H",       superheater.su_H)
    _row("Salt → SH:ex_H",       superheater.ex_H)
    _row("Salt → Eva:ex_H",      eva.ex_H)
    _row("Salt → Eco:ex_H",      eco.ex_H)


def print_efficiency_sh(perf):
    print("=== Superheated Rankine — Performance ===")
    print(f"  W_pump      : {perf['W_pump']:.2f} kW/(kg/s)")
    print(f"  W_turb      : {perf['W_turb']:.2f} kW/(kg/s)")
    print(f"  W_net       : {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Q_eco       : {perf['Q_eco']:.2f} kW/(kg/s)")
    print(f"  Q_eva       : {perf['Q_eva']:.2f} kW/(kg/s)")
    print(f"  Q_sh        : {perf['Q_sh']:.2f} kW/(kg/s)")
    print(f"  Q_boiler    : {perf['Q_boiler']:.2f} kW/(kg/s)  (eco+eva+sh)")
    print(f"  Q_cond      : {perf['Q_cond']:.2f} kW/(kg/s)")
    print(f"  eta_th      : {perf['eta_th'] * 100:.2f} %")
    print(f"  x_turb_ex   : {perf['x_turb_ex']:.4f}  [-]")
    print("  --- First law check ---")
    balance = perf['Q_boiler'] - perf['Q_cond']
    print(f"  Q_boil - Q_cond : {balance:.4f} kW")
    print(f"  W_net           : {perf['W_net']:.4f} kW")
    print(f"  Discrepancy     : {abs(perf['W_net'] - balance):.6f} kW")
    print("=========================================")


def scale_to_power_sh(W_net_target_MW, perf, PRINT_SCALE=True):
    scale = (W_net_target_MW * 1000) / perf['W_net']
    if PRINT_SCALE:
        print("=== Scaling to target net power ===")
        print(f"  W_net_target         : {W_net_target_MW:.2f} MW")
        print(f"  W_net_ref (1 kg/s)   : {perf['W_net']:.2f} kW")
        print(f"  Scale factor         : {scale:.4f}")
        print(f"  m_dot_steam required : {scale:.4f} kg/s")
        print(f"  W_pump               : {perf['W_pump']   * scale / 1000:.2f} MW")
        print(f"  W_turb               : {perf['W_turb']   * scale / 1000:.2f} MW")
        print(f"  W_net                : {perf['W_net']    * scale / 1000:.2f} MW")
        print(f"  Q_boiler             : {perf['Q_boiler'] * scale / 1000:.2f} MW")
        print(f"  Q_cond               : {perf['Q_cond']   * scale / 1000:.2f} MW")
        print(f"  eta_th               : {perf['eta_th'] * 100:.2f} %")
        print("====================================")
    return scale


#%%
# =============================================================================
# CONFIG 3 — REHEAT
# Pump → Eco → Eva → SH → Turb_HP → RH → Turb_LP → Condenser → Pump
# Salt series: HSource → SH(H) → RH(H) → Eva(H) → Eco(H)
# =============================================================================

def rankine_reheat(eta_pump, eta_turb,
                   eta_eco, eta_eva, eta_sh, eta_rh,
                   pinch_hx, pinch_cond,
                   HSource, CSource,
                   P_low, P_high, P_reheat, m_dot_st,
                   n_disc=10):
    """
    Reheat steam Rankine cycle.

    Steam side : Pump → Eco → Eva → SH → Turb_HP → RH → Turb_LP → Condenser → Pump
    Salt side  : HSource → SH(H) → RH(H) → Eva(H) → Eco(H)  (series counter-flow)

    All HX use HexCstEffDisc (pinch enforced) except Condenser (HexCstPinch).

    Parameters
    ----------
    eta_pump, eta_turb           : isentropic efficiencies [-]
    eta_eco, eta_eva, eta_sh, eta_rh : eta_max for each HX [-]
    pinch_hx                     : Pinch_min for Eco, Eva, SH, RH [K]
    pinch_cond                   : Pinch for Condenser [K]
    HSource                      : SolarSaltConnector — hot salt (to SH)
    CSource                      : MassConnector      — cooling water
    P_low, P_high                : LP and HP pressures [Pa]
    P_reheat                     : reheat pressure [Pa] (recommend √(P_high·P_low))
    m_dot_st                     : steam mass flow rate [kg/s]
    n_disc                       : discretisation segments
    """
    T_sat_hi   = PropsSI('T', 'P', P_high,   'Q', 1, 'Water')
    T_sat_lo   = PropsSI('T', 'P', P_low,    'Q', 0, 'Water')
    T_salt_su  = HSource.T
    m_dot_salt = HSource.m_dot
    p_salt     = HSource.p

    # Intermediate salt connectors (series: SH → RH → Eva → Eco)
    Salt_rh_su  = _salt_connector(T_sat_hi + 50, p_salt, m_dot_salt)
    Salt_eva_su = _salt_connector(T_sat_hi + 30, p_salt, m_dot_salt)
    Salt_eco_su = _salt_connector(T_sat_hi +  5, p_salt, m_dot_salt)

    cycle     = RecursiveCircuit('Water')
    pump      = PumpCstEff()
    eco       = _make_disc_hx(eta_max=eta_eco, pinch_min=pinch_hx, n_disc=n_disc)
    eva       = _make_disc_hx(eta_max=eta_eva, pinch_min=pinch_hx, n_disc=n_disc)
    sh        = _make_disc_hx(eta_max=eta_sh,  pinch_min=pinch_hx, n_disc=n_disc)
    rh        = _make_disc_hx(eta_max=eta_rh,  pinch_min=pinch_hx, n_disc=n_disc)
    turb_hp   = ExpanderCstEff()
    turb_lp   = ExpanderCstEff()
    condenser = _make_pinch_cond(pinch=pinch_cond)

    pump.set_parameters(eta_is=eta_pump)
    turb_hp.set_parameters(eta_is=eta_turb)
    turb_lp.set_parameters(eta_is=eta_turb)

    cycle.add_component(pump,      "Pump")
    cycle.add_component(eco,       "Economiser")
    cycle.add_component(eva,       "Evaporator")
    cycle.add_component(sh,        "Superheater")
    cycle.add_component(turb_hp,   "Turb_HP")
    cycle.add_component(rh,        "Reheater")
    cycle.add_component(turb_lp,   "Turb_LP")
    cycle.add_component(condenser, "Condenser")

    cycle.link_components("Pump",        "m-ex",   "Economiser",  "m-su_C")
    cycle.link_components("Economiser",  "m-ex_C", "Evaporator",  "m-su_C")
    cycle.link_components("Evaporator",  "m-ex_C", "Superheater", "m-su_C")
    cycle.link_components("Superheater", "m-ex_C", "Turb_HP",     "m-su")
    cycle.link_components("Turb_HP",     "m-ex",   "Reheater",    "m-su_C")
    cycle.link_components("Reheater",    "m-ex_C", "Turb_LP",     "m-su")
    cycle.link_components("Turb_LP",     "m-ex",   "Condenser",   "m-su_H")
    cycle.link_components("Condenser",   "m-ex_H", "Pump",        "m-su")

    cycle.add_source("SaltSource_SH",  HSource,     cycle.components["Superheater"], "m-su_H")
    cycle.add_source("SaltSource_RH",  Salt_rh_su,  cycle.components["Reheater"],    "m-su_H")
    cycle.add_source("SaltSource_EVA", Salt_eva_su, cycle.components["Evaporator"],  "m-su_H")
    cycle.add_source("SaltSource_ECO", Salt_eco_su, cycle.components["Economiser"],  "m-su_H")
    cycle.add_source("ColdSource",     CSource,     cycle.components["Condenser"],   "m-su_C")

    cycle.set_fixed_properties(target="Pump:ex",    p=P_high)
    cycle.set_fixed_properties(target="Turb_HP:ex", p=P_reheat)
    cycle.set_fixed_properties(target="Turb_LP:ex", p=P_low)

    # Analytical guesses
    T_HP_su_g = T_salt_su - pinch_hx
    T_LP_su_g = T_salt_su - pinch_hx

    h_sat_lo     = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    s_pump_su    = PropsSI('S', 'P', P_low,  'Q', 0, 'Water')
    h_pump_ex_is = PropsSI('H', 'P', P_high, 'S', s_pump_su, 'Water')
    h_pump_ex    = h_sat_lo + (h_pump_ex_is - h_sat_lo) / eta_pump
    T_pump_ex    = PropsSI('T', 'H', h_pump_ex, 'P', P_high, 'Water')

    s_HP_su    = PropsSI('S', 'T', T_HP_su_g, 'P', P_high, 'Water')
    h_HP_su    = PropsSI('H', 'T', T_HP_su_g, 'P', P_high, 'Water')
    h_HP_ex_is = PropsSI('H', 'P', P_reheat, 'S', s_HP_su, 'Water')
    h_HP_ex    = h_HP_su - eta_turb * (h_HP_su - h_HP_ex_is)
    T_HP_ex    = PropsSI('T', 'H', h_HP_ex, 'P', P_reheat, 'Water')

    s_LP_su    = PropsSI('S', 'T', T_LP_su_g, 'P', P_reheat, 'Water')
    h_LP_su    = PropsSI('H', 'T', T_LP_su_g, 'P', P_reheat, 'Water')
    h_LP_ex_is = PropsSI('H', 'P', P_low, 'S', s_LP_su, 'Water')
    h_LP_ex    = h_LP_su - eta_turb * (h_LP_su - h_LP_ex_is)
    T_LP_ex    = PropsSI('T', 'H', h_LP_ex, 'P', P_low, 'Water')

    cycle.set_cycle_guess(target='Pump:su',          p=P_low,    SC=3,           m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',          p=P_high,   T=T_pump_ex)
    cycle.set_cycle_guess(target='Economiser:ex_C',  p=P_high,   T=T_sat_hi - 5)
    cycle.set_cycle_guess(target='Evaporator:ex_C',  p=P_high,   T=T_sat_hi, x=1)
    cycle.set_cycle_guess(target='Superheater:ex_C', p=P_high,   T=T_HP_su_g)
    cycle.set_cycle_guess(target='Turb_HP:su',       p=P_high,   T=T_HP_su_g,   m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_HP:ex',       p=P_reheat, T=T_HP_ex)
    cycle.set_cycle_guess(target='Reheater:ex_C',    p=P_reheat, T=T_LP_su_g)
    cycle.set_cycle_guess(target='Turb_LP:su',       p=P_reheat, T=T_LP_su_g,   m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_LP:ex',       p=P_low,    T=T_LP_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',   p=P_low,    T=T_sat_lo - 2)

    for target, var in [
        ('Turb_HP:ex',       'h'), ('Turb_HP:ex',       'p'),
        ('Turb_LP:ex',       'h'), ('Turb_LP:ex',       'p'),
        ('Pump:ex',          'h'), ('Pump:ex',          'p'),
        ('Superheater:ex_C', 'h'), ('Superheater:ex_C', 'p'),
        ('Reheater:ex_C',    'h'), ('Reheater:ex_C',    'p'),
        ('Evaporator:ex_C',  'h'), ('Evaporator:ex_C',  'p'),
        ('Economiser:ex_C',  'h'), ('Economiser:ex_C',  'p'),
        ('Condenser:ex_H',   'h'), ('Condenser:ex_H',   'p'),
    ]:
        cycle.set_residual_variable(target=target, variable=var, tolerance=1e-6)

    return cycle, pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser


def compute_cycle_performance_rh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser):
    m_dot = turb_hp.su.m_dot

    W_pump  = m_dot * (pump.ex.h    - pump.su.h)              / 1000
    W_HP    = m_dot * (turb_hp.su.h - turb_hp.ex.h)           / 1000
    W_LP    = m_dot * (turb_lp.su.h - turb_lp.ex.h)           / 1000
    W_turb  = W_HP + W_LP
    W_net   = W_turb - W_pump
    Q_eco   = m_dot * (eco.ex_C.h   - eco.su_C.h)             / 1000
    Q_eva   = m_dot * (eva.ex_C.h   - eva.su_C.h)             / 1000
    Q_sh    = m_dot * (sh.ex_C.h    - sh.su_C.h)              / 1000
    Q_rh    = m_dot * (rh.ex_C.h    - rh.su_C.h)              / 1000
    Q_boil  = Q_eco + Q_eva + Q_sh + Q_rh
    Q_cond  = m_dot * (condenser.su_H.h - condenser.ex_H.h)   / 1000
    eta_th  = W_net / Q_boil
    x_LP_ex = PropsSI('Q', 'H', turb_lp.ex.h, 'P', turb_lp.ex.p, 'Water')

    return {
        'W_pump':   W_pump,
        'W_HP':     W_HP,
        'W_LP':     W_LP,
        'W_turb':   W_turb,
        'W_net':    W_net,
        'Q_eco':    Q_eco,
        'Q_eva':    Q_eva,
        'Q_sh':     Q_sh,
        'Q_rh':     Q_rh,
        'Q_boiler': Q_boil,
        'Q_cond':   Q_cond,
        'eta_th':   eta_th,
        'x_LP_ex':  x_LP_ex,
    }


def print_states_rh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser):
    def _row(label, conn):
        print(f"  {label:<35} T={conn.T-273.15:7.2f} °C   p={conn.p/1e5:7.3f} bar   h={conn.h:10.2f} J/kg")
    print("=== Steam states ===")
    _row("Pump:su",               pump.su)
    _row("Pump:ex / Eco:su_C",    pump.ex)
    _row("Eco:ex_C / Eva:su_C",   eco.ex_C)
    _row("Eva:ex_C / SH:su_C",    eva.ex_C)
    _row("SH:ex_C / Turb_HP:su",  sh.ex_C)
    _row("Turb_HP:ex / RH:su_C",  turb_hp.ex)
    _row("RH:ex_C / Turb_LP:su",  rh.ex_C)
    _row("Turb_LP:ex",            turb_lp.ex)
    _row("Condenser:ex_H",        condenser.ex_H)
    print("\n=== Salt states (hot side) ===")
    _row("Salt → SH:su_H",        sh.su_H)
    _row("Salt → SH:ex_H",        sh.ex_H)
    _row("Salt → RH:ex_H",        rh.ex_H)
    _row("Salt → Eva:ex_H",       eva.ex_H)
    _row("Salt → Eco:ex_H",       eco.ex_H)


def print_efficiency_rh(perf):
    print("=== Reheat Rankine — Performance ===")
    print(f"  W_pump      : {perf['W_pump']:.2f} kW/(kg/s)")
    print(f"  W_HP        : {perf['W_HP']:.2f} kW/(kg/s)")
    print(f"  W_LP        : {perf['W_LP']:.2f} kW/(kg/s)")
    print(f"  W_turb      : {perf['W_turb']:.2f} kW/(kg/s)")
    print(f"  W_net       : {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Q_eco       : {perf['Q_eco']:.2f} kW/(kg/s)")
    print(f"  Q_eva       : {perf['Q_eva']:.2f} kW/(kg/s)")
    print(f"  Q_sh        : {perf['Q_sh']:.2f} kW/(kg/s)")
    print(f"  Q_rh        : {perf['Q_rh']:.2f} kW/(kg/s)")
    print(f"  Q_boiler    : {perf['Q_boiler']:.2f} kW/(kg/s)  (eco+eva+sh+rh)")
    print(f"  Q_cond      : {perf['Q_cond']:.2f} kW/(kg/s)")
    print(f"  eta_th      : {perf['eta_th'] * 100:.2f} %")
    print(f"  x_LP_ex     : {perf['x_LP_ex']:.4f}  [-]  (steam quality at LP turbine exit)")
    print("  --- First law check ---")
    balance = perf['Q_boiler'] - perf['Q_cond']
    print(f"  Q_boil - Q_cond : {balance:.4f} kW")
    print(f"  W_net           : {perf['W_net']:.4f} kW")
    print(f"  Discrepancy     : {abs(perf['W_net'] - balance):.6f} kW")
    print("=====================================")


def scale_to_power_rh(W_net_target_MW, perf, PRINT_SCALE=True):
    scale = (W_net_target_MW * 1000) / perf['W_net']
    if PRINT_SCALE:
        print("=== Scaling to target net power ===")
        print(f"  W_net_target         : {W_net_target_MW:.2f} MW")
        print(f"  W_net_ref (1 kg/s)   : {perf['W_net']:.2f} kW")
        print(f"  Scale factor         : {scale:.4f}")
        print(f"  m_dot_steam required : {scale:.4f} kg/s")
        print(f"  W_pump               : {perf['W_pump']  * scale / 1000:.2f} MW")
        print(f"  W_HP                 : {perf['W_HP']    * scale / 1000:.2f} MW")
        print(f"  W_LP                 : {perf['W_LP']    * scale / 1000:.2f} MW")
        print(f"  W_net                : {perf['W_net']   * scale / 1000:.2f} MW")
        print(f"  Q_boiler             : {perf['Q_boiler']* scale / 1000:.2f} MW")
        print(f"  Q_cond               : {perf['Q_cond']  * scale / 1000:.2f} MW")
        print(f"  eta_th               : {perf['eta_th'] * 100:.2f} %")
        print("====================================")
    return scale


#%%
# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    study_case   = "Superheated"   # "Simple" | "Superheated" | "Reheat"
    PRINT        = True
    DETAIL       = True
    PRINT_SCALE  = True
    W_net_target = 100        # [MW]

    # --- Component efficiencies ---
    eta_pump = 0.75
    eta_turb = 0.85

    # --- HX parameters ---
    eta_hx     = 0.95   # eta_max for all HexCstEffDisc
    pinch_hx   = 5.0    # Pinch_min for Eco, Eva, SH, RH  [K]
    pinch_cond = 5.0    # Pinch for Condenser              [K]

    # --- Cycle pressures ---
    P_high = 160e5     # Pa
    P_low  = 0.10e5    # Pa

    # --- Steam mass flow rate (specific basis for scaling) ---
    m_dot_st = 1.0     # kg/s

    # --- Hot source: Solar Salt ---
    T_salt_su = 565 + 273.15   # K
    P_salt    = 2e5            # Pa

    HSource = SolarSaltConnector()
    HSource.set_properties(T=T_salt_su, p=P_salt, m_dot=50.0)
    print(f"Hot source : Solar Salt at {T_salt_su - 273.15:.1f} °C")

    # --- Cold source: cooling water ---
    T_cw_su = 20 + 273.15
    P_cw    = 3e5

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=T_cw_su, P=P_cw, m_dot=100.0)

    # -------------------------------------------------------------------------
    if study_case == "Simple":
        cycle, pump, boiler, turbine, condenser = \
            rankine_simple(
                eta_pump, eta_turb,
                eta_hx, pinch_hx, pinch_cond,
                HSource, CSource,
                P_low, P_high, m_dot_st,
            )
        cycle.plot_flag = False
        cycle.solve()
        perf = compute_cycle_performance(pump, boiler, turbine, condenser)
        if PRINT:
            if DETAIL: print_states(pump, boiler, turbine, condenser)
            print_efficiency(perf)
        scale = scale_to_power(W_net_target, perf, PRINT_SCALE)

    # -------------------------------------------------------------------------
    elif study_case == "Superheated":
        cycle, pump, eco, eva, superheater, turbine, condenser = \
            rankine_superheated(
                eta_pump, eta_turb,
                eta_hx, eta_hx, eta_hx,
                pinch_hx, pinch_cond,
                HSource, CSource,
                P_low, P_high, m_dot_st,
            )
        cycle.plot_flag = False
        cycle.solve()
        perf = compute_cycle_performance_sh(pump, eco, eva, superheater, turbine, condenser)
        if PRINT:
            if DETAIL: print_states_sh(pump, eco, eva, superheater, turbine, condenser)
            print_efficiency_sh(perf)
        scale = scale_to_power_sh(W_net_target, perf, PRINT_SCALE)

    # -------------------------------------------------------------------------
    elif study_case == "Reheat":
        P_reheat = np.sqrt(P_high * P_low)
        print(f"P_reheat (geometric mean) : {P_reheat/1e5:.2f} bar")
        cycle, pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser = \
            rankine_reheat(
                eta_pump, eta_turb,
                eta_hx, eta_hx, eta_hx, eta_hx,
                pinch_hx, pinch_cond,
                HSource, CSource,
                P_low, P_high, P_reheat, m_dot_st,
            )
        cycle.plot_flag = False
        cycle.solve()
        perf = compute_cycle_performance_rh(pump, eco, eva, sh, rh,
                                            turb_hp, turb_lp, condenser)
        if PRINT:
            if DETAIL: print_states_rh(pump, eco, eva, sh, rh,
                                       turb_hp, turb_lp, condenser)
            print_efficiency_rh(perf)
        scale = scale_to_power_rh(W_net_target, perf, PRINT_SCALE)