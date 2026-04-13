# -*- coding: utf-8 -*-
"""
Created on Apr 01 2026
Updated on Apr 02 2026

@author: gregoire.hendrix

Steam Rankine cycle — Solar Salt heat source

Available configurations (choose in __main__):
    study_case = "Simple"       # Pump → Boiler → Turbine → Condenser
    study_case = "Superheated"  # Pump → Eco → Eva → SH → Turbine → Condenser
                                #   salt series: SH → Eva → Eco  (outer convergence loop)
    study_case = "Reheat"       # Pump → Eco → Eva → SH → Turb_HP → RH → Turb_LP → Condenser
                                #   salt split: (1-f)*m_salt → SH → Eva → Eco (outer loop)
                                #               f*m_salt     → RH (fresh salt, parallel)
                                #   2D sweep on f_rh × P_reheat

HX components:
    - HexCstEffDisc : Boiler, Eco, Eva, SH, RH  (eta_max + Pinch_min enforced internally)
    - HexCstPinch   : Condenser  (HX_type='condenser', pinch-based with phase change)

Salt series fix:
    RecursiveCircuit treats each add_source as an independent fixed boundary condition.
    For series salt routing (SH → Eva → Eco), the intermediate connector temperatures
    must be iterated externally until sh.ex_H.T → Salt_eva_su.T and
    eva.ex_H.T → Salt_eco_su.T converge.  The outer loop below implements this.
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.expander.expander_csteff       import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc
from labothappy.component.heat_exchanger.hex_cstpinch    import HexCstPinch
from labothappy.component.pump.pump_csteff               import PumpCstEff

import warnings
warnings.filterwarnings('ignore')
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
# Single salt source — no outer loop needed.
# =============================================================================

def rankine_simple(eta_pump, eta_turb,
                   eta_boiler, pinch_boiler, pinch_cond,
                   HSource, CSource,
                   P_low, P_high, m_dot_st,
                   T_turb_su_guess=None, T_turb_ex_guess=None,
                   n_disc=10):
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

    cycle.set_cycle_guess(target='Pump:su',        p=P_low,  SC=3,               m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',        p=P_high, T=T_pump_ex)
    cycle.set_cycle_guess(target='Turbine_St:su',  p=P_high, T=T_turb_su_guess,  m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turbine_St:ex',  p=P_low,  T=T_turb_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H', p=P_low,  T=T_sat_lo - 2)

    for target, var in [
        ('Turbine_St:ex', 'h'), ('Turbine_St:ex', 'p'),
        ('Pump:ex',       'h'), ('Pump:ex',       'p'),
        ('Boiler:ex_C',   'h'), ('Boiler:ex_C',   'p'),
        ('Condenser:ex_H','h'), ('Condenser:ex_H','p'),
    ]:
        cycle.set_residual_variable(target=target, variable=var, tolerance=1e-6)

    return cycle, pump, boiler, turbine, condenser


def compute_cycle_performance(pump, boiler, turbine, condenser):
    m_dot     = turbine.su.m_dot
    W_pump    = m_dot * (pump.ex.h        - pump.su.h)          / 1000
    W_turb    = m_dot * (turbine.su.h     - turbine.ex.h)       / 1000
    W_net     = W_turb - W_pump
    Q_boil    = m_dot * (boiler.ex_C.h    - boiler.su_C.h)      / 1000
    Q_cond    = m_dot * (condenser.su_H.h - condenser.ex_H.h)   / 1000
    eta_th    = W_net / Q_boil
    x_turb_ex = PropsSI('Q', 'H', turbine.ex.h, 'P', turbine.ex.p, 'Water')
    return {
        'W_pump': W_pump, 'W_turb': W_turb, 'W_net': W_net,
        'Q_boiler': Q_boil, 'Q_cond': Q_cond,
        'eta_th': eta_th, 'x_turb_ex': x_turb_ex,
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
# CONFIG 2 — SUPERHEATED  (with outer salt-series convergence loop)
# Pump → Eco → Eva → SH → Turbine → Condenser → Pump
# Salt series (counter-flow): HSource(565°C) → SH → Eva → Eco
#
# Outer loop:
#   After each inner solve, update:
#     Salt_eva_su.T ← superheater.ex_H.T   (salt leaving SH)
#     Salt_eco_su.T ← eva.ex_H.T           (salt leaving Eva)
#   Rebuild and re-solve until max(ΔT) < tol_outer [K].
# =============================================================================

def _build_superheated_circuit(eta_pump, eta_turb,
                                eta_eco, eta_eva, eta_sh,
                                pinch_hx, pinch_cond,
                                HSource, Salt_eva_su, Salt_eco_su, CSource,
                                P_low, P_high, m_dot_st, n_disc):
    """Build and return a fresh RecursiveCircuit for the superheated cycle."""
    T_sat_hi  = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo  = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_salt_su = HSource.T

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

    cycle.set_cycle_guess(target='Pump:su',          p=P_low,  SC=3,               m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',          p=P_high, T=T_pump_ex)
    cycle.set_cycle_guess(target='Economiser:ex_C',  p=P_high, T=T_sat_hi - 5)
    cycle.set_cycle_guess(target='Evaporator:ex_C',  p=P_high, T=T_sat_hi,    x=1)
    cycle.set_cycle_guess(target='Superheater:ex_C', p=P_high, T=T_turb_su_guess)
    cycle.set_cycle_guess(target='Turbine_St:su',    p=P_high, T=T_turb_su_guess,  m_dot=m_dot_st)
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


def rankine_superheated(eta_pump, eta_turb,
                        eta_eco, eta_eva, eta_sh,
                        pinch_hx, pinch_cond,
                        HSource, CSource,
                        P_low, P_high, m_dot_st,
                        n_disc=10, tol_outer=0.1, max_outer=20):
    """
    Superheated Rankine with outer salt-series convergence loop.

    Intermediate salt temperatures (SH exit → Eva inlet, Eva exit → Eco inlet)
    are iterated until convergence, ensuring physically consistent series HT.
    """
    T_sat_hi   = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    p_salt     = HSource.p
    m_dot_salt = HSource.m_dot

    T_eva_su = HSource.T - 10.0   # initial guess: just below SH inlet
    T_eco_su = T_sat_hi + 5.0     # initial guess: just above saturation

    for outer_iter in range(max_outer):
        Salt_eva_su = _salt_connector(T_eva_su, p_salt, m_dot_salt)
        Salt_eco_su = _salt_connector(T_eco_su, p_salt, m_dot_salt)

        cycle, pump, eco, eva, superheater, turbine, condenser = \
            _build_superheated_circuit(
                eta_pump, eta_turb, eta_eco, eta_eva, eta_sh,
                pinch_hx, pinch_cond,
                HSource, Salt_eva_su, Salt_eco_su, CSource,
                P_low, P_high, m_dot_st, n_disc,
            )
        cycle.plot_flag = False
        cycle.solve()

        T_eva_su_new = superheater.ex_H.T   # salt leaving SH → entering Eva
        T_eco_su_new = eva.ex_H.T           # salt leaving Eva → entering Eco

        err = max(abs(T_eva_su_new - T_eva_su), abs(T_eco_su_new - T_eco_su))
        T_eva_su = T_eva_su_new
        T_eco_su = T_eco_su_new

        if err < tol_outer:
            print(f"  [SH outer loop] converged in {outer_iter + 1} iteration(s)  "
                  f"(err={err:.4f} K)")
            break
    else:
        print(f"  [SH outer loop] WARNING: did not converge after {max_outer} iterations "
              f"(err={err:.4f} K)")

    return cycle, pump, eco, eva, superheater, turbine, condenser


def compute_cycle_performance_sh(pump, eco, eva, superheater, turbine, condenser):
    m_dot   = turbine.su.m_dot
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
        'W_pump': W_pump, 'W_turb': W_turb, 'W_net': W_net,
        'Q_eco': Q_eco, 'Q_eva': Q_eva, 'Q_sh': Q_sh,
        'Q_boiler': Q_boil, 'Q_cond': Q_cond,
        'eta_th': eta_th, 'x_turb_ex': x_turb_ex,
    }


def print_states_sh(pump, eco, eva, superheater, turbine, condenser):
    def _row(label, conn):
        print(f"  {label:<38} T={conn.T-273.15:7.2f} °C   p={conn.p/1e5:7.3f} bar   h={conn.h:10.2f} J/kg")
    print("=== Steam states ===")
    _row("Pump:su",              pump.su)
    _row("Pump:ex / Eco:su_C",   pump.ex)
    _row("Eco:ex_C / Eva:su_C",  eco.ex_C)
    _row("Eva:ex_C / SH:su_C",   eva.ex_C)
    _row("SH:ex_C / Turb:su",    superheater.ex_C)
    _row("Turbine:ex",           turbine.ex)
    _row("Condenser:ex_H",       condenser.ex_H)
    print("\n=== Salt states (series: SH → Eva → Eco) ===")
    _row("Salt → SH:su_H",       superheater.su_H)
    _row("Salt → SH:ex_H / Eva:su_H",  superheater.ex_H)
    _row("Salt → Eva:ex_H / Eco:su_H", eva.ex_H)
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
# CONFIG 3 — REHEAT  (parallel salt split + outer boiler-branch convergence)
# Steam: Pump → Eco → Eva → SH → Turb_HP → RH → Turb_LP → Condenser → Pump
#
# Salt routing:
#   (1 - f_rh) * m_dot_salt  →  SH → Eva → Eco  (series, outer convergence loop)
#         f_rh * m_dot_salt  →  RH only           (fresh salt at T_salt_su, fixed)
#
# 2D sweep on f_rh × P_reheat run from __main__.
# =============================================================================

def _build_reheat_circuit(eta_pump, eta_turb,
                          eta_eco, eta_eva, eta_sh, eta_rh,
                          pinch_hx, pinch_cond,
                          HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh, CSource,
                          P_low, P_high, P_reheat, m_dot_st, n_disc):
    """Build and return a fresh RecursiveCircuit for the reheat split cycle."""
    T_sat_hi  = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo  = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_salt_su = HSource_sh.T

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

    cycle.add_source("SaltSource_SH",  HSource_sh,  cycle.components["Superheater"], "m-su_H")
    cycle.add_source("SaltSource_EVA", Salt_eva_su, cycle.components["Evaporator"],  "m-su_H")
    cycle.add_source("SaltSource_ECO", Salt_eco_su, cycle.components["Economiser"],  "m-su_H")
    cycle.add_source("SaltSource_RH",  HSource_rh,  cycle.components["Reheater"],    "m-su_H")
    cycle.add_source("ColdSource",     CSource,     cycle.components["Condenser"],   "m-su_C")

    cycle.set_fixed_properties(target="Pump:ex",    p=P_high)
    cycle.set_fixed_properties(target="Turb_HP:ex", p=P_reheat)
    cycle.set_fixed_properties(target="Turb_LP:ex", p=P_low)

    T_HP_su_g = T_salt_su - pinch_hx
    T_LP_su_g = T_salt_su - pinch_hx   # fresh salt → RH: same T as SH inlet

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
    cycle.set_cycle_guess(target='Evaporator:ex_C',  p=P_high,   T=T_sat_hi,    x=1)
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


def rankine_reheat(eta_pump, eta_turb,
                   eta_eco, eta_eva, eta_sh, eta_rh,
                   pinch_hx, pinch_cond,
                   T_salt_su, P_salt, m_dot_salt_total,
                   f_rh, CSource,
                   P_low, P_high, P_reheat, m_dot_st,
                   n_disc=10, tol_outer=0.1, max_outer=20):
    """
    Reheat Rankine — parallel salt split with outer boiler-branch convergence.,

    Parameters
    ----------
    T_salt_su        : salt supply temperature [K]
    P_salt           : salt pressure [Pa]
    m_dot_salt_total : total salt mass flow rate [kg/s]
    f_rh             : fraction of salt routed to RH (fresh, parallel) [0–1]
    tol_outer        : convergence tolerance for intermediate salt T [K]
    max_outer        : maximum outer iterations
    """
    T_sat_hi     = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    m_dot_boiler = m_dot_salt_total * (1.0 - f_rh)
    m_dot_rh     = m_dot_salt_total * f_rh

    # RH branch: always fresh salt at T_salt_su — no outer loop needed
    HSource_rh = _salt_connector(T_salt_su, P_salt, m_dot_rh)

    # Boiler branch: initial guesses for intermediate salt temperatures
    T_eva_su = T_salt_su - 10.0   # just below SH inlet salt T
    T_eco_su = T_sat_hi + 5.0     # just above saturation

    for outer_iter in range(max_outer):
        HSource_sh  = _salt_connector(T_salt_su, P_salt, m_dot_boiler)
        Salt_eva_su = _salt_connector(T_eva_su,  P_salt, m_dot_boiler)
        Salt_eco_su = _salt_connector(T_eco_su,  P_salt, m_dot_boiler)

        cycle, pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser = \
            _build_reheat_circuit(
                eta_pump, eta_turb, eta_eco, eta_eva, eta_sh, eta_rh,
                pinch_hx, pinch_cond,
                HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh, CSource,
                P_low, P_high, P_reheat, m_dot_st, n_disc,
            )
        cycle.plot_flag = False
        cycle.solve()

        T_eva_su_new = sh.ex_H.T    # salt leaving SH → entering Eva
        T_eco_su_new = eva.ex_H.T   # salt leaving Eva → entering Eco

        T_eva_su = T_eva_su_new
        T_eco_su = T_eco_su_new


    return cycle, pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser


def compute_cycle_performance_rh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser):
    m_dot   = turb_hp.su.m_dot
    W_pump  = m_dot * (pump.ex.h    - pump.su.h)            / 1000
    W_HP    = m_dot * (turb_hp.su.h - turb_hp.ex.h)         / 1000
    W_LP    = m_dot * (turb_lp.su.h - turb_lp.ex.h)         / 1000
    W_turb  = W_HP + W_LP
    W_net   = W_turb - W_pump
    Q_eco   = m_dot * (eco.ex_C.h   - eco.su_C.h)           / 1000
    Q_eva   = m_dot * (eva.ex_C.h   - eva.su_C.h)           / 1000
    Q_sh    = m_dot * (sh.ex_C.h    - sh.su_C.h)            / 1000
    Q_rh    = m_dot * (rh.ex_C.h    - rh.su_C.h)            / 1000
    Q_boil  = Q_eco + Q_eva + Q_sh + Q_rh
    Q_cond  = m_dot * (condenser.su_H.h - condenser.ex_H.h) / 1000
    eta_th  = W_net / Q_boil
    h_lp = turb_lp.ex.h
    p_lp = turb_lp.ex.p
    h_sat_v = PropsSI('H', 'P', p_lp, 'Q', 1, 'Water')
    if h_lp > h_sat_v:
        x_LP_ex = 1.0 + (h_lp - h_sat_v) / h_sat_v  # indicateur superchauffe
    else:
        x_LP_ex = PropsSI('Q', 'H', h_lp, 'P', p_lp, 'Water')
    return {
        'W_pump': W_pump, 'W_HP': W_HP, 'W_LP': W_LP,
        'W_turb': W_turb, 'W_net': W_net,
        'Q_eco': Q_eco, 'Q_eva': Q_eva, 'Q_sh': Q_sh, 'Q_rh': Q_rh,
        'Q_boiler': Q_boil, 'Q_cond': Q_cond,
        'eta_th': eta_th, 'x_LP_ex': x_LP_ex,
    }


def print_states_rh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser, f_rh):
    def _row(label, conn):
        print(f"  {label:<42} T={conn.T-273.15:7.2f} °C   p={conn.p/1e5:7.3f} bar   h={conn.h:10.2f} J/kg")
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
    print(f"\n=== Salt states  (split: f_rh={f_rh:.2f}) ===")
    print(f"  Boiler branch  ({(1-f_rh)*100:.0f}% of salt, series SH→Eva→Eco)")
    _row("  Salt → SH:su_H",           sh.su_H)
    _row("  Salt → SH:ex_H / Eva:su_H",sh.ex_H)
    _row("  Salt → Eva:ex_H / Eco:su_H",eva.ex_H)
    _row("  Salt → Eco:ex_H",          eco.ex_H)
    print(f"  RH branch  ({f_rh*100:.0f}% of salt, fresh at T_salt_su)")
    _row("  Salt → RH:su_H",           rh.su_H)
    _row("  Salt → RH:ex_H",           rh.ex_H)


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
    print(f"  x_LP_ex     : {perf['x_LP_ex']:.4f}  [-]  (steam quality at LP exit)")
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

    study_case   = "Reheat"    # "Simple" | "Superheated" | "Reheat"
    PRINT        = True
    DETAIL       = True
    PRINT_SCALE  = False
    W_net_target = 100         # [MW]

    # --- Component efficiencies ---
    eta_pump = 0.75
    eta_turb = 0.85

    # --- HX parameters ---
    eta_hx     = 0.95
    pinch_hx   = 5.0    # [K]
    pinch_cond = 5.0    # [K]

    # --- Cycle pressures ---
    P_high = 160e5     # Pa
    P_low  = 0.10e5    # Pa

    # --- Steam mass flow rate (specific basis) ---
    m_dot_st = 1.0     # kg/s

    # --- Hot source: Solar Salt ---
    T_salt_su      = 565 + 273.15   # K
    P_salt         = 2e5            # Pa
    m_dot_salt_tot = 100*5.99           # kg/s

    HSource = SolarSaltConnector()
    HSource.set_properties(T=T_salt_su, p=P_salt, m_dot=m_dot_salt_tot)
    print(f"Hot source : Solar Salt at {T_salt_su - 273.15:.1f} °C")

    # --- Cold source: cooling water ---
    T_cw_su = 20 + 273.15
    P_cw    = 3e5
    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=T_cw_su, P=P_cw, m_dot=100.0)

    # -------------------------------------------------------------------------
    if study_case == "Simple":
        cycle, pump, boiler, turbine, condenser = \
            rankine_simple(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                           HSource, CSource, P_low, P_high, m_dot_st)
        cycle.plot_flag = False
        cycle.solve()
        perf = compute_cycle_performance(pump, boiler, turbine, condenser)
        if PRINT:
            if DETAIL: print_states(pump, boiler, turbine, condenser)
            print_efficiency(perf)
        scale = scale_to_power(W_net_target, perf, PRINT_SCALE)

    # -------------------------------------------------------------------------
    elif study_case == "Superheated":
        # rankine_superheated() calls solve() internally inside the outer loop
        cycle, pump, eco, eva, superheater, turbine, condenser = \
            rankine_superheated(eta_pump, eta_turb, eta_hx, eta_hx, eta_hx,
                                pinch_hx, pinch_cond, HSource, CSource,
                                P_low, P_high, m_dot_st)
        perf = compute_cycle_performance_sh(pump, eco, eva, superheater, turbine, condenser)
        if PRINT:
            if DETAIL: print_states_sh(pump, eco, eva, superheater, turbine, condenser)
            print_efficiency_sh(perf)
        scale = scale_to_power_sh(W_net_target, perf, PRINT_SCALE)

    # -------------------------------------------------------------------------
    elif study_case == "Reheat":
        # 2D sweep: f_rh × P_reheat
        # rankine_reheat() calls solve() internally inside the outer loop
        f_rh_values   = [0.10, 0.20, 0.30, 0.40, 0.50]
        P_reheat_bars = [10, 15, 20, 25, 30]

        results = []

        print(f"\n=== Reheat sweep  "
              f"(P_high={P_high/1e5:.0f} bar, T_salt={T_salt_su-273.15:.0f}°C) ===")
        print(f"  {'f_rh':>6} {'P_rh[bar]':>10} {'eta_th[%]':>11} {'W_net[kW]':>11} "
              f"{'x_LP_ex':>9} {'T_LP_su[°C]':>12} {'Q_rh[kW]':>10}")
        print("  " + "-" * 75)

        for P_rh_bar in P_reheat_bars:
            for f_rh in f_rh_values:
                try:
                    cycle, pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser = \
                        rankine_reheat(
                            eta_pump, eta_turb,
                            eta_hx, eta_hx, eta_hx, eta_hx,
                            pinch_hx, pinch_cond,
                            T_salt_su, P_salt, m_dot_salt_tot,
                            f_rh, CSource,
                            P_low, P_high, P_rh_bar * 1e5, m_dot_st,
                        )
                    perf = compute_cycle_performance_rh(pump, eco, eva, sh, rh,
                                                        turb_hp, turb_lp, condenser)
                    results.append({
                        'f_rh':     f_rh,
                        'P_rh_bar': P_rh_bar,
                        'eta_th':   perf['eta_th'] * 100,
                        'W_net':    perf['W_net'],
                        'x_LP_ex':  perf['x_LP_ex'],
                        'T_LP_su':  turb_lp.su.T - 273.15,
                        'Q_rh':     perf['Q_rh'],
                        'pump': pump, 'eco': eco, 'eva': eva, 'sh': sh, 'rh': rh,
                        'turb_hp': turb_hp, 'turb_lp': turb_lp,
                        'condenser': condenser, 'perf': perf,
                    })
                    print(f"  {f_rh:>6.2f} {P_rh_bar:>10.0f} "
                          f"{perf['eta_th']*100:>11.2f} {perf['W_net']:>11.2f} "
                          f"{perf['x_LP_ex']:>9.4f} {turb_lp.su.T-273.15:>12.2f} "
                          f"{perf['Q_rh']:>10.2f}")
                except Exception as e:
                    print(f"  f_rh={f_rh:.2f}, P_rh={P_rh_bar} bar → failed: {e}")

        # Global optimum
        best = max(results, key=lambda r: r['eta_th'])
        print(f"\n→ Global optimum: f_rh={best['f_rh']:.2f}, "
              f"P_reheat={best['P_rh_bar']} bar  →  "
              f"η_th={best['eta_th']:.2f} %  x_LP_ex={best['x_LP_ex']:.4f}")

        if PRINT and DETAIL:
            print_states_rh(best['pump'], best['eco'], best['eva'],
                            best['sh'], best['rh'],
                            best['turb_hp'], best['turb_lp'],
                            best['condenser'], best['f_rh'])
            print_efficiency_rh(best['perf'])
        if PRINT_SCALE:
            scale = scale_to_power_rh(W_net_target, best['perf'], PRINT_SCALE)

        # --- Plot ---
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        ax1, ax2  = axes
        colors    = plt.cm.viridis(np.linspace(0.15, 0.85, len(P_reheat_bars)))

        for i, P_rh_bar in enumerate(P_reheat_bars):
            subset   = [r for r in results if r['P_rh_bar'] == P_rh_bar]
            f_vals   = [r['f_rh']    for r in subset]
            eta_vals = [r['eta_th']  for r in subset]
            x_vals   = [r['x_LP_ex'] for r in subset]
            ax1.plot(f_vals, eta_vals, 'o-', color=colors[i], label=f'{P_rh_bar} bar')
            ax2.plot(f_vals, x_vals,   's-', color=colors[i], label=f'{P_rh_bar} bar')

        ax1.axvline(best['f_rh'], color='k', linestyle='--', alpha=0.4,
                    label=f"Optimum f_rh={best['f_rh']:.2f}")
        ax1.set_xlabel('f_rh (salt fraction to RH) [-]')
        ax1.set_ylabel('η_th [%]')
        ax1.set_title('Thermal efficiency vs salt split')
        ax1.legend(title='P_reheat', fontsize=8)
        ax1.grid(True, alpha=0.3)

        ax2.axhline(0.88, color='red',  linestyle='--', alpha=0.6, label='x = 0.88 limit')
        ax2.axhline(1.00, color='gray', linestyle=':',  alpha=0.4, label='x = 1 (sat.)')
        ax2.set_xlabel('f_rh (salt fraction to RH) [-]')
        ax2.set_ylabel('x_LP_ex [-]')
        ax2.set_title('LP turbine exit quality vs salt split')
        ax2.legend(title='P_reheat', fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.suptitle(f'Reheat — P_high={P_high/1e5:.0f} bar, '
                     f'T_salt={T_salt_su-273.15:.0f}°C, η_turb={eta_turb}',
                     fontsize=11)
        plt.tight_layout()
        plt.savefig('reheat_sweep.png', dpi=150)
        plt.show()