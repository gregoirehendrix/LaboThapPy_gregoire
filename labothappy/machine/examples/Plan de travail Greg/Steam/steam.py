# -*- coding: utf-8 -*-
"""
Created on Apr 01 2026

@author: gregoire.hendrix

Steam Rankine cycle — Solar Salt heat source

Available configurations (choose in __main__):
    study_case = "Simple"       # Pump → Boiler → Turbine → Condenser → Pump

Future configurations to add:
    study_case = "Recuperated"  # with feedwater heater / preheater
    study_case = "Reheat"       # with reheat stage
"""

#%%
import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.expander.expander_csteff  import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from labothappy.component.pump.pump_csteff          import PumpCstEff


#%%
def rankine_simple(eta_pump, eta_turb, eta_boiler, eta_cond,
                   HSource, CSource,
                   P_low, P_high, m_dot_st,
                   T_turb_su_guess=None, T_turb_ex_guess=None):
    """
    Basic steam Rankine cycle.
    Pump → Boiler (hot: Solar Salt) → Turbine_St → Condenser (cold: water) → Pump

    Parameters
    ----------
    eta_pump, eta_turb   : isentropic efficiencies [-]
    eta_boiler, eta_cond : HX effectivenesses [-]
    HSource              : SolarSaltConnector — hot source
    CSource              : MassConnector     — cold source (cooling water)
    P_low, P_high        : cycle pressures [Pa]
    m_dot_st             : steam mass flow rate [kg/s]
    T_turb_su_guess      : turbine inlet temperature guess [K]
    T_turb_ex_guess      : turbine exit temperature guess [K]
    """

    # Default guesses if not provided
    T_sat_hi = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')

    if T_turb_su_guess is None:
        T_turb_su_guess = T_sat_hi + 30    # slightly superheated
    if T_turb_ex_guess is None:
        T_turb_ex_guess = T_sat_lo + 5

    # --- Build circuit ---
    cycle     = RecursiveCircuit('Water')
    pump      = PumpCstEff()
    boiler    = HexCstEff()
    turbine   = ExpanderCstEff()
    condenser = HexCstEff()

    pump.set_parameters(eta_is=eta_pump)
    turbine.set_parameters(eta_is=eta_turb)
    boiler.set_parameters(eta=eta_boiler)
    condenser.set_parameters(eta=eta_cond)

    cycle.add_component(pump,      "Pump")
    cycle.add_component(boiler,    "Boiler")
    cycle.add_component(turbine,   "Turbine_St")
    cycle.add_component(condenser, "Condenser")

    # Pump → Boiler → Turbine → Condenser → Pump
    cycle.link_components("Pump",       "m-ex",   "Boiler",     "m-su_C")
    cycle.link_components("Boiler",     "m-ex_C", "Turbine_St", "m-su")
    cycle.link_components("Turbine_St", "m-ex",   "Condenser",  "m-su_H")
    cycle.link_components("Condenser",  "m-ex_H", "Pump",       "m-su")

    # --- Sources ---
    cycle.add_source("HotSource",  HSource, cycle.components["Boiler"],    "m-su_H")
    cycle.add_source("ColdSource", CSource, cycle.components["Condenser"], "m-su_C")

    # --- Fixed pressures ---
    cycle.set_fixed_properties(target="Pump:ex",       p=P_high)
    cycle.set_fixed_properties(target="Turbine_St:ex", p=P_low)

    # --- Analytical guesses ---
    h_sat_lo     = PropsSI('H', 'P', P_low,  'Q', 0,   'Water')
    s_pump_su    = PropsSI('S', 'P', P_low,  'Q', 0,   'Water')
    h_pump_ex_is = PropsSI('H', 'P', P_high, 'S', s_pump_su, 'Water')
    h_pump_ex    = h_sat_lo + (h_pump_ex_is - h_sat_lo) / eta_pump
    T_pump_ex    = PropsSI('T', 'H', h_pump_ex, 'P', P_high, 'Water')

    s_turb_su    = PropsSI('S', 'T', T_turb_su_guess, 'P', P_high, 'Water')
    h_turb_su    = PropsSI('H', 'T', T_turb_su_guess, 'P', P_high, 'Water')
    h_turb_ex_is = PropsSI('H', 'P', P_low, 'S', s_turb_su, 'Water')
    h_turb_ex    = h_turb_su - eta_turb * (h_turb_su - h_turb_ex_is)
    T_turb_ex    = PropsSI('T', 'H', h_turb_ex, 'P', P_low, 'Water')

    cycle.set_cycle_guess(target='Pump:su',       p=P_low,  SC=3,               m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',       p=P_high, T=T_pump_ex)
    cycle.set_cycle_guess(target='Turbine_St:su', p=P_high, T=T_turb_su_guess,  m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turbine_St:ex', p=P_low,  T=T_turb_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',p=P_low,  T=T_sat_lo - 2)

    # --- Residual variables ---
    cycle.set_residual_variable(target='Turbine_St:ex', variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Turbine_St:ex', variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Pump:ex',       variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Pump:ex',       variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Boiler:ex_C',   variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Boiler:ex_C',   variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Condenser:ex_H',variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Condenser:ex_H',variable='p', tolerance=1e-6)

    return cycle, pump, boiler, turbine, condenser


#%%
def compute_cycle_performance(pump, boiler, turbine, condenser):
    m_dot = turbine.su.m_dot

    W_pump  = m_dot * (pump.ex.h   - pump.su.h)   / 1000   # kW
    W_turb  = m_dot * (turbine.su.h - turbine.ex.h) / 1000  # kW
    W_net   = W_turb - W_pump                               # kW
    Q_boil  = m_dot * (boiler.ex_C.h   - boiler.su_C.h)   / 1000   # kW
    Q_cond  = m_dot * (condenser.su_H.h - condenser.ex_H.h) / 1000  # kW
    eta_th  = W_net / Q_boil

    # Steam quality at turbine exit
    x_turb_ex = PropsSI('Q', 'H', turbine.ex.h, 'P', turbine.ex.p, 'Water')

    return {
        'W_pump':     W_pump,
        'W_turb':     W_turb,
        'W_net':      W_net,
        'Q_boiler':   Q_boil,
        'Q_cond':     Q_cond,
        'eta_th':     eta_th,
        'x_turb_ex':  x_turb_ex,
    }


#%%
def print_states(pump, boiler, turbine, condenser):
    print("=== Pump:su ===")
    print(f"  T     = {pump.su.T - 273.15:.2f} °C")
    print(f"  p     = {pump.su.p / 1e5:.3f} bar")
    print(f"  h     = {pump.su.h:.2f} J/kg")
    print(f"  m_dot = {pump.su.m_dot:.4f} kg/s")
    print("=== Pump:ex ===")
    print(f"  T     = {pump.ex.T - 273.15:.2f} °C")
    print(f"  p     = {pump.ex.p / 1e5:.2f} bar")
    print(f"  h     = {pump.ex.h:.2f} J/kg")
    print("=== Boiler:ex_C (turbine inlet) ===")
    print(f"  T     = {boiler.ex_C.T - 273.15:.2f} °C")
    print(f"  p     = {boiler.ex_C.p / 1e5:.2f} bar")
    print(f"  h     = {boiler.ex_C.h:.2f} J/kg")
    print("=== Turbine_St:ex ===")
    print(f"  T     = {turbine.ex.T - 273.15:.2f} °C")
    print(f"  p     = {turbine.ex.p / 1e5:.3f} bar")
    print(f"  h     = {turbine.ex.h:.2f} J/kg")
    print("=== Condenser:ex_H (pump inlet) ===")
    print(f"  T     = {condenser.ex_H.T - 273.15:.2f} °C")
    print(f"  p     = {condenser.ex_H.p / 1e5:.3f} bar")
    print(f"  h     = {condenser.ex_H.h:.2f} J/kg")


#%%
def print_efficiency(perf):
    print("=== Rankine Cycle Efficiency ===")
    print(f"  W_pump      : {perf['W_pump']:.2f} kW/kg/s")
    print(f"  W_turb      : {perf['W_turb']:.2f} kW/kg/s")
    print(f"  W_net       : {perf['W_net']:.2f} kW/kg/s")
    print(f"  Q_boiler    : {perf['Q_boiler']:.2f} kW/kg/s")
    print(f"  Q_cond      : {perf['Q_cond']:.2f} kW/kg/s")
    print(f"  eta_th      : {perf['eta_th'] * 100:.2f} %")
    print(f"  x_turb_ex   : {perf['x_turb_ex']:.4f}  [-]  (steam quality at turbine exit)")
    print("  --- First law check ---")
    balance = perf['Q_boiler'] - perf['Q_cond']
    print(f"  Q_boil - Q_cond : {balance:.4f} kW")
    print(f"  W_net           : {perf['W_net']:.4f} kW")
    print(f"  Discrepancy     : {abs(perf['W_net'] - balance):.6f} kW")
    print("================================")


#%%
def scale_to_power(W_net_target_MW, perf, PRINT_SCALE=True):
    """Scale all quantities to a target net electrical power [MW]."""
    scale           = (W_net_target_MW * 1000) / perf['W_net']
    W_net_scaled    = perf['W_net']    * scale / 1000
    Q_boiler_scaled = perf['Q_boiler'] * scale / 1000
    Q_cond_scaled   = perf['Q_cond']   * scale / 1000

    if PRINT_SCALE:
        print("=== Scaling to target net power ===")
        print(f"  W_net_target         : {W_net_target_MW:.2f} MW")
        print(f"  W_net_ref (1 kg/s)   : {perf['W_net']:.2f} kW")
        print(f"  Scale factor         : {scale:.4f}")
        print(f"  m_dot_steam required : {scale:.4f} kg/s")
        print(f"  W_pump               : {perf['W_pump'] * scale / 1000:.2f} MW")
        print(f"  W_turb               : {perf['W_turb'] * scale / 1000:.2f} MW")
        print(f"  W_net                : {W_net_scaled:.2f} MW")
        print(f"  Q_boiler             : {Q_boiler_scaled:.2f} MW")
        print(f"  Q_cond               : {Q_cond_scaled:.2f} MW")
        print(f"  eta_th               : {perf['eta_th'] * 100:.2f} %")
        print("====================================")
    return scale


#%%
if __name__ == "__main__":

    study_case  = "Simple"    # "Simple"  (more to come)
    PRINT       = True
    DETAIL      = True
    PRINT_SCALE = True
    W_net_target = 100          # [MW]

    # --- Component efficiencies ---
    eta_pump  = 0.75
    eta_turb  = 0.85
    eta_boil  = 0.95
    eta_cond  = 0.95

    # --- Cycle pressures ---
    P_high    = 160e5          # Pa
    P_low     = 0.10e5        # Pa

    # --- Steam mass flow rate (specific basis for scaling) ---
    m_dot_st  = 1.0           # kg/s

    # --- Hot source: Solar Salt ---
    T_salt_su = 565 + 273.15  # K
    P_salt    = 2e5           # Pa

    HSource = SolarSaltConnector()
    HSource.set_properties(T=T_salt_su, p=P_salt, m_dot=50.0)
    print(f"Hot source : Solar Salt at {T_salt_su - 273.15:.1f} °C")

    # --- Cold source: cooling water ---
    T_cw_su = 20 + 273.15     # K
    P_cw    = 3e5             # Pa

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=T_cw_su, P=P_cw, m_dot=100.0)

    # --- Turbine guesses ---
    T_turb_su_guess = T_salt_su - 10
    T_turb_ex_guess = PropsSI('T', 'P', P_low, 'Q', 1, 'Water') + 5

    # --- Run ---
    if study_case == "Simple":
        cycle, pump, boiler, turbine, condenser = \
            rankine_simple(eta_pump, eta_turb, eta_boil, eta_cond,
                           HSource, CSource,
                           P_low, P_high, m_dot_st,
                           T_turb_su_guess=T_turb_su_guess,
                           T_turb_ex_guess=T_turb_ex_guess)
        cycle.solve()
        perf = compute_cycle_performance(pump, boiler, turbine, condenser)
        if PRINT:
            if DETAIL: print_states(pump, boiler, turbine, condenser)
            print_efficiency(perf)
        scale = scale_to_power(W_net_target, perf, PRINT_SCALE)