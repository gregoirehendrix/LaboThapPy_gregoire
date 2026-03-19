# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:09:29 2026

@author: gregoire.hendrix
"""

from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff 
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff


def sCO2_basic(eta_cp, eta_tb, eta_heater, eta_cooler,
               HSource, CSource, P_low, P_high, T_c_su, m_dot_CO2):
    cycle = RecursiveCircuit('CO2')
    
    #%% Create components
    compressor = CompressorCstEff()
    expander   = ExpanderCstEff()
    gas_heater = HexCstEff()
    gas_cooler = HexCstEff()
    
    #%% Set parameters
    compressor.set_parameters(eta_is=eta_cp)
    expander.set_parameters(eta_is=eta_tb)
    gas_heater.set_parameters(eta=eta_heater)
    gas_cooler.set_parameters(eta=eta_cooler)
    
    #%% Add components to cycle
    cycle.add_component(compressor, "Compressor")
    cycle.add_component(gas_heater, "Gas_Heater")
    cycle.add_component(expander,   "Expander")
    cycle.add_component(gas_cooler, "Gas_Cooler")
    
    #%% Link components
    cycle.link_components("Compressor", "m-ex",   "Gas_Heater", "m-su_C")
    cycle.link_components("Gas_Heater", "m-ex_C", "Expander",   "m-su")
    cycle.link_components("Expander",   "m-ex",   "Gas_Cooler", "m-su_H")
    cycle.link_components("Gas_Cooler", "m-ex_H", "Compressor", "m-su")
    
    #%% External sources
    cycle.add_source("Hot_air",  HSource, cycle.components["Gas_Heater"], "m-su_H")
    cycle.add_source("Cold_air", CSource, cycle.components["Gas_Cooler"], "m-su_C")
    
    #%% Initial guesses
    cycle.set_cycle_guess(target='Compressor:su', m_dot=m_dot_CO2, T=T_c_su,       p=P_low)
    cycle.set_cycle_guess(target='Compressor:ex',                  T=T_c_su + 50,  p=P_high)
    cycle.set_cycle_guess(target='Expander:su',                    T=500 + 273.15, p=P_high)
    cycle.set_cycle_guess(target='Expander:ex',                    T=200 + 273.15, p=P_low)
    
    #%% Residuals
    cycle.set_residual_variable(target='Expander:ex',     variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Expander:ex',     variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Compressor:ex',   variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Compressor:ex',   variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Heater:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Heater:ex_C', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H', variable='p', tolerance=1e-3)

    return cycle, compressor, expander, gas_heater, gas_cooler


def sCO2_recuperated(eta_cp, eta_tb, eta_heater, eta_cooler, eta_recuperator,
                     HSource, CSource, P_low, P_high, T_c_su, m_dot_CO2):
    cycle = RecursiveCircuit('CO2')
    
    #%% Create components
    compressor  = CompressorCstEff()
    expander    = ExpanderCstEff()
    gas_heater  = HexCstEff()
    gas_cooler  = HexCstEff()
    recuperator = HexCstEff()
    
    #%% Set parameters
    compressor.set_parameters(eta_is=eta_cp)
    expander.set_parameters(eta_is=eta_tb)
    gas_heater.set_parameters(eta=eta_heater)
    gas_cooler.set_parameters(eta=eta_cooler)
    recuperator.set_parameters(eta=eta_recuperator)
    
    #%% Add components to cycle
    cycle.add_component(compressor,  "Compressor")
    cycle.add_component(recuperator, "Recuperator")
    cycle.add_component(gas_heater,  "Gas_Heater")
    cycle.add_component(expander,    "Expander")
    cycle.add_component(gas_cooler,  "Gas_Cooler")
    
    #%% Link components
    # High pressure side: Compressor -> Recuperator (cold) -> Gas_Heater -> Expander
    cycle.link_components("Compressor",  "m-ex",   "Recuperator", "m-su_C")
    cycle.link_components("Recuperator", "m-ex_C", "Gas_Heater",  "m-su_C")
    cycle.link_components("Gas_Heater",  "m-ex_C", "Expander",    "m-su")
    # Low pressure side: Expander -> Recuperator (hot) -> Gas_Cooler -> Compressor
    cycle.link_components("Expander",    "m-ex",   "Recuperator", "m-su_H")
    cycle.link_components("Recuperator", "m-ex_H", "Gas_Cooler",  "m-su_H")
    cycle.link_components("Gas_Cooler",  "m-ex_H", "Compressor",  "m-su")
    
    #%% External sources
    cycle.add_source("Hot_air",  HSource, cycle.components["Gas_Heater"], "m-su_H")
    cycle.add_source("Cold_air", CSource, cycle.components["Gas_Cooler"], "m-su_C")
    
    #%% Initial guesses
    cycle.set_cycle_guess(target='Compressor:su',    m_dot=m_dot_CO2, T=T_c_su,        p=P_low)
    cycle.set_cycle_guess(target='Compressor:ex',                     T=T_c_su + 50,   p=P_high)
    cycle.set_cycle_guess(target='Expander:su',                       T=500 + 273.15,  p=P_high)
    cycle.set_cycle_guess(target='Expander:ex',                       T=200 + 273.15,  p=P_low)
    # Recuperator inlet guesses — required to unblock the recursive solver
    cycle.set_cycle_guess(target='Recuperator:su_C', m_dot=m_dot_CO2, T=T_c_su + 50,  p=P_high)
    cycle.set_cycle_guess(target='Recuperator:su_H', m_dot=m_dot_CO2, T=200 + 273.15, p=P_low)
    # Recuperator outlet guesses
    cycle.set_cycle_guess(target='Recuperator:ex_C',                  T=180 + 273.15, p=P_high)
    cycle.set_cycle_guess(target='Recuperator:ex_H',                  T=60  + 273.15, p=P_low)
    
    #%% Residuals
    cycle.set_residual_variable(target='Expander:ex',      variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Expander:ex',      variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Compressor:ex',    variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Compressor:ex',    variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Heater:ex_C',  variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Heater:ex_C',  variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H',  variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H',  variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_H', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_H', variable='p', tolerance=1e-3)

    return cycle, compressor, expander, gas_heater, gas_cooler, recuperator


def compute_cycle_performance(compressor, expander, gas_heater, gas_cooler, recuperator=None):
    """
    Compute cycle performance using connector enthalpies (reliable values).
    Returns a dict with all key performance indicators in kW.
    """
    m_dot = compressor.su.m_dot  # kg/s

    # Work computed from enthalpies — more reliable than W.W_dot
    W_comp = m_dot * (compressor.ex.h - compressor.su.h) / 1000   # kW
    W_exp  = m_dot * (expander.su.h   - expander.ex.h)   / 1000   # kW
    W_net  = W_exp - W_comp                                         # kW

    # Heat exchanged computed from enthalpies
    Q_heater = m_dot * (gas_heater.ex_C.h - gas_heater.su_C.h) / 1000   # kW
    Q_cooler = m_dot * (gas_cooler.su_H.h - gas_cooler.ex_H.h) / 1000   # kW

    eta = W_net / Q_heater

    result = {
        'W_comp':   W_comp,
        'W_exp':    W_exp,
        'W_net':    W_net,
        'Q_heater': Q_heater,
        'Q_cooler': Q_cooler,
        'eta':      eta,
    }

    if recuperator is not None:
        Q_recup = m_dot * (recuperator.ex_C.h - recuperator.su_C.h) / 1000   # kW
        result['Q_recup'] = Q_recup

    return result


def print_states(compressor, expander, gas_cooler, recuperator=None):
    """Print the thermodynamic states at each key node of the cycle."""
    print("=== Compressor:su ===")
    print(f"  T     = {compressor.su.T - 273.15:.2f} °C")
    print(f"  p     = {compressor.su.p / 1e5:.2f} bar")
    print(f"  h     = {compressor.su.h:.2f} J/kg")
    print(f"  m_dot = {compressor.su.m_dot:.2f} kg/s")
    print("=== Compressor:ex ===")
    print(f"  T     = {compressor.ex.T - 273.15:.2f} °C")
    print(f"  p     = {compressor.ex.p / 1e5:.2f} bar")
    print("=== Expander:su ===")
    print(f"  T     = {expander.su.T - 273.15:.2f} °C")
    print(f"  p     = {expander.su.p / 1e5:.2f} bar")
    print("=== Expander:ex ===")
    print(f"  T     = {expander.ex.T - 273.15:.2f} °C")
    print(f"  p     = {expander.ex.p / 1e5:.2f} bar")
    print("=== Gas_Cooler:su_H ===")
    print(f"  T     = {gas_cooler.su_H.T - 273.15:.2f} °C")
    print(f"  p     = {gas_cooler.su_H.p / 1e5:.2f} bar")
    print("=== Gas_Cooler:ex_H ===")
    print(f"  T     = {gas_cooler.ex_H.T - 273.15:.2f} °C")
    print(f"  p     = {gas_cooler.ex_H.p / 1e5:.2f} bar")
    print(f"  h     = {gas_cooler.ex_H.h:.2f} J/kg")
    print("=== Gas_Cooler:su_C ===")
    print(f"  T     = {gas_cooler.su_C.T - 273.15:.2f} °C")
    print(f"  p     = {gas_cooler.su_C.p / 1e5:.2f} bar")
    if recuperator is not None:
        print("=== Recuperator:su_H ===")
        print(f"  T     = {recuperator.su_H.T - 273.15:.2f} °C")
        print(f"  p     = {recuperator.su_H.p / 1e5:.2f} bar")
        print("=== Recuperator:ex_H ===")
        print(f"  T     = {recuperator.ex_H.T - 273.15:.2f} °C")
        print(f"  p     = {recuperator.ex_H.p / 1e5:.2f} bar")
        print("=== Recuperator:su_C ===")
        print(f"  T     = {recuperator.su_C.T - 273.15:.2f} °C")
        print(f"  p     = {recuperator.su_C.p / 1e5:.2f} bar")
        print("=== Recuperator:ex_C ===")
        print(f"  T     = {recuperator.ex_C.T - 273.15:.2f} °C")
        print(f"  p     = {recuperator.ex_C.p / 1e5:.2f} bar")


def print_efficiency(perf):
    """Print the cycle energy balance and thermal efficiency from a performance dict."""
    print("=== Cycle Efficiency ===")
    print(f"  W_comp   : {perf['W_comp']:.2f} kW")
    print(f"  W_exp    : {perf['W_exp']:.2f} kW")
    print(f"  W_net    : {perf['W_net']:.2f} kW")
    print(f"  Q_heater : {perf['Q_heater']:.2f} kW")
    print(f"  Q_cooler : {perf['Q_cooler']:.2f} kW")
    if 'Q_recup' in perf:
        print(f"  Q_recup  : {perf['Q_recup']:.2f} kW")
    print(f"  eta      : {perf['eta'] * 100:.2f} %")
    # First law check
    balance = perf['Q_heater'] - perf['Q_cooler']
    print(f"  --- First law check ---")
    print(f"  Q_heater - Q_cooler : {balance:.2f} kW")
    print(f"  W_net               : {perf['W_net']:.2f} kW")
    print(f"  Discrepancy         : {abs(perf['W_net'] - balance):.2f} kW")
    print("========================")


def scale_to_power(W_net_target_MW, perf):
    """
    Scale the mass flow rate to reach the target net power output.
    The cycle must have been solved with m_dot_CO2 = 1 kg/s as reference.
    All heat and power values are linearly proportional to mass flow rate.
    Uses enthalpie-based performance dict for reliable values.
    """
    W_net_target_kW = W_net_target_MW * 1000

    # Scale factor to reach target net power
    scale = W_net_target_kW / perf['W_net']

    print("=== Scaling to target net power ===")
    print(f"  W_net_target       : {W_net_target_MW:.2f} MW")
    print(f"  W_net_ref (1 kg/s) : {perf['W_net']:.2f} kW")
    print(f"  Scale factor       : {scale:.4f}")
    print(f"  m_dot_CO2 required : {scale:.4f} kg/s")
    print(f"  W_comp             : {perf['W_comp']   * scale / 1000:.2f} MW")
    print(f"  W_exp              : {perf['W_exp']    * scale / 1000:.2f} MW")
    print(f"  W_net              : {perf['W_net']    * scale / 1000:.2f} MW")
    print(f"  Q_heater           : {perf['Q_heater'] * scale / 1000:.2f} MW")
    print(f"  Q_cooler           : {perf['Q_cooler'] * scale / 1000:.2f} MW")
    if 'Q_recup' in perf:
        print(f"  Q_recup            : {perf['Q_recup'] * scale / 1000:.2f} MW")
    print(f"  eta                : {perf['eta'] * 100:.2f} %")
    # First law check at scaled power
    balance_scaled = (perf['Q_heater'] - perf['Q_cooler']) * scale / 1000
    W_net_scaled   = perf['W_net'] * scale / 1000
    print(f"  --- First law check ---")
    print(f"  Q_heater - Q_cooler : {balance_scaled:.2f} MW")
    print(f"  W_net               : {W_net_scaled:.2f} MW")
    print(f"  Discrepancy         : {abs(W_net_scaled - balance_scaled):.4f} MW")
    print("====================================")

    return scale


#%%
if __name__ == "__main__":

    study_case   = "Recuperated"   # "Simple" or "Recuperated"
    PRINT        = True
    W_net_target = 100             # [MW] — target net power output

    #%% Common inputs
    P_low_guess  = 80  * 1e5      # [Pa]
    P_high       = 240 * 1e5      # [Pa]
    T_c_su_guess = 40  + 273.15   # [K]
    m_dot_CO2    = 1              # [kg/s] — reference mass flow rate for scaling

    eta_heater   = 0.968382
    eta_cooler   = 0.9
    eta_tb       = 0.9
    eta_cp       = 0.8

    # Hot source (hot air / flue gas)
    T_hot_su  = 565 + 273.15      # [K]
    p_hot_su  = 1.2 * 1e5         # [Pa]
    m_dot_hot = 20                 # [kg/s]

    # Cold source (ambient air)
    T_cold_su  = 25 + 273.15      # [K]
    p_cold_su  = 101325           # [Pa]
    m_dot_cold = 15               # [kg/s]

    HSource = MassConnector()
    HSource.set_properties(fluid='air', T=T_hot_su,  p=p_hot_su,  m_dot=m_dot_hot)

    CSource = MassConnector()
    CSource.set_properties(fluid='air', T=T_cold_su, p=p_cold_su, m_dot=m_dot_cold)

    #%% Solve and scale
    if study_case == "Simple":
        cycle, compressor, expander, gas_heater, gas_cooler = \
            sCO2_basic(eta_cp, eta_tb, eta_heater, eta_cooler,
                       HSource, CSource, P_low_guess, P_high, T_c_su_guess, m_dot_CO2)
        cycle.solve()
        cycle.plot_cycle_Ts()
        perf = compute_cycle_performance(compressor, expander, gas_heater, gas_cooler)
        if PRINT:
            print_states(compressor, expander, gas_cooler)
            print_efficiency(perf)
        scale_to_power(W_net_target, perf)

    elif study_case == "Recuperated":
        eta_recuperator = 0.85
        cycle, compressor, expander, gas_heater, gas_cooler, recuperator = \
            sCO2_recuperated(eta_cp, eta_tb, eta_heater, eta_cooler, eta_recuperator,
                             HSource, CSource, P_low_guess, P_high, T_c_su_guess, m_dot_CO2)
        cycle.solve()
        cycle.plot_cycle_Ts()
        perf = compute_cycle_performance(compressor, expander, gas_heater, gas_cooler, recuperator)
        if PRINT:
            print_states(compressor, expander, gas_cooler, recuperator)
            print_efficiency(perf)
        scale_to_power(W_net_target, perf)