# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:09:29 2026
@author: gregoire.hendrix
"""

#%%
import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc
from labothappy.component.tank.tank_spliter import TankSpliter
from labothappy.component.tank.tank_mixer import TankMixer

#%%
def compute_salt_mdot(Q_heater_kW, T_salt_in_K, T_salt_out_K):
    """Solar salt mass flow rate [kg/s] from energy balance (Zavoico 2001)."""
    T_in_C  = T_salt_in_K  - 273.15
    T_out_C = T_salt_out_K - 273.15
    delta_h = 1443.0 * (T_in_C - T_out_C) + 0.086 * (T_in_C**2 - T_out_C**2)
    return (Q_heater_kW * 1000) / delta_h


#%%
def sCO2_basic(eta_cp, eta_tb, eta_heater, eta_cooler,
               HSource, CSource, P_low, P_high, T_c_su, m_dot_CO2,
               pinch_heater=5, n_disc=20,
               T_exp_su_guess=None, T_exp_ex_guess=None):
    if T_exp_su_guess is None: T_exp_su_guess = 598.8 + 273.15
    if T_exp_ex_guess is None: T_exp_ex_guess = 488.4 + 273.15

    cycle = RecursiveCircuit('CO2')

    compressor = CompressorCstEff()
    expander   = ExpanderCstEff()
    gas_heater = HexCstEffDisc()
    gas_cooler = HexCstEff()

    compressor.set_parameters(eta_is=eta_cp)
    expander.set_parameters(eta_is=eta_tb)
    gas_heater.set_parameters(**{'eta_max': eta_heater, 'n_disc': n_disc, 'Pinch_min': pinch_heater})
    gas_cooler.set_parameters(eta=eta_cooler)

    cycle.add_component(compressor, "Compressor")
    cycle.add_component(gas_heater, "Gas_Heater")
    cycle.add_component(expander,   "Expander")
    cycle.add_component(gas_cooler, "Gas_Cooler")

    cycle.link_components("Compressor", "m-ex",   "Gas_Heater", "m-su_C")
    cycle.link_components("Gas_Heater", "m-ex_C", "Expander",   "m-su")
    cycle.link_components("Expander",   "m-ex",   "Gas_Cooler", "m-su_H")
    cycle.link_components("Gas_Cooler", "m-ex_H", "Compressor", "m-su")

    cycle.add_source("Hot_source",  HSource, cycle.components["Gas_Heater"], "m-su_H")
    cycle.add_source("Cold_source", CSource, cycle.components["Gas_Cooler"], "m-su_C")

    cycle.set_cycle_guess(target='Compressor:su', m_dot=m_dot_CO2, T=T_c_su,          p=P_low)
    cycle.set_cycle_guess(target='Compressor:ex',                  T=T_c_su + 50,     p=P_high)
    cycle.set_cycle_guess(target='Expander:su',   m_dot=m_dot_CO2, T=T_exp_su_guess,  p=P_high)
    cycle.set_cycle_guess(target='Expander:ex',                    T=T_exp_ex_guess,  p=P_low)

    cycle.set_residual_variable(target='Expander:ex',     variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Expander:ex',     variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Compressor:ex',   variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Compressor:ex',   variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Heater:ex_C', variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Heater:ex_C', variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H', variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H', variable='p', tolerance=1e-6)

    return cycle, compressor, expander, gas_heater, gas_cooler

#%%
def sCO2_recuperated(eta_cp, eta_tb, eta_heater, eta_cooler, eta_recuperator,
                     HSource, CSource, P_low, P_high, T_c_su, m_dot_CO2,
                     pinch_heater=5, pinch_recup=2, n_disc=20,
                     T_exp_su_guess=None, T_exp_ex_guess=None):
    if T_exp_su_guess is None: T_exp_su_guess = 598.8 + 273.15
    if T_exp_ex_guess is None: T_exp_ex_guess = 488.4 + 273.15

    cycle = RecursiveCircuit('CO2')

    compressor  = CompressorCstEff()
    expander    = ExpanderCstEff()
    gas_heater  = HexCstEffDisc()
    gas_cooler  = HexCstEff()
    recuperator = HexCstEffDisc()

    compressor.set_parameters(eta_is=eta_cp)
    expander.set_parameters(eta_is=eta_tb)
    gas_heater.set_parameters(**{'eta_max': eta_heater,      'n_disc': n_disc, 'Pinch_min': pinch_heater})
    gas_cooler.set_parameters(eta=eta_cooler)
    recuperator.set_parameters(**{'eta_max': eta_recuperator, 'n_disc': n_disc, 'Pinch_min': pinch_recup})

    cycle.add_component(compressor,  "Compressor")
    cycle.add_component(recuperator, "Recuperator")
    cycle.add_component(gas_heater,  "Gas_Heater")
    cycle.add_component(expander,    "Expander")
    cycle.add_component(gas_cooler,  "Gas_Cooler")

    cycle.link_components("Compressor",  "m-ex",   "Recuperator", "m-su_C")
    cycle.link_components("Recuperator", "m-ex_C", "Gas_Heater",  "m-su_C")
    cycle.link_components("Gas_Heater",  "m-ex_C", "Expander",    "m-su")
    cycle.link_components("Expander",    "m-ex",   "Recuperator", "m-su_H")
    cycle.link_components("Recuperator", "m-ex_H", "Gas_Cooler",  "m-su_H")
    cycle.link_components("Gas_Cooler",  "m-ex_H", "Compressor",  "m-su")

    cycle.add_source("Hot_source",  HSource, cycle.components["Gas_Heater"], "m-su_H")
    cycle.add_source("Cold_source", CSource, cycle.components["Gas_Cooler"], "m-su_C")

    cycle.set_cycle_guess(target='Compressor:su',    m_dot=m_dot_CO2, T=T_c_su,          p=P_low)
    cycle.set_cycle_guess(target='Compressor:ex',                     T=T_c_su + 50,     p=P_high)
    cycle.set_cycle_guess(target='Expander:su',      m_dot=m_dot_CO2, T=T_exp_su_guess,  p=P_high)
    cycle.set_cycle_guess(target='Expander:ex',                       T=T_exp_ex_guess,  p=P_low)
    cycle.set_cycle_guess(target='Recuperator:su_C', m_dot=m_dot_CO2, T=T_c_su + 50,    p=P_high)
    cycle.set_cycle_guess(target='Recuperator:su_H', m_dot=m_dot_CO2, T=T_exp_ex_guess, p=P_low)
    cycle.set_cycle_guess(target='Recuperator:ex_C',                  T=180 + 273.15,   p=P_high)
    cycle.set_cycle_guess(target='Recuperator:ex_H',                  T=60  + 273.15,   p=P_low)

    cycle.set_residual_variable(target='Expander:ex',      variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Expander:ex',      variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Compressor:ex',    variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Compressor:ex',    variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Heater:ex_C',  variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Heater:ex_C',  variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H',  variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H',  variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Recuperator:ex_H', variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Recuperator:ex_H', variable='p', tolerance=1e-6)

    return cycle, compressor, expander, gas_heater, gas_cooler, recuperator

#%%
def sCO2_recompression(eta_cp, eta_rc, eta_tb, eta_heater, eta_cooler,
                       eta_recup_lt, eta_recup_ht, split_ratio,
                       HSource, CSource, P_low, P_high, T_c_su, m_dot_CO2,
                       pinch_heater=5, pinch_recup=2, n_disc=20,
                       T_exp_su_guess=None, T_exp_ex_guess=None):
    """
    Recompression sCO2 Brayton cycle.
    Gas_Heater / RecupLT / RecupHT : HexCstEffDisc (supports SolarSalt)
    Gas_Cooler                     : HexCstEff
    split_ratio    : fraction alpha directed to Gas_Cooler + main Compressor.
    T_exp_su_guess : expander inlet guess [K] — derived from hot source T.
    T_exp_ex_guess : expander outlet guess [K] — rough estimate.
    """
    if T_exp_su_guess is None: T_exp_su_guess = 598.8 + 273.15
    if T_exp_ex_guess is None: T_exp_ex_guess = 488.4 + 273.15

    cycle = RecursiveCircuit('CO2')
    rep_spliter = [split_ratio, 1 - split_ratio]

    compressor   = CompressorCstEff()
    recompressor = CompressorCstEff()
    expander     = ExpanderCstEff()
    gas_heater   = HexCstEffDisc()
    gas_cooler   = HexCstEff()
    recup_lt     = HexCstEffDisc()
    recup_ht     = HexCstEffDisc()
    spliter      = TankSpliter(outlet_repartition=rep_spliter)
    mixer        = TankMixer(n_inlets=2)

    compressor.set_parameters(eta_is=eta_cp)
    recompressor.set_parameters(eta_is=eta_rc)
    expander.set_parameters(eta_is=eta_tb)
    gas_heater.set_parameters(**{'eta_max': eta_heater,  'n_disc': n_disc, 'Pinch_min': pinch_heater})
    gas_cooler.set_parameters(eta=eta_cooler)
    recup_lt.set_parameters(**{'eta_max':  eta_recup_lt, 'n_disc': n_disc, 'Pinch_min': pinch_recup})
    recup_ht.set_parameters(**{'eta_max':  eta_recup_ht, 'n_disc': n_disc, 'Pinch_min': pinch_recup})

    cycle.add_component(compressor,   "Compressor")
    cycle.add_component(recompressor, "Recompressor")
    cycle.add_component(expander,     "Expander")
    cycle.add_component(gas_heater,   "Gas_Heater")
    cycle.add_component(gas_cooler,   "Gas_Cooler")
    cycle.add_component(recup_lt,     "RecupLT")
    cycle.add_component(recup_ht,     "RecupHT")
    cycle.add_component(spliter,      "Spliter")
    cycle.add_component(mixer,        "Mixer")

    # Low pressure side
    cycle.link_components("Expander",     "m-ex",   "RecupHT",      "m-su_H")
    cycle.link_components("RecupHT",      "m-ex_H", "RecupLT",      "m-su_H")
    cycle.link_components("RecupLT",      "m-ex_H", "Spliter",      "m-su")
    # Branch alpha -> Gas_Cooler -> Compressor -> Mixer
    cycle.link_components("Spliter",      "m-ex_1", "Gas_Cooler",   "m-su_H")
    cycle.link_components("Gas_Cooler",   "m-ex_H", "Compressor",   "m-su")
    cycle.link_components("Compressor",   "m-ex",   "Mixer",        "m-su_1")
    # Branch (1-alpha) -> Recompressor -> Mixer
    cycle.link_components("Spliter",      "m-ex_2", "Recompressor", "m-su")
    cycle.link_components("Recompressor", "m-ex",   "Mixer",        "m-su_2")
    # High pressure side
    cycle.link_components("Mixer",        "m-ex",   "RecupLT",      "m-su_C")
    cycle.link_components("RecupLT",      "m-ex_C", "RecupHT",      "m-su_C")
    cycle.link_components("RecupHT",      "m-ex_C", "Gas_Heater",   "m-su_C")
    cycle.link_components("Gas_Heater",   "m-ex_C", "Expander",     "m-su")

    cycle.add_source("Hot_source",  HSource, cycle.components["Gas_Heater"], "m-su_H")
    cycle.add_source("Cold_source", CSource, cycle.components["Gas_Cooler"], "m-su_C")

    # Intermediate temperature guesses scaled to actual operating conditions
    T_recupHT_ex_H = T_c_su + 165         # ~202°C at base conditions
    T_recupHT_ex_C = T_exp_su_guess - 160  # ~440°C at base conditions

    cycle.set_cycle_guess(target='Expander:su',      m_dot=m_dot_CO2,                   T=T_exp_su_guess,  p=P_high)
    cycle.set_cycle_guess(target='Expander:ex',                                          T=T_exp_ex_guess,  p=P_low)
    cycle.set_cycle_guess(target='Compressor:su',    m_dot=m_dot_CO2 * split_ratio,      T=T_c_su,          p=P_low)
    cycle.set_cycle_guess(target='Compressor:ex',                                        T=T_c_su + 50,     p=P_high)
    cycle.set_cycle_guess(target='Recompressor:su',  m_dot=m_dot_CO2 * (1-split_ratio),  T=T_c_su + 60,     p=P_low)
    cycle.set_cycle_guess(target='Recompressor:ex',                                      T=T_c_su + 110,    p=P_high)
    cycle.set_cycle_guess(target='RecupHT:su_H',     m_dot=m_dot_CO2,                   T=T_exp_ex_guess,  p=P_low)
    cycle.set_cycle_guess(target='RecupHT:ex_H',     m_dot=m_dot_CO2,                   T=T_recupHT_ex_H,  p=P_low)
    cycle.set_cycle_guess(target='RecupHT:su_C',     m_dot=m_dot_CO2,                   T=T_c_su + 100,    p=P_high)
    cycle.set_cycle_guess(target='RecupHT:ex_C',                                         T=T_recupHT_ex_C,  p=P_high)
    cycle.set_cycle_guess(target='RecupLT:su_H',     m_dot=m_dot_CO2,                   T=T_recupHT_ex_H,  p=P_low)
    cycle.set_cycle_guess(target='RecupLT:ex_H',                                         T=T_c_su + 50,     p=P_low)
    cycle.set_cycle_guess(target='RecupLT:su_C',     m_dot=m_dot_CO2,                   T=T_c_su + 50,     p=P_high)
    cycle.set_cycle_guess(target='RecupLT:ex_C',                                         T=T_c_su + 150,    p=P_high)

    cycle.set_residual_variable(target='Expander:ex',       variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Expander:ex',       variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Compressor:ex',     variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Compressor:ex',     variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Heater:ex_C',   variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Heater:ex_C',   variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H',   variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Gas_Cooler:ex_H',   variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='RecupHT:su_H',      variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='RecupHT:su_H',      variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='RecupHT:ex_H',      variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='RecupHT:ex_H',      variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='RecupLT:ex_H',      variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='RecupLT:ex_H',      variable='p', tolerance=1e-6)
    cycle.set_residual_variable(target='Spliter:ex_2',      variable='h', tolerance=1e-6)
    cycle.set_residual_variable(target='Spliter:ex_2',      variable='p', tolerance=1e-6)

    return cycle, compressor, recompressor, expander, gas_heater, gas_cooler, recup_lt, recup_ht, spliter, mixer


#%%
def compute_cycle_performance(compressor, expander, gas_heater, gas_cooler,
                               recompressor=None, recuperator=None,
                               recup_lt=None, recup_ht=None):
    m_dot_cp  = compressor.su.m_dot
    m_dot_tot = expander.su.m_dot

    W_comp   = m_dot_cp  * (compressor.ex.h - compressor.su.h) / 1000
    W_exp    = m_dot_tot * (expander.su.h   - expander.ex.h)   / 1000
    W_recomp = (recompressor.su.m_dot * (recompressor.ex.h - recompressor.su.h) / 1000
                if recompressor is not None else 0.0)

    W_net    = W_exp - W_comp - W_recomp
    Q_heater = m_dot_tot * (gas_heater.ex_C.h - gas_heater.su_C.h) / 1000
    Q_cooler = (Q_heater - W_net if recompressor is not None
                else m_dot_cp * (gas_cooler.su_H.h - gas_cooler.ex_H.h) / 1000)

    result = {
        'W_comp': W_comp, 'W_recomp': W_recomp, 'W_exp': W_exp,
        'W_net': W_net, 'Q_heater': Q_heater, 'Q_cooler': Q_cooler,
        'eta': W_net / Q_heater,
    }
    if recuperator is not None:
        result['Q_recup'] = m_dot_tot * (recuperator.ex_C.h - recuperator.su_C.h) / 1000
    if recup_lt is not None and recup_ht is not None:
        result['Q_recup_lt'] = m_dot_tot * (recup_lt.ex_C.h - recup_lt.su_C.h) / 1000
        result['Q_recup_ht'] = m_dot_tot * (recup_ht.ex_C.h - recup_ht.su_C.h) / 1000
    return result

#%%
def print_states(compressor, expander, gas_cooler,
                 recompressor=None, recuperator=None,
                 recup_lt=None, recup_ht=None):
    print("=== Compressor:su ===")
    print(f"  T     = {compressor.su.T - 273.15:.2f} °C")
    print(f"  p     = {compressor.su.p / 1e5:.2f} bar")
    print(f"  h     = {compressor.su.h:.2f} J/kg")
    print(f"  m_dot = {compressor.su.m_dot:.4f} kg/s")
    print("=== Compressor:ex ===")
    print(f"  T     = {compressor.ex.T - 273.15:.2f} °C")
    print(f"  p     = {compressor.ex.p / 1e5:.2f} bar")
    if recompressor is not None:
        print("=== Recompressor:su ===")
        print(f"  T     = {recompressor.su.T - 273.15:.2f} °C")
        print(f"  p     = {recompressor.su.p / 1e5:.2f} bar")
        print(f"  m_dot = {recompressor.su.m_dot:.4f} kg/s")
        print("=== Recompressor:ex ===")
        print(f"  T     = {recompressor.ex.T - 273.15:.2f} °C")
        print(f"  p     = {recompressor.ex.p / 1e5:.2f} bar")
    print("=== Expander:su ===")
    print(f"  T     = {expander.su.T - 273.15:.2f} °C")
    print(f"  p     = {expander.su.p / 1e5:.2f} bar")
    print("=== Expander:ex ===")
    print(f"  T     = {expander.ex.T - 273.15:.2f} °C")
    print(f"  p     = {expander.ex.p / 1e5:.2f} bar")
    print("=== Gas_Cooler:su_H / ex_H ===")
    print(f"  su_H T = {gas_cooler.su_H.T - 273.15:.2f} °C")
    print(f"  ex_H T = {gas_cooler.ex_H.T - 273.15:.2f} °C")
    if recuperator is not None:
        print("=== Recuperator:su_H / ex_H / ex_C ===")
        print(f"  su_H T = {recuperator.su_H.T - 273.15:.2f} °C")
        print(f"  ex_H T = {recuperator.ex_H.T - 273.15:.2f} °C")
        print(f"  ex_C T = {recuperator.ex_C.T - 273.15:.2f} °C")
    if recup_lt is not None and recup_ht is not None:
        print("=== RecupLT:su_H / ex_H / su_C / ex_C ===")
        print(f"  su_H T = {recup_lt.su_H.T - 273.15:.2f} °C")
        print(f"  ex_H T = {recup_lt.ex_H.T - 273.15:.2f} °C")
        print(f"  su_C T = {recup_lt.su_C.T - 273.15:.2f} °C")
        print(f"  ex_C T = {recup_lt.ex_C.T - 273.15:.2f} °C")
        print("=== RecupHT:su_H / ex_H / su_C / ex_C ===")
        print(f"  su_H T = {recup_ht.su_H.T - 273.15:.2f} °C")
        print(f"  ex_H T = {recup_ht.ex_H.T - 273.15:.2f} °C")
        print(f"  su_C T = {recup_ht.su_C.T - 273.15:.2f} °C")
        print(f"  ex_C T = {recup_ht.ex_C.T - 273.15:.2f} °C")

#%%
def print_efficiency(perf):
    print("=== Cycle Efficiency ===")
    print(f"  W_comp   : {perf['W_comp']:.2f} kW")
    if perf.get('W_recomp', 0) > 0:
        print(f"  W_recomp : {perf['W_recomp']:.2f} kW")
    print(f"  W_exp    : {perf['W_exp']:.2f} kW")
    print(f"  W_net    : {perf['W_net']:.2f} kW")
    print(f"  Q_heater : {perf['Q_heater']:.2f} kW")
    print(f"  Q_cooler : {perf['Q_cooler']:.2f} kW")
    if 'Q_recup'    in perf: print(f"  Q_recup    : {perf['Q_recup']:.2f} kW")
    if 'Q_recup_lt' in perf:
        print(f"  Q_recup_LT : {perf['Q_recup_lt']:.2f} kW")
        print(f"  Q_recup_HT : {perf['Q_recup_ht']:.2f} kW")
    print(f"  eta      : {perf['eta'] * 100:.2f} %")
    balance = perf['Q_heater'] - perf['Q_cooler']
    print("  --- First law check ---")
    print(f"  Q_heater - Q_cooler : {balance:.2f} kW")
    print(f"  W_net               : {perf['W_net']:.2f} kW")
    print(f"  Discrepancy         : {abs(perf['W_net'] - balance):.4f} kW")
    print("========================")


def print_efficiency_compact(perf):
    print(f"  eta : {perf['eta'] * 100:.2f} %")

#%%
def scale_to_power(W_net_target_MW, perf, PRINT_SCALE=True, First_Law_Check=False):
    """Scale all quantities to target net power. Reference run at m_dot_CO2 = 1 kg/s."""
    scale           = (W_net_target_MW * 1000) / perf['W_net']
    W_net_scaled    = perf['W_net']    * scale / 1000
    Q_heater_scaled = perf['Q_heater'] * scale / 1000
    Q_cooler_scaled = perf['Q_cooler'] * scale / 1000

    if PRINT_SCALE:
        print("=== Scaling to target net power ===")
        print(f"  W_net_target       : {W_net_target_MW:.2f} MW")
        print(f"  W_net_ref (1 kg/s) : {perf['W_net']:.2f} kW")
        print(f"  Scale factor       : {scale:.4f}")
        print(f"  m_dot_CO2 required : {scale:.4f} kg/s")
        print(f"  W_comp             : {perf['W_comp']   * scale / 1000:.2f} MW")
        if perf.get('W_recomp', 0) > 0:
            print(f"  W_recomp           : {perf['W_recomp'] * scale / 1000:.2f} MW")
        print(f"  W_exp              : {perf['W_exp']    * scale / 1000:.2f} MW")
        print(f"  W_net              : {W_net_scaled:.2f} MW")
        print(f"  Q_heater           : {Q_heater_scaled:.2f} MW")
        print(f"  Q_cooler           : {Q_cooler_scaled:.2f} MW")
        if 'Q_recup'    in perf: print(f"  Q_recup            : {perf['Q_recup']    * scale / 1000:.2f} MW")
        if 'Q_recup_lt' in perf:
            print(f"  Q_recup_LT         : {perf['Q_recup_lt'] * scale / 1000:.2f} MW")
            print(f"  Q_recup_HT         : {perf['Q_recup_ht'] * scale / 1000:.2f} MW")
        print(f"  eta (thermo)       : {perf['eta'] * 100:.2f} %")
        if First_Law_Check:
            print("  --- First law check ---")
            print(f"  Q_heater - Q_cooler : {Q_heater_scaled - Q_cooler_scaled:.4f} MW")
            print(f"  W_net               : {W_net_scaled:.4f} MW")
            print(f"  Discrepancy         : {abs(W_net_scaled-(Q_heater_scaled-Q_cooler_scaled)):.6f} MW")
        print("====================================")
    return scale

#%%
def print_net_electric(perf, scale, W_net_target_MW=5.0):
    """
    Apply manufacturer losses: Generator 97%, Mechanical 6%, AirCooler 115 kW, BOP 55.1 kW.
    Scales to W_net_elec = W_net_target_MW for a correct m_dot comparison.
    """
    eta_gen        = 0.97
    eta_mech       = 0.94
    W_parasitic_kW = 115 + 55.1

    W_net_thermo_target = (W_net_target_MW * 1000 + W_parasitic_kW) / (eta_gen * eta_mech)
    scale_elec          = W_net_thermo_target / perf['W_net']
    Q_heater_elec       = perf['Q_heater'] * scale_elec / 1000
    eta_net_elec        = W_net_target_MW / Q_heater_elec

    print("=== Net electrical efficiency ===")
    print(f"  W_net thermo needed    : {W_net_thermo_target:.2f} kW")
    print(f"  m_dot_CO2 for {W_net_target_MW:.0f} MWe  : {scale_elec:.2f} kg/s")
    print(f"  Q_heater               : {Q_heater_elec:.3f} MW")
    print(f"  eta net electrical     : {eta_net_elec*100:.2f}%  (manufacturer = 37.2%)")
    print("=================================")
    return scale_elec


#%%
if __name__ == "__main__":

    study_case      = "Recompression"   # "Simple" | "Recuperated" | "Recompression"
    PRINT           = True
    DETAIL          = True
    PRINT_SCALE     = True
    First_Law_Check = False
    VALIDATION      = False
    W_net_target    = 5                 # [MW]

    # HX discretisation
    PINCH_RECUP = 2    # [K]
    N_DISC      = 20

    # Pressures from Hanwha
    P_low_guess = 89.6  * 1e5   # [Pa]
    P_high      = 215.6 * 1e5   # [Pa]

    # CO2 reference mass flow rate for scaling
    m_dot_CO2 = 1   # [kg/s]

    # efficiencies
    eta_tb = 0.9
    eta_cp = 0.8
    eta_rc = 0.8

    # Split ratio from PFD (Hanwha): m_dot_MC = 47.68, m_dot_total = 68.11 kg/s
    split_ratio = 47.68 / 68.11

    T_hot_su     = 600 + 273.15   # [K] 
    T_salt_limit = 600 + 273.15   # [K]

    if T_hot_su <= T_salt_limit:
        HSource      = SolarSaltConnector()
        HSource.set_properties(T=T_hot_su, p=2e5, m_dot=500)
        PINCH_HEATER = 1.2                  # pinch_min & eta_max, la première contrainte atteinte est limitante, donc si eta = 0.8 et à ce moment là le pinch est > 1.2, alors le solver s'arrete là
        eta_heater   = 0.999                #
        print(f"Hot source : Solar Salt at {T_hot_su - 273.15:.1f}°C")
    else:
        HSource      = MassConnector()
        HSource.set_properties(fluid='air', T=T_hot_su, p=1.2e5, m_dot=500)
        PINCH_HEATER = 5
        eta_heater   = 0.85
        print(f"Hot source : Air at {T_hot_su - 273.15:.1f}°C")

    
    T_exp_su_guess = T_hot_su - PINCH_HEATER   # [K] 
    T_exp_ex_guess = T_exp_su_guess - 100       # [K] 

    T_cold_su    = 37.1 + 273.15   # [K]
    T_c_su_guess = 37.1 + 273.15   # [K]
    eta_cooler   = 0.999

    CSource = MassConnector()
    CSource.set_properties(fluid='air', T=T_cold_su, p=101325, m_dot=500)


    if study_case == "Simple":
        cycle, compressor, expander, gas_heater, gas_cooler = \
            sCO2_basic(eta_cp, eta_tb, eta_heater, eta_cooler,
                       HSource, CSource, P_low_guess, P_high, T_c_su_guess, m_dot_CO2,
                       pinch_heater=PINCH_HEATER, n_disc=N_DISC,
                       T_exp_su_guess=T_exp_su_guess, T_exp_ex_guess=T_exp_ex_guess)
        cycle.solve()
        cycle.plot_cycle_Ts()
        perf = compute_cycle_performance(compressor, expander, gas_heater, gas_cooler)
        if PRINT:
            if DETAIL: print_states(compressor, expander, gas_cooler)
            print_efficiency(perf)
        scale = scale_to_power(W_net_target, perf, PRINT_SCALE, First_Law_Check)
        #print_net_electric(perf, scale, W_net_target)

    elif study_case == "Recuperated":
        eta_recuperator = 0.85
        cycle, compressor, expander, gas_heater, gas_cooler, recuperator = \
            sCO2_recuperated(eta_cp, eta_tb, eta_heater, eta_cooler, eta_recuperator,
                             HSource, CSource, P_low_guess, P_high, T_c_su_guess, m_dot_CO2,
                             pinch_heater=PINCH_HEATER, pinch_recup=PINCH_RECUP, n_disc=N_DISC,
                             T_exp_su_guess=T_exp_su_guess, T_exp_ex_guess=T_exp_ex_guess)
        cycle.solve()
        cycle.plot_cycle_Ts()
        perf = compute_cycle_performance(compressor, expander, gas_heater, gas_cooler,
                                         recuperator=recuperator)
        if PRINT:
            if DETAIL: print_states(compressor, expander, gas_cooler, recuperator=recuperator)
            print_efficiency(perf)
        scale = scale_to_power(W_net_target, perf, PRINT_SCALE, First_Law_Check)
        #print_net_electric(perf, scale, W_net_target)

    elif study_case == "Recompression":
        eta_recup_lt = 0.85
        eta_recup_ht = 0.85

        cycle, compressor, recompressor, expander, gas_heater, gas_cooler, \
            recup_lt, recup_ht, spliter, mixer = \
            sCO2_recompression(eta_cp, eta_rc, eta_tb, eta_heater, eta_cooler,
                               eta_recup_lt, eta_recup_ht, split_ratio,
                               HSource, CSource, P_low_guess, P_high, T_c_su_guess, m_dot_CO2,
                               pinch_heater=PINCH_HEATER, pinch_recup=PINCH_RECUP, n_disc=N_DISC,
                               T_exp_su_guess=T_exp_su_guess, T_exp_ex_guess=T_exp_ex_guess)
        cycle.solve()
        cycle.plot_cycle_Ts()
        perf = compute_cycle_performance(compressor, expander, gas_heater, gas_cooler,
                                         recompressor=recompressor,
                                         recup_lt=recup_lt, recup_ht=recup_ht)
        if PRINT:
            if DETAIL:
                print_states(compressor, expander, gas_cooler,
                             recompressor=recompressor,
                             recup_lt=recup_lt, recup_ht=recup_ht)
            print_efficiency(perf)
        scale = scale_to_power(W_net_target, perf, PRINT_SCALE, First_Law_Check)
        print_net_electric(perf, scale, W_net_target)

        if VALIDATION:
            print("\n=== Validation vs manufacturer ===")
            print(f"  Compressor:su  T = {compressor.su.T-273.15:.1f}°C   (P7  = 37.1°C)")
            print(f"  Expander:su    T = {expander.su.T-273.15:.1f}°C  (P1  = 598.8°C)")
            print(f"  Expander:ex    T = {expander.ex.T-273.15:.1f}°C  (P2  = 488.4°C)")
            print(f"  RecupHT:ex_H   T = {recup_ht.ex_H.T-273.15:.1f}°C  (P4  = 202.4°C)")
            print(f"  RecupLT:ex_H   T = {recup_lt.ex_H.T-273.15:.1f}°C  (P5  = ~86°C)")
            print(f"  RecupLT:ex_C   T = {recup_lt.ex_C.T-273.15:.1f}°C  (P11 = 183.7°C)")
            print(f"  RecupHT:ex_C   T = {recup_ht.ex_C.T-273.15:.1f}°C  (P15 = 440.9°C)")
            if T_hot_su <= T_salt_limit:
                T_salt_ex = SolarSaltConnector._T_from_h(
                    HSource.h - (perf['Q_heater'] * scale * 1000) / (HSource.m_dot * scale))
                print(f"  T_salt_out     = {T_salt_ex - 273.15:.1f}°C  (P14 = 440.9°C)")
            print("==================================")