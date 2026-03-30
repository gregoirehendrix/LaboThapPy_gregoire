# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:38:29 2026

@author: gregoire.hendrix

Basic Steam Rankine Cycle — RecursiveCircuit
Topology : Pump -> SteamGenerator (HexCstEff) -> Expander -> Condenser (HexCstEff)
Hot source : thermal oil (Cp = 2200 J/kg/K, custom fluid)
Cold sink  : cooling water
Working fluid : Water
"""

#%%
import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.pump.pump_csteff import PumpCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

#%%
def make_oil_connector(T_K, p_Pa, m_dot_kgs, Cp_Jkg):
    """
    Build a MassConnector for a custom oil (not in CoolProp).
    Attributes are assigned directly to bypass CoolProp lookup.
    HexCstEff Path 3 reads su_H.cp, su_H.T, su_H.h directly.
    """
    oil = MassConnector()
    # bypass set_properties — 'OilHTF' is not a CoolProp fluid
    oil.fluid        = 'OilHTF'
    oil.T            = T_K
    oil.p            = p_Pa
    oil.m_dot        = m_dot_kgs
    oil.h            = Cp_Jkg * T_K   # reference: h = 0 at 0 K
    oil.cp           = Cp_Jkg
    oil.state_known  = True
    oil.completely_known = True
    return oil

#%%
def rankine_basic(eta_pump, eta_turb, eta_sg, eta_cd,
                  HSource, CSource,
                  P_low, P_high, m_dot_wf,
                  T_pump_su_guess=None, T_turb_su_guess=None):

    if T_pump_su_guess is None:
        T_pump_su_guess = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    if T_turb_su_guess is None:
        T_turb_su_guess = PropsSI('T', 'P', P_high, 'Q', 1, 'Water') + 20

    cycle = RecursiveCircuit('Water')

    pump     = PumpCstEff()
    sg       = HexCstEff()
    expander = ExpanderCstEff()
    cd       = HexCstEff()

    pump.set_parameters(eta_is=eta_pump)
    sg.set_parameters(eta=eta_sg)
    expander.set_parameters(eta_is=eta_turb)
    cd.set_parameters(eta=eta_cd)

    cycle.add_component(pump,     "Pump")
    cycle.add_component(sg,       "SG")
    cycle.add_component(expander, "Expander")
    cycle.add_component(cd,       "CD")

    # Working fluid loop
    cycle.link_components("Pump",     "m-ex",   "SG",       "m-su_C")
    cycle.link_components("SG",       "m-ex_C", "Expander", "m-su")
    cycle.link_components("Expander", "m-ex",   "CD",       "m-su_H")
    cycle.link_components("CD",       "m-ex_H", "Pump",     "m-su")

    # External sources
    cycle.add_source("Hot_source",  HSource, cycle.components["SG"], "m-su_H")
    cycle.add_source("Cold_source", CSource, cycle.components["CD"], "m-su_C")

    # Pump inlet: saturated liquid at P_low — fully defines the starting state
    cycle.set_cycle_guess(target='Pump:su', m_dot=m_dot_wf, T=T_pump_su_guess, p=P_low)

    # Pump outlet pressure (fixed — not overwritten by solver iterations)
    cycle.set_fixed_properties(target='Pump:ex', p=P_high)

    # Expander inlet guess
    cycle.set_cycle_guess(target='Expander:su', m_dot=m_dot_wf, T=T_turb_su_guess, p=P_high)

    # Expander outlet pressure (fixed)
    cycle.set_fixed_properties(target='Expander:ex', p=P_low)

    return cycle, pump, expander, sg, cd

#%%
def compute_cycle_performance(pump, expander, sg, cd):
    m_dot  = expander.su.m_dot

    W_pump = m_dot * (pump.ex.h     - pump.su.h)     / 1000   # kW
    W_exp  = m_dot * (expander.su.h - expander.ex.h) / 1000   # kW
    W_net  = W_exp - W_pump

    Q_sg   = m_dot * (sg.ex_C.h - sg.su_C.h)  / 1000   # kW
    Q_cd   = m_dot * (cd.su_H.h - cd.ex_H.h)  / 1000   # kW

    eta    = W_net / Q_sg

    return {
        'W_pump': W_pump, 'W_exp': W_exp, 'W_net': W_net,
        'Q_sg': Q_sg, 'Q_cd': Q_cd, 'eta': eta,
    }

#%%
def print_states(pump, expander, sg, cd):
    print("=== Pump:su ===")
    print(f"  T     = {pump.su.T - 273.15:.2f} °C")
    print(f"  p     = {pump.su.p / 1e5:.2f} bar")
    print(f"  h     = {pump.su.h:.2f} J/kg")
    print(f"  m_dot = {pump.su.m_dot:.4f} kg/s")
    print("=== Pump:ex ===")
    print(f"  T     = {pump.ex.T - 273.15:.2f} °C")
    print(f"  p     = {pump.ex.p / 1e5:.2f} bar")
    print(f"  h     = {pump.ex.h:.2f} J/kg")
    print("=== SG — water side (su_C / ex_C) ===")
    print(f"  su_C T = {sg.su_C.T - 273.15:.2f} °C")
    print(f"  ex_C T = {sg.ex_C.T - 273.15:.2f} °C")
    print(f"  ex_C p = {sg.ex_C.p / 1e5:.2f} bar")
    print("=== SG — oil side (su_H / ex_H) ===")
    print(f"  su_H T = {sg.su_H.T - 273.15:.2f} °C")
    print(f"  ex_H T = {sg.ex_H.T - 273.15:.2f} °C")
    print("=== Expander:su ===")
    print(f"  T     = {expander.su.T - 273.15:.2f} °C")
    print(f"  p     = {expander.su.p / 1e5:.2f} bar")
    print(f"  h     = {expander.su.h:.2f} J/kg")
    print("=== Expander:ex ===")
    print(f"  T     = {expander.ex.T - 273.15:.2f} °C")
    print(f"  p     = {expander.ex.p / 1e5:.2f} bar")
    print(f"  h     = {expander.ex.h:.2f} J/kg")
    print("=== CD — steam side (su_H / ex_H) ===")
    print(f"  su_H T = {cd.su_H.T - 273.15:.2f} °C")
    print(f"  ex_H T = {cd.ex_H.T - 273.15:.2f} °C")
    print("=== CD — cooling water (su_C / ex_C) ===")
    print(f"  su_C T = {cd.su_C.T - 273.15:.2f} °C")
    print(f"  ex_C T = {cd.ex_C.T - 273.15:.2f} °C")

#%%
def print_efficiency(perf):
    print("=== Cycle Performance ===")
    print(f"  W_pump   : {perf['W_pump']:.2f} kW")
    print(f"  W_exp    : {perf['W_exp']:.2f} kW")
    print(f"  W_net    : {perf['W_net']:.2f} kW")
    print(f"  Q_sg     : {perf['Q_sg']:.2f} kW")
    print(f"  Q_cd     : {perf['Q_cd']:.2f} kW")
    print(f"  eta      : {perf['eta'] * 100:.2f} %")
    balance = perf['Q_sg'] - perf['Q_cd']
    print("  --- First law check ---")
    print(f"  Q_sg - Q_cd : {balance:.2f} kW")
    print(f"  W_net       : {perf['W_net']:.2f} kW")
    print(f"  Discrepancy : {abs(perf['W_net'] - balance):.4f} kW")
    print("=========================")

#%%
def scale_to_power(W_net_target_MW, perf, PRINT_SCALE=True):
    """Scale all quantities to target net power. Reference run at m_dot_wf = 1 kg/s."""
    scale        = (W_net_target_MW * 1000) / perf['W_net']
    W_net_scaled = perf['W_net'] * scale / 1000
    Q_sg_scaled  = perf['Q_sg']  * scale / 1000
    Q_cd_scaled  = perf['Q_cd']  * scale / 1000

    if PRINT_SCALE:
        print("=== Scaling to target net power ===")
        print(f"  W_net_target      : {W_net_target_MW:.2f} MW")
        print(f"  W_net_ref (1kg/s) : {perf['W_net']:.2f} kW")
        print(f"  Scale factor      : {scale:.4f}")
        print(f"  m_dot_wf required : {scale:.4f} kg/s")
        print(f"  W_pump            : {perf['W_pump'] * scale / 1000:.2f} MW")
        print(f"  W_exp             : {perf['W_exp']  * scale / 1000:.2f} MW")
        print(f"  W_net             : {W_net_scaled:.2f} MW")
        print(f"  Q_sg              : {Q_sg_scaled:.2f} MW")
        print(f"  Q_cd              : {Q_cd_scaled:.2f} MW")
        print(f"  eta (thermo)      : {perf['eta'] * 100:.2f} %")
        print("====================================")
    return scale

#%%
if __name__ == "__main__":

    PRINT        = True
    DETAIL       = True
    PRINT_SCALE  = True
    W_net_target = 5   # [MW]

    # Cycle pressures
    P_low  = 0.10e5   # [Pa] — condensation
    P_high = 50e5     # [Pa] — evaporation

    # Working fluid reference mass flow rate (scaled afterwards)
    m_dot_wf = 1.0   # [kg/s]

    # Efficiencies
    eta_pump = 0.75
    eta_turb = 0.85
    eta_sg   = 0.90
    eta_cd   = 0.90

    # Thermal oil (hot source)
    Cp_oil    = 2200            # [J/kg/K]
    T_oil_su  = 300 + 273.15   # [K]
    m_dot_oil = 2.0             # [kg/s]
    P_oil     = 2e5             # [Pa]

    HSource = make_oil_connector(T_oil_su, P_oil, m_dot_oil, Cp_oil)

    # Cooling water (cold sink)
    T_cw_su  = 20 + 273.15   # [K]
    P_cw     = 2e5            # [Pa]
    m_dot_cw = 5.0            # [kg/s]

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=T_cw_su, p=P_cw, m_dot=m_dot_cw)

    # Temperature guesses
    T_pump_su_guess = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_turb_su_guess = PropsSI('T', 'P', P_high, 'Q', 1, 'Water') + 20

    cycle, pump, expander, sg, cd = rankine_basic(
        eta_pump, eta_turb, eta_sg, eta_cd,
        HSource, CSource,
        P_low, P_high, m_dot_wf,
        T_pump_su_guess=T_pump_su_guess,
        T_turb_su_guess=T_turb_su_guess,
    )

    cycle.solve()

    perf = compute_cycle_performance(pump, expander, sg, cd)

    if PRINT:
        if DETAIL:
            print_states(pump, expander, sg, cd)
        print_efficiency(perf)

    scale = scale_to_power(W_net_target, perf, PRINT_SCALE)