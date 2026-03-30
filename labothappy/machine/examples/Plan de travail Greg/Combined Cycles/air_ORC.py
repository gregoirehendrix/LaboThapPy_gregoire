# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 2026

@author: gregoire.hendrix

Combined Open Air Brayton + ORC cycle

Available configurations (choose in __main__):
  study_case_brayton : "Simple" | "Recuperated"
  study_case_orc     : "Simple" | "Recuperated"

Brayton topologies:
  Simple      : AirInlet → Compressor → Heater → Turbine_Br → ORC evap. (hot)
  Recuperated : AirInlet → Compressor → Recuperator (cold) → Heater → Turbine_Br
                         → Recuperator (hot) → ORC evap. (hot)

ORC topologies:
  Simple      : Pump → Evaporator → Expander → Condenser → Pump
  Recuperated : Pump → Recup_ORC (cold) → Evaporator → Expander
                     → Recup_ORC (hot)  → Condenser  → Pump
"""

# =============================================================================
# IMPORTS
# =============================================================================
from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.machine.circuit_it  import IterativeCircuit
from labothappy.connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI

from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff     import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff    import HexCstEff
from labothappy.component.pump.pump_csteff             import PumpCstEff


#%%
# =============================================================================
# SHARED PARAMETERS
# =============================================================================
# --- Brayton ---
fluid_br    = 'Air'
T_amb       = 13.2  + 273.15    # K  — ambient / compressor inlet
P_atm       = 1.01325e5         # Pa
PR          = 11                # pressure ratio [-]
eta_comp    = 0.80              # compressor isentropic efficiency [-]
eta_turb_br = 0.90              # Brayton turbine isentropic efficiency [-]
eta_heater  = 1.0               # heater effectiveness [-]
eta_recup   = 0.85              # Brayton recuperator effectiveness [-]

T_hot_su    = 1000 + 273.15     # K  — hot source temperature
P_hot       = 1e5               # Pa
m_dot_hot   = 10.0              # kg/s

# --- ORC ---
fluid_orc    = 'R245fa'
eta_pump_orc = 0.75
eta_exp_orc  = 0.85
eta_evap_orc = 0.90
eta_cond_orc = 0.90
eta_recup_orc= 0.85             # ORC recuperator effectiveness [-]

P_high_orc   = 20e5            # Pa — ORC high pressure
P_low_orc    = 2e5             # Pa — ORC low pressure
m_dot_orc    = 1.0             # kg/s

# Cold source (condenser)
T_cold_su   = T_amb
P_cold      = P_atm
m_dot_cold  = 10.0             # kg/s



#%%
# =============================================================================
# BRAYTON FUNCTIONS
# =============================================================================
def brayton_simple():
    """Open Air Brayton — no recuperator. Sequential solve via IterativeCircuit."""

    cycle      = IterativeCircuit(fluid=fluid_br)
    compressor = CompressorCstEff()
    heater     = HexCstEff()
    turbine    = ExpanderCstEff()

    compressor.set_parameters(eta_is=eta_comp)
    turbine.set_parameters(eta_is=eta_turb_br)
    heater.set_parameters(eta=eta_heater)

    cycle.add_component(compressor, "Compressor")
    cycle.add_component(heater,     "Heater")
    cycle.add_component(turbine,    "Turbine_Br")

    cycle.link_components("Compressor", "m-ex",   "Heater",     "m-su_C")
    cycle.link_components("Heater",     "m-ex_C", "Turbine_Br", "m-su")

    air_in = MassConnector()
    cycle.add_source("AirInlet", air_in, cycle.components["Compressor"], "m-su")
    cycle.set_source_properties(target="AirInlet",
                                fluid=fluid_br, T=T_amb, P=P_atm, m_dot=1.0)

    hot_source = MassConnector()
    cycle.add_source("HotSource", hot_source, cycle.components["Heater"], "m-su_H")
    cycle.set_source_properties(target="HotSource",
                                fluid=fluid_br, T=T_hot_su, P=P_hot, m_dot=m_dot_hot)

    cycle.set_cycle_input(target="Compressor:ex", p=P_atm * PR)
    cycle.set_cycle_input(target="Turbine_Br:ex", p=P_atm)

    cycle._build_solve_order()
    for name in cycle.solve_start_components:
        cycle.components[name].solve()

    return cycle, compressor, heater, turbine, None   # None = no recuperator

#%%
def brayton_recuperated():
    """Open Air Brayton with recuperator. RecursiveCircuit."""

    cycle       = RecursiveCircuit(fluid_br)
    compressor  = CompressorCstEff()
    heater      = HexCstEff()
    turbine     = ExpanderCstEff()
    recuperator = HexCstEff()

    compressor.set_parameters(eta_is=eta_comp)
    turbine.set_parameters(eta_is=eta_turb_br)
    heater.set_parameters(eta=eta_heater)
    recuperator.set_parameters(eta=eta_recup)

    cycle.add_component(compressor,  "Compressor")
    cycle.add_component(recuperator, "Recuperator")
    cycle.add_component(heater,      "Heater")
    cycle.add_component(turbine,     "Turbine_Br")

    cycle.link_components("Compressor",  "m-ex",   "Recuperator", "m-su_C")
    cycle.link_components("Recuperator", "m-ex_C", "Heater",      "m-su_C")
    cycle.link_components("Heater",      "m-ex_C", "Turbine_Br",  "m-su")
    cycle.link_components("Turbine_Br",  "m-ex",   "Recuperator", "m-su_H")

    air_in = MassConnector()
    cycle.add_source("AirInlet", air_in, cycle.components["Compressor"], "m-su")
    cycle.set_source_properties(target="AirInlet",
                                fluid=fluid_br, T=T_amb, P=P_atm, m_dot=1.0)

    hot_source = MassConnector()
    cycle.add_source("HotSource", hot_source, cycle.components["Heater"], "m-su_H")
    cycle.set_source_properties(target="HotSource",
                                fluid=fluid_br, T=T_hot_su, P=P_hot, m_dot=m_dot_hot)

    cycle.set_fixed_properties(target="Compressor:ex", p=P_atm * PR)
    cycle.set_fixed_properties(target="Turbine_Br:ex", p=P_atm)

    # Analytical guesses (no solved component yet)
    h_cs    = PropsSI('H', 'T', T_amb,    'P', P_atm,     fluid_br)
    s_cs    = PropsSI('S', 'T', T_amb,    'P', P_atm,     fluid_br)
    h_cx_is = PropsSI('H', 'P', P_atm*PR, 'S', s_cs,      fluid_br)
    h_cx    = h_cs + (h_cx_is - h_cs) / eta_comp
    T_cx    = PropsSI('T', 'H', h_cx,    'P', P_atm*PR,   fluid_br)

    h_ts    = PropsSI('H', 'T', T_hot_su, 'P', P_atm*PR,  fluid_br)
    s_ts    = PropsSI('S', 'T', T_hot_su, 'P', P_atm*PR,  fluid_br)
    h_tx_is = PropsSI('H', 'P', P_atm,   'S', s_ts,       fluid_br)
    h_tx    = h_ts - eta_turb_br * (h_ts - h_tx_is)
    T_tx    = PropsSI('T', 'H', h_tx,    'P', P_atm,      fluid_br)

    T_rec_xC = T_cx + eta_recup * (T_tx - T_cx)
    T_rec_xH = T_tx - eta_recup * (T_tx - T_cx)

    cycle.set_cycle_guess(target='Recuperator:su_C', T=T_cx,      p=P_atm*PR, m_dot=1.0)
    cycle.set_cycle_guess(target='Recuperator:su_H', T=T_tx,      p=P_atm,    m_dot=1.0)
    cycle.set_cycle_guess(target='Recuperator:ex_C', T=T_rec_xC,  p=P_atm*PR)
    cycle.set_cycle_guess(target='Recuperator:ex_H', T=T_rec_xH,  p=P_atm)

    cycle.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance=1e-3)

    cycle.mute_print()
    cycle.solve(max_iter=50)

    return cycle, compressor, heater, turbine, recuperator

#%%
# =============================================================================
# ORC FUNCTIONS
# =============================================================================

def _orc_guesses():
    """Return temperature guesses derived from ORC pressure levels."""
    T_sat_hi = PropsSI('T', 'P', P_high_orc, 'Q', 1, fluid_orc)
    T_sat_lo = PropsSI('T', 'P', P_low_orc,  'Q', 1, fluid_orc)
    return T_sat_hi + 10, T_sat_lo + 20, T_sat_lo - 3
    # returns: T_exp_su_guess, T_exp_ex_guess, T_pump_su_guess

#%%
def orc_simple(T_hot_source, P_hot_source=None):
    """
    Simple ORC — no recuperator.
    Pump → Evaporator → Expander → Condenser → Pump
    """
    if P_hot_source is None:
        P_hot_source = P_atm

    cycle      = RecursiveCircuit(fluid_orc)
    pump       = PumpCstEff()
    evap       = HexCstEff()
    expander   = ExpanderCstEff()
    condenser  = HexCstEff()

    pump.set_parameters(eta_is=eta_pump_orc)
    expander.set_parameters(eta_is=eta_exp_orc)
    evap.set_parameters(eta=eta_evap_orc)
    condenser.set_parameters(eta=eta_cond_orc)

    cycle.add_component(pump,      "Pump")
    cycle.add_component(evap,      "Evaporator")
    cycle.add_component(expander,  "Expander_ORC")
    cycle.add_component(condenser, "Condenser")

    cycle.link_components("Pump",         "m-ex",   "Evaporator",   "m-su_C")
    cycle.link_components("Evaporator",   "m-ex_C", "Expander_ORC", "m-su")
    cycle.link_components("Expander_ORC", "m-ex",   "Condenser",    "m-su_H")
    cycle.link_components("Condenser",    "m-ex_H", "Pump",         "m-su")

    hot = MassConnector()
    cycle.add_source("BraytonExhaust", hot, cycle.components["Evaporator"], "m-su_H")
    cycle.set_source_properties(target="BraytonExhaust",
                                fluid=fluid_br, T=T_hot_source, P=P_hot_source, m_dot=1.0)

    cold = MassConnector()
    cycle.add_source("ColdSource", cold, cycle.components["Condenser"], "m-su_C")
    cycle.set_source_properties(target="ColdSource",
                                fluid=fluid_br, T=T_cold_su, P=P_cold, m_dot=m_dot_cold)

    T_es, T_ex, T_ps = _orc_guesses()

    cycle.set_cycle_guess(target='Pump:su',         p=P_low_orc,  SC=3,       m_dot=m_dot_orc)
    cycle.set_cycle_guess(target='Pump:ex',         p=P_high_orc, T=T_ps+2)
    cycle.set_cycle_guess(target='Expander_ORC:su', p=P_high_orc, T=T_es,     m_dot=m_dot_orc)
    cycle.set_cycle_guess(target='Expander_ORC:ex', p=P_low_orc,  T=T_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',  p=P_low_orc,  T=T_ps+5)

    cycle.set_residual_variable(target='Expander_ORC:ex', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Expander_ORC:ex', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Pump:ex',         variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Pump:ex',         variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Evaporator:ex_C', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Condenser:ex_H',  variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Condenser:ex_H',  variable='p', tolerance=1e-3)

    cycle.mute_print()
    cycle.solve(max_iter=50)

    return cycle, pump, evap, expander, condenser, None

#%% fct pas encore
def orc_recuperated(T_hot_source, P_hot_source=None):
    """
    ORC with recuperator.
    Pump → Recup_ORC (cold) → Evaporator → Expander → Recup_ORC (hot) → Condenser → Pump
    """
    if P_hot_source is None:
        P_hot_source = P_atm

    cycle      = RecursiveCircuit(fluid_orc)
    pump       = PumpCstEff()
    recup      = HexCstEff()
    evap       = HexCstEff()
    expander   = ExpanderCstEff()
    condenser  = HexCstEff()

    pump.set_parameters(eta_is=eta_pump_orc)
    expander.set_parameters(eta_is=eta_exp_orc)
    evap.set_parameters(eta=eta_evap_orc)
    condenser.set_parameters(eta=eta_cond_orc)
    recup.set_parameters(eta=eta_recup_orc)

    cycle.add_component(pump,      "Pump")
    cycle.add_component(recup,     "Recup_ORC")
    cycle.add_component(evap,      "Evaporator")
    cycle.add_component(expander,  "Expander_ORC")
    cycle.add_component(condenser, "Condenser")

    cycle.link_components("Pump",         "m-ex",   "Recup_ORC",    "m-su_C")
    cycle.link_components("Recup_ORC",    "m-ex_C", "Evaporator",   "m-su_C")
    cycle.link_components("Evaporator",   "m-ex_C", "Expander_ORC", "m-su")
    cycle.link_components("Expander_ORC", "m-ex",   "Recup_ORC",    "m-su_H")
    cycle.link_components("Recup_ORC",    "m-ex_H", "Condenser",    "m-su_H")
    cycle.link_components("Condenser",    "m-ex_H", "Pump",         "m-su")

    hot = MassConnector()
    cycle.add_source("BraytonExhaust", hot, cycle.components["Evaporator"], "m-su_H")
    cycle.set_source_properties(target="BraytonExhaust",
                                fluid=fluid_br, T=T_hot_source, P=P_hot_source, m_dot=1.0)

    cold = MassConnector()
    cycle.add_source("ColdSource", cold, cycle.components["Condenser"], "m-su_C")
    cycle.set_source_properties(target="ColdSource",
                                fluid=fluid_br, T=T_cold_su, P=P_cold, m_dot=m_dot_cold)

    T_es, T_ex, T_ps = _orc_guesses()
    T_rxC = T_ps + eta_recup_orc * (T_ex - T_ps)
    T_rxH = T_ex - eta_recup_orc * (T_ex - T_ps)

    cycle.set_cycle_guess(target='Pump:su',         p=P_low_orc,  SC=3,       m_dot=m_dot_orc)
    cycle.set_cycle_guess(target='Pump:ex',         p=P_high_orc, T=T_ps+2)
    cycle.set_cycle_guess(target='Recup_ORC:su_C',  p=P_high_orc, T=T_ps+2,   m_dot=m_dot_orc)
    cycle.set_cycle_guess(target='Recup_ORC:su_H',  p=P_low_orc,  T=T_ex,     m_dot=m_dot_orc)
    cycle.set_cycle_guess(target='Recup_ORC:ex_C',  p=P_high_orc, T=T_rxC)
    cycle.set_cycle_guess(target='Recup_ORC:ex_H',  p=P_low_orc,  T=T_rxH)
    cycle.set_cycle_guess(target='Expander_ORC:su', p=P_high_orc, T=T_es,     m_dot=m_dot_orc)
    cycle.set_cycle_guess(target='Expander_ORC:ex', p=P_low_orc,  T=T_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',  p=P_low_orc,  T=T_ps+5)

    cycle.set_residual_variable(target='Expander_ORC:ex', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Expander_ORC:ex', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Pump:ex',         variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Pump:ex',         variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Evaporator:ex_C', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Condenser:ex_H',  variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Condenser:ex_H',  variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recup_ORC:su_H',  variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recup_ORC:su_H',  variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recup_ORC:ex_C',  variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recup_ORC:ex_C',  variable='p', tolerance=1e-3)

    cycle.mute_print()
    cycle.solve(max_iter=50)

    return cycle, pump, evap, expander, condenser, recup

#%%
# =============================================================================
# PRINT FUNCTIONS
# =============================================================================

def print_brayton(compressor, heater, turbine, recuperator=None):
    tag = " + Recuperator" if recuperator is not None else ""
    print(f"{'='*52}")
    print(f"  Brayton Cycle{tag}")
    print(f"{'='*52}")
    print(f"  PR             = {PR:.1f}  [-]")
    print(f"  TIT            = {turbine.su.T - 273.15:.1f}  [°C]")
    print(f"  T_ex_comp      = {compressor.ex.T - 273.15:.1f}  [°C]")
    print(f"  T_ex_turb      = {turbine.ex.T - 273.15:.1f}  [°C]")
    if recuperator is not None:
        print(f"  T_recup_ex_H   = {recuperator.ex_H.T - 273.15:.1f}  [°C]  ← feeds ORC")
        print(f"  T_recup_ex_C   = {recuperator.ex_C.T - 273.15:.1f}  [°C]  (heater inlet)")
        print(f"  Q_recup        = {recuperator.Q.Q_dot/1e3:.2f}  [kW/kg/s air]")
    W_net = turbine.W.W_dot - compressor.W.W_dot
    Q_in  = heater.Q.Q_dot
    print(f"  W_comp         = {compressor.W.W_dot/1e3:.2f}  [kW/kg/s air]")
    print(f"  W_turb         = {turbine.W.W_dot/1e3:.2f}  [kW/kg/s air]")
    print(f"  W_net          = {W_net/1e3:.2f}  [kW/kg/s air]")
    print(f"  Q_in           = {Q_in/1e3:.2f}  [kW/kg/s air]")
    print(f"  eta_Br         = {W_net/Q_in*100:.1f}  [%]")

#%%
def print_combined(compressor, heater, turbine, recuperator,
                   pump, evap, expander, condenser, recup_orc):

    br_tag  = " + Recuperator" if recuperator is not None else ""
    orc_tag = " + Recuperator" if recup_orc   is not None else ""

    W_net_br  = turbine.W.W_dot - compressor.W.W_dot
    Q_in_br   = heater.Q.Q_dot
    eta_br    = W_net_br / Q_in_br

    W_net_orc = expander.W.W_dot - pump.W.W_dot
    Q_evap    = evap.Q.Q_dot
    eta_orc   = W_net_orc / Q_evap

    W_net_comb = W_net_br + W_net_orc
    eta_comb   = W_net_comb / Q_in_br

    print(f"\n{'='*52}")
    print(f"  Combined — Brayton{br_tag} + ORC{orc_tag}")
    print(f"{'='*52}")
    print(f"  {'─'*48}")
    print(f"  BRAYTON{br_tag}")
    print(f"  {'─'*48}")
    print(f"  W_net_Br       = {W_net_br/1e3:7.2f}  [kW/kg/s air]")
    print(f"  Q_in_Br        = {Q_in_br/1e3:7.2f}  [kW/kg/s air]")
    print(f"  eta_Br         = {eta_br*100:7.1f}  [%]")
    print(f"  T_ex_turb_Br   = {turbine.ex.T - 273.15:7.1f}  [°C]")
    if recuperator is not None:
        print(f"  T_ORC_hot_su   = {recuperator.ex_H.T - 273.15:7.1f}  [°C]  (→ ORC evap.)")
    print(f"  {'─'*48}")
    print(f"  ORC ({fluid_orc}){orc_tag}")
    print(f"  {'─'*48}")
    if recup_orc is not None:
        print(f"  T_recup_ex_C   = {recup_orc.ex_C.T - 273.15:7.1f}  [°C]  (evap. inlet)")
        print(f"  T_recup_ex_H   = {recup_orc.ex_H.T - 273.15:7.1f}  [°C]  (cond. inlet)")
        print(f"  Q_recup_ORC    = {recup_orc.Q.Q_dot/1e3:7.2f}  [kW/kg/s ORC]")
    print(f"  T_evap_su      = {evap.su_C.T - 273.15:7.1f}  [°C]  (evap. cold inlet)")
    print(f"  T_evap_ex      = {evap.ex_C.T - 273.15:7.1f}  [°C]  (expander inlet)")
    print(f"  T_exp_ex       = {expander.ex.T - 273.15:7.1f}  [°C]")
    print(f"  T_cond_ex      = {condenser.ex_H.T - 273.15:7.1f}  [°C]  (pump inlet)")
    print(f"  P_high_ORC     = {expander.su.p/1e5:7.2f}  [bar]")
    print(f"  P_low_ORC      = {pump.su.p/1e5:7.2f}  [bar]")
    print(f"  W_exp_ORC      = {expander.W.W_dot/1e3:7.2f}  [kW/kg/s ORC]")
    print(f"  W_pump_ORC     = {pump.W.W_dot/1e3:7.2f}  [kW/kg/s ORC]")
    print(f"  W_net_ORC      = {W_net_orc/1e3:7.2f}  [kW/kg/s ORC]")
    print(f"  Q_evap_ORC     = {Q_evap/1e3:7.2f}  [kW/kg/s ORC]")
    print(f"  eta_ORC        = {eta_orc*100:7.1f}  [%]")
    print(f"  {'─'*48}")
    print(f"  COMBINED")
    print(f"  {'─'*48}")
    print(f"  W_net_combined = {W_net_comb/1e3:7.2f}  [kW/kg/s air]")
    print(f"  Q_in           = {Q_in_br/1e3:7.2f}  [kW/kg/s air]")
    print(f"  eta_combined   = {eta_comb*100:7.1f}  [%]")
    print(f"  Gain vs Br     = +{(eta_comb - eta_br)*100:5.1f}  pp")
    print(f"  {'='*48}")

#%%
# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":

    study_case_brayton = "Simple"   # "Simple" | "Recuperated"
    study_case_orc     = "Simple"   # "Simple" | "Recuperated"

    # --- Step 1: Brayton ---
    if study_case_brayton == "Simple":
        cycle_br, compressor, heater, turbine, recuperator = brayton_simple()
        T_orc_hot = turbine.ex.T
    elif study_case_brayton == "Recuperated":
        cycle_br, compressor, heater, turbine, recuperator = brayton_recuperated()
        T_orc_hot = recuperator.ex_H.T

    print_brayton(compressor, heater, turbine, recuperator)

    # --- Step 2: ORC ---
    if study_case_orc == "Simple":
        cycle_orc, pump, evap, expander, condenser, recup_orc = \
            orc_simple(T_orc_hot)
    elif study_case_orc == "Recuperated":
        cycle_orc, pump, evap, expander, condenser, recup_orc = \
            orc_recuperated(T_orc_hot)

    orc_ok = cycle_orc.converged or all(
        rv.converged for rv in cycle_orc.res_vars.values())

    # --- Step 3: Combined results ---
    if orc_ok:
        print_combined(compressor, heater, turbine, recuperator,
                       pump, evap, expander, condenser, recup_orc)
    else:
        print("\n  ORC did not converge. Check initial guesses and parameters.")