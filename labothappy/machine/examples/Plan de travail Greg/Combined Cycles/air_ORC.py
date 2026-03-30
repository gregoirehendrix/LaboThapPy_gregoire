# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 2026

@author: gregoire.hendrix

Combined Open Air Brayton + ORC cycle (first version — no recuperator)

Topology:
  AIR BRAYTON (open, sequential solve):
    AirInlet → Compressor → Heater (hot air source) → Turbine_Br → Evaporator (hot side) → exhaust

  ORC (R245fa, RecursiveCircuit):
    Pump → Evaporator (cold side) → Expander_ORC → Condenser (hot side) → Pump
    Cold source: ambient air (MassConnector)

Connection:
  Turbine_Br exhaust air feeds the hot side of the ORC Evaporator directly.
"""

# =============================================================================
# IMPORTS
# =============================================================================
from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.machine.circuit_it import IterativeCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector

from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from labothappy.component.pump.pump_csteff import PumpCstEff

# =============================================================================
# PARAMETERS
# =============================================================================

# --- Brayton ---
fluid_br    = 'Air'
T_amb       = 13.2 + 273.15     # K  — compressor inlet / cold source
P_atm       = 1.01325e5         # Pa
PR          = 11                # pressure ratio
eta_comp    = 0.80              # compressor isentropic efficiency
eta_turb_br = 0.90              # Brayton turbine isentropic efficiency
eta_heater  = 1.0               # heater effectiveness

T_hot_su    = 1000 + 273.15     # K  — hot air source temperature
P_hot       = 1e5               # Pa — hot source pressure
m_dot_hot   = 10.0              # kg/s — hot source mass flow rate (arbitrary, scaled later)

# --- ORC ---
fluid_orc       = 'R245fa'
eta_pump_orc    = 0.75
eta_exp_orc     = 0.85
eta_evap_orc    = 0.90          # evaporator effectiveness
eta_cond_orc    = 0.90          # condenser effectiveness

P_high_orc      = 20e5         # Pa — ORC high pressure (initial guess)
P_low_orc       = 2e5          # Pa — ORC low pressure (initial guess)
m_dot_orc       = 1.0          # kg/s — ORC working fluid mass flow rate

# Cold source (condenser)
T_cold_su   = T_amb             # K  — ambient air temperature
P_cold      = P_atm             # Pa
m_dot_cold  = 10.0              # kg/s — cold air mass flow rate

# =============================================================================
# STEP 1 — BRAYTON CYCLE (sequential, open cycle)
# =============================================================================
print("=" * 52)
print("  Step 1 — Open Air Brayton Cycle")
print("=" * 52)

brayton     = IterativeCircuit(fluid=fluid_br)
compressor  = CompressorCstEff()
heater      = HexCstEff()
turbine_br  = ExpanderCstEff()

compressor.set_parameters(eta_is=eta_comp)
turbine_br.set_parameters(eta_is=eta_turb_br)
heater.set_parameters(eta=eta_heater)

brayton.add_component(compressor, "Compressor")
brayton.add_component(heater,     "Heater")
brayton.add_component(turbine_br, "Turbine_Br")

brayton.link_components("Compressor", "m-ex",   "Heater",     "m-su_C")
brayton.link_components("Heater",     "m-ex_C", "Turbine_Br", "m-su")

# Air inlet source
air_in = MassConnector()
brayton.add_source("AirInlet", air_in, brayton.components["Compressor"], "m-su")
brayton.set_source_properties(target="AirInlet",
                               fluid=fluid_br, T=T_amb, P=P_atm, m_dot=1.0)

# Hot source (hot air)
hot_source = MassConnector()
brayton.add_source("HotSource", hot_source, brayton.components["Heater"], "m-su_H")
brayton.set_source_properties(target="HotSource",
                               fluid=fluid_br, T=T_hot_su, P=P_hot, m_dot=m_dot_hot)

# Fixed pressures
brayton.set_cycle_input(target="Compressor:ex", p=P_atm * PR)
brayton.set_cycle_input(target="Turbine_Br:ex", p=P_atm)

# Sequential solve
brayton._build_solve_order()
for name in brayton.solve_start_components:
    brayton.components[name].solve()

# Brayton results
W_comp   = compressor.W.W_dot
W_turb_br = turbine_br.W.W_dot
W_net_br  = W_turb_br - W_comp
Q_in_br   = heater.Q.Q_dot
eta_br    = W_net_br / Q_in_br

T_br_ex   = turbine_br.ex.T   # K — Brayton exhaust temperature (feeds ORC evaporator)
h_br_ex   = turbine_br.ex.h   # J/kg

print(f"  PR             = {PR:.1f}  [-]")
print(f"  TIT            = {turbine_br.su.T - 273.15:.1f}  [°C]")
print(f"  T_ex_comp      = {compressor.ex.T - 273.15:.1f}  [°C]")
print(f"  T_ex_turb_Br   = {T_br_ex - 273.15:.1f}  [°C]  ← feeds ORC evaporator")
print(f"  W_comp         = {W_comp/1e3:.2f}  [kW per kg/s air]")
print(f"  W_turb_Br      = {W_turb_br/1e3:.2f}  [kW per kg/s air]")
print(f"  W_net_Br       = {W_net_br/1e3:.2f}  [kW per kg/s air]")
print(f"  Q_in_Br        = {Q_in_br/1e3:.2f}  [kW per kg/s air]")
print(f"  eta_Br         = {eta_br*100:.1f}  [%]")

# =============================================================================
# STEP 2 — ORC CYCLE (RecursiveCircuit)
# =============================================================================
print()
print("=" * 52)
print("  Step 2 — ORC Cycle (R245fa)")
print("=" * 52)

orc          = RecursiveCircuit(fluid_orc)
pump_orc     = PumpCstEff()
evaporator   = HexCstEff()
expander_orc = ExpanderCstEff()
condenser    = HexCstEff()

pump_orc.set_parameters(eta_is=eta_pump_orc)
expander_orc.set_parameters(eta_is=eta_exp_orc)
evaporator.set_parameters(eta=eta_evap_orc)
condenser.set_parameters(eta=eta_cond_orc)

orc.add_component(pump_orc,     "Pump")
orc.add_component(evaporator,   "Evaporator")
orc.add_component(expander_orc, "Expander_ORC")
orc.add_component(condenser,    "Condenser")

# ORC loop links
orc.link_components("Pump",         "m-ex",   "Evaporator",   "m-su_C")
orc.link_components("Evaporator",   "m-ex_C", "Expander_ORC", "m-su")
orc.link_components("Expander_ORC", "m-ex",   "Condenser",    "m-su_H")
orc.link_components("Condenser",    "m-ex_H", "Pump",         "m-su")

# --- Hot source: Brayton exhaust air → ORC evaporator hot side ---
brayton_exhaust = MassConnector()
orc.add_source("BraytonExhaust", brayton_exhaust,
               orc.components["Evaporator"], "m-su_H")
orc.set_source_properties(target="BraytonExhaust",
                           fluid=fluid_br,
                           T=T_br_ex,
                           P=P_atm,
                           m_dot=1.0)   # per kg/s air — scaled later

# --- Cold source: ambient air → ORC condenser cold side ---
cold_source = MassConnector()
orc.add_source("ColdSource", cold_source,
               orc.components["Condenser"], "m-su_C")
orc.set_source_properties(target="ColdSource",
                           fluid=fluid_br,
                           T=T_cold_su,
                           P=P_cold,
                           m_dot=m_dot_cold)

# --- R245fa realistic pressure levels ---
# T_crit = 427 K, P_crit = 36.5 bar
# P_high ~ 20 bar → T_sat ~ 150°C — fluid superheated by Brayton exhaust (480°C)
# P_low  ~  2 bar → T_sat ~  50°C — condensation above T_amb (13°C)
from CoolProp.CoolProp import PropsSI as _PSI
T_sat_high      = _PSI('T', 'P', P_high_orc, 'Q', 1, fluid_orc)
T_sat_low       = _PSI('T', 'P', P_low_orc,  'Q', 1, fluid_orc)
T_exp_su_guess  = T_sat_high + 10   # slightly superheated at expander inlet
T_exp_ex_guess  = T_sat_low  + 20   # superheated at expander exit
T_pump_su_guess = T_sat_low  - 3    # subcooled at pump inlet (SC=3 K)

# --- Initial guesses (RecursiveCircuit: set_cycle_guess) ---
orc.set_cycle_guess(target='Pump:su',         p=P_low_orc,  SC=3,               m_dot=m_dot_orc)
orc.set_cycle_guess(target='Pump:ex',         p=P_high_orc, T=T_pump_su_guess+2)
orc.set_cycle_guess(target='Expander_ORC:su', p=P_high_orc, T=T_exp_su_guess,   m_dot=m_dot_orc)
orc.set_cycle_guess(target='Expander_ORC:ex', p=P_low_orc,  T=T_exp_ex_guess)
orc.set_cycle_guess(target='Condenser:ex_H',  p=P_low_orc,  T=T_pump_su_guess+5)

# --- Residual variables (RecursiveCircuit: target and variable are separate args) ---
orc.set_residual_variable(target='Expander_ORC:ex', variable='h', tolerance=1e-6)
orc.set_residual_variable(target='Expander_ORC:ex', variable='p', tolerance=1e-6)
orc.set_residual_variable(target='Pump:ex',         variable='h', tolerance=1e-6)
orc.set_residual_variable(target='Pump:ex',         variable='p', tolerance=1e-6)
orc.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance=1e-6)
orc.set_residual_variable(target='Evaporator:ex_C', variable='p', tolerance=1e-6)
orc.set_residual_variable(target='Condenser:ex_H',  variable='h', tolerance=1e-6)
orc.set_residual_variable(target='Condenser:ex_H',  variable='p', tolerance=1e-6)

# Solve ORC
orc.mute_print()
orc.solve(max_iter=50)

# =============================================================================
# STEP 3 — COMBINED RESULTS
# =============================================================================
print()
print("=" * 52)
print("  Step 3 — Combined Cycle Results")
print("=" * 52)

if orc.converged or all(rv.converged for rv in orc.res_vars.values()):
    W_pump_orc   = pump_orc.W.W_dot        # W per kg/s ORC fluid
    W_exp_orc    = expander_orc.W.W_dot    # W per kg/s ORC fluid
    W_net_orc    = W_exp_orc - W_pump_orc  # W per kg/s ORC fluid
    Q_evap       = evaporator.Q.Q_dot      # W per kg/s ORC fluid
    Q_cond       = condenser.Q.Q_dot       # W per kg/s ORC fluid
    eta_orc      = W_net_orc / Q_evap

    # m_dot_orc per kg/s air: Q_evap (ORC) must equal Q available from Brayton exhaust
    # Q_avail_br = m_dot_air * (h_br_ex - h_br_ex_cooled)
    # We match: m_dot_orc * Q_evap_sp = m_dot_air * Q_avail_br_sp
    # With m_dot_air = 1 kg/s (specific basis), Q_avail = W_net_br + Q_in_br - W_net_br = Q_in_br - W_net_br
    # Actually: Q_available from exhaust = m_dot_air * (h_br_ex - h_cold_ref)
    # We use the evaporator energy balance: Q_evap_sp * m_dot_orc_sp = Q_br_ex_sp
    Q_br_exhaust_sp = evaporator.Q.Q_dot   # W — already balanced by HexCstEff with m_dot=1 air
    # m_dot ratio: how much ORC fluid per kg/s air
    m_dot_orc_per_air = Q_br_exhaust_sp / Q_evap  # = 1 by construction, but explicit for clarity

    # Combined efficiency (per kg/s air basis)
    W_net_combined = W_net_br + W_net_orc * m_dot_orc_per_air
    eta_combined   = W_net_combined / Q_in_br

    print(f"  {'─'*48}")
    print(f"  BRAYTON")
    print(f"  {'─'*48}")
    print(f"  W_net_Br       = {W_net_br/1e3:7.2f}  [kW per kg/s air]")
    print(f"  Q_in_Br        = {Q_in_br/1e3:7.2f}  [kW per kg/s air]")
    print(f"  eta_Br         = {eta_br*100:7.1f}  [%]")
    print(f"  T_ex_turb_Br   = {T_br_ex - 273.15:7.1f}  [°C]")
    print(f"  {'─'*48}")
    print(f"  ORC (R245fa)")
    print(f"  {'─'*48}")
    print(f"  T_evap_su_ORC  = {evaporator.su_C.T - 273.15:7.1f}  [°C]  (pump exit)")
    print(f"  T_evap_ex_ORC  = {evaporator.ex_C.T - 273.15:7.1f}  [°C]  (expander inlet)")
    print(f"  T_exp_ex_ORC   = {expander_orc.ex.T - 273.15:7.1f}  [°C]  (condenser inlet)")
    print(f"  T_cond_ex_ORC  = {condenser.ex_H.T - 273.15:7.1f}  [°C]  (pump inlet)")
    print(f"  P_high_ORC     = {expander_orc.su.p/1e5:7.2f}  [bar]")
    print(f"  P_low_ORC      = {pump_orc.su.p/1e5:7.2f}  [bar]")
    print(f"  W_exp_ORC      = {W_exp_orc/1e3:7.2f}  [kW per kg/s ORC]")
    print(f"  W_pump_ORC     = {W_pump_orc/1e3:7.2f}  [kW per kg/s ORC]")
    print(f"  W_net_ORC      = {W_net_orc/1e3:7.2f}  [kW per kg/s ORC]")
    print(f"  Q_evap_ORC     = {Q_evap/1e3:7.2f}  [kW per kg/s ORC]")
    print(f"  eta_ORC        = {eta_orc*100:7.1f}  [%]")
    print(f"  m_dot_ORC/air  = {m_dot_orc_per_air:7.4f}  [kg/s ORC per kg/s air]")
    print(f"  {'─'*48}")
    print(f"  COMBINED")
    print(f"  {'─'*48}")
    print(f"  W_net_combined = {W_net_combined/1e3:7.2f}  [kW per kg/s air]")
    print(f"  Q_in           = {Q_in_br/1e3:7.2f}  [kW per kg/s air]")
    print(f"  eta_combined   = {eta_combined*100:7.1f}  [%]")
    print(f"  Gain vs Br     = +{(eta_combined - eta_br)*100:5.1f}  pp")
    print(f"  {'='*48}")
else:
    print("  ORC did not converge. Check initial guesses and parameters.")