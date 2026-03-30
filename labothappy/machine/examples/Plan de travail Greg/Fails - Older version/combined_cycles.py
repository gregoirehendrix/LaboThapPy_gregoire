# -*- coding: utf-8 -*-
"""
Combined Air Brayton + Steam Rankine — scipy optimizer
Topology : Air (inlet) -> [Compressor] -> [Heater] -> [Turbine] -> [HEX1: Air/Oil] -> Oil -> [SG: Oil/Water] -> Rankine

Free variables : [T_oil_su_C, T_oil_retour_C, P_high_bar]
m_dot_oil  : derived from HEX1 energy balance
m_dot_wf   : scaled so that Q_sg = Q_HEX1 (energy consistency)
Objective  : maximize eta_combined = (W_Brayton + W_Rankine_scaled) / Q_in_Brayton
"""

#%%
import numpy as np
from scipy.optimize import minimize
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_it import IterativeCircuit
from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.pump.pump_csteff import PumpCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

#%%
# ---------------------------------------------------------------------------
# Fixed parameters
# ---------------------------------------------------------------------------
fluid        = 'Air'
T_amb        = 13.2 + 273.15   # K
P_atm        = 1.01325e5       # Pa
PR           = 11
eta_comp     = 0.8
eta_turb_br  = 0.9
eta_HX       = 1.0
P_salt       = 1e5             # Pa
T_salt_cold  = 290 + 273.15   # K
T_salt_limit = 565 + 273.15   # K
T_hot_su     = 850 + 273.15   # K

Cp_oil   = 2200.0      # J/kg/K
P_air_ex = P_atm       # Pa — turbine exhaust at atmospheric pressure
P_oil    = 2e5         # Pa
eta_HEX1 = 0.9
eta_sg   = 0.9
eta_cd   = 0.9
eta_pump = 0.75
eta_turb = 0.85

T_cw_su  = 20 + 273.15   # K
P_cw     = 2e5            # Pa
m_dot_cw = 5.0            # kg/s
P_low    = 0.10e5         # Pa
PINCH_MIN = 10.0          # K

#%%
def run_brayton(PR, T_hot_su, eta_comp, eta_turb_br, eta_HX,
                T_amb, P_atm, P_salt, T_salt_cold, T_salt_limit,
                print_results=False):
    """
    Run open Air Brayton cycle.
    Returns dict with Brayton outputs, or None if failed.
    """
    if T_hot_su <= T_salt_limit:
        hot_source = SolarSaltConnector()
        hot_source.set_properties(T=T_hot_su, p=P_salt, m_dot=10.0)
        hot_fluid = 'SolarSalt'
    else:
        hot_source = MassConnector()
        hot_source.set_properties(fluid='Air', T=T_hot_su, p=P_salt, m_dot=10.0)
        hot_fluid = 'Air'

    cycle = IterativeCircuit(fluid='Air')

    Compressor = CompressorCstEff()
    Heater     = HexCstEff()
    Turbine    = ExpanderCstEff()

    Compressor.set_parameters(eta_is=eta_comp)
    Turbine.set_parameters(eta_is=eta_turb_br)
    Heater.set_parameters(eta=eta_HX)

    cycle.add_component(Compressor, "Compressor")
    cycle.add_component(Heater,     "Heater")
    cycle.add_component(Turbine,    "Turbine")

    cycle.link_components("Compressor", "m-ex",   "Heater",  "m-su_C")
    cycle.link_components("Heater",     "m-ex_C", "Turbine", "m-su")

    air_in = MassConnector()
    cycle.add_source("AirInlet", air_in, cycle.components["Compressor"], "m-su")
    cycle.set_source_properties(target="AirInlet", fluid='Air', T=T_amb, P=P_atm, m_dot=1.0)

    cycle.add_source("HotSource", hot_source, cycle.components["Heater"], "m-su_H")
    cycle.set_source_properties(target="HotSource", fluid=hot_fluid, T=T_hot_su, P=P_salt, m_dot=10.0)

    cycle.set_cycle_input(target="Compressor:ex", p=P_atm * PR)
    cycle.set_cycle_input(target="Turbine:ex",    p=P_atm)

    cycle._build_solve_order()
    for name in cycle.solve_start_components:
        cycle.components[name].solve()

    W_comp = Compressor.W.W_dot
    W_turb = Turbine.W.W_dot
    W_net  = W_turb - W_comp
    Q_in   = Heater.Q.Q_dot
    eta_th = W_net / Q_in

    if T_hot_su <= T_salt_limit:
        Cp_salt_avg = SolarSaltConnector._Cp((T_hot_su + T_salt_cold) / 2)
        m_dot_hot   = Q_in / (Cp_salt_avg * (T_hot_su - T_salt_cold))
        hot_label   = "m_dot_salt"
    else:
        T_hot_ex  = Heater.ex_H.T
        h_hot_in  = PropsSI('H', 'T', T_hot_su, 'P', P_salt, 'Air')
        h_hot_ex  = PropsSI('H', 'T', T_hot_ex, 'P', P_salt, 'Air')
        m_dot_hot = Q_in / (h_hot_in - h_hot_ex)
        hot_label = "m_dot_air_hot"

    if print_results:
        print(f"{'─'*45}")
        print(f"  Simple Open Air Brayton Cycle")
        print(f"  Hot source : {hot_fluid} at {T_hot_su - 273.15:.1f} °C")
        print(f"{'─'*45}")
        print(f"  PR             = {PR:.1f}  [-]")
        print(f"  TIT            = {Turbine.su.T - 273.15:.1f}  [°C]")
        print(f"  T_ex_comp      = {Compressor.ex.T - 273.15:.1f}  [°C]")
        print(f"  T_ex_turb      = {Turbine.ex.T - 273.15:.1f}  [°C]")
        print(f"  dT_turb        = {Turbine.su.T - Turbine.ex.T:.1f}  [°C]")
        print(f"  T_ex_heater_H  = {Heater.ex_H.T - 273.15:.1f}  [°C]")
        print(f"  W_comp         = {W_comp/1e3:.2f}  [kW per kg/s air]")
        print(f"  W_turb         = {W_turb/1e3:.2f}  [kW per kg/s air]")
        print(f"  W_net          = {W_net/1e3:.2f}  [kW per kg/s air]")
        print(f"  Q_in           = {Q_in/1e3:.2f}  [kW per kg/s air]")
        print(f"  eta_th         = {eta_th*100:.1f}  [%]")
        print(f"  {hot_label:<15}= {m_dot_hot:.4f}  [kg/s per kg/s air]")
        print(f"{'─'*45}")

    return {
        'T_ex_turb':    Turbine.ex.T,
        'W_net':        W_net,
        'Q_in':         Q_in,
        'eta':          eta_th,
        'T_ex_comp':    Compressor.ex.T,
        'TIT':          Turbine.su.T,
        'm_dot_hot':    m_dot_hot,
    }

#%%
def make_oil_connector(T_K, p_Pa, m_dot_kgs, Cp_Jkg):
    oil = MassConnector()
    oil.fluid        = 'OilHTF'
    oil.T            = T_K
    oil.p            = p_Pa
    oil.m_dot        = m_dot_kgs
    oil.h            = Cp_Jkg * T_K
    oil.cp           = Cp_Jkg
    oil.state_known  = True
    oil.completely_known = True
    return oil

#%%
def run_hex1(T_oil_retour_K, T_oil_su_K, T_ex_turb_K):
    """
    HEX1: Air (hot, from Brayton turbine exhaust) -> Oil (cold), counterflow.
    Returns (m_dot_oil, Q_HEX1, T_air_ex, pinch) or None if infeasible.
    """
    h_air_su = PropsSI('H', 'T', T_ex_turb_K,   'P', P_air_ex, 'Air')
    h_air_ex = PropsSI('H', 'T', T_oil_retour_K, 'P', P_air_ex, 'Air')

    Q_air_max = h_air_su - h_air_ex   # per kg/s air
    Q_HEX1    = eta_HEX1 * Q_air_max

    if Q_HEX1 <= 0:
        return None

    dT_oil = T_oil_su_K - T_oil_retour_K
    if dT_oil <= 0:
        return None

    m_dot_oil = Q_HEX1 / (Cp_oil * dT_oil)

    h_air_ex_actual = h_air_su - Q_HEX1
    T_air_ex        = PropsSI('T', 'H', h_air_ex_actual, 'P', P_air_ex, 'Air')

    pinch_hot_end  = T_ex_turb_K   - T_oil_su_K
    pinch_cold_end = T_air_ex      - T_oil_retour_K
    pinch          = min(pinch_hot_end, pinch_cold_end)

    return m_dot_oil, Q_HEX1, T_air_ex, pinch

#%%
def run_rankine_unit(T_oil_su_K, m_dot_oil, P_high):
    """
    Run Rankine at m_dot_wf = 1 kg/s reference.
    Returns perf dict (specific quantities per kg/s wf) or None if infeasible.
    """
    T_sat_high = PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water')
    if T_sat_high >= T_oil_su_K - PINCH_MIN:
        return None

    T_pump_su_guess = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_turb_su_guess = PropsSI('T', 'P', P_high, 'Q', 1, 'Water') + 20

    HSource = make_oil_connector(T_oil_su_K, P_oil, m_dot_oil, Cp_oil)
    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=T_cw_su, p=P_cw, m_dot=m_dot_cw)

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

    cycle.link_components("Pump",     "m-ex",   "SG",       "m-su_C")
    cycle.link_components("SG",       "m-ex_C", "Expander", "m-su")
    cycle.link_components("Expander", "m-ex",   "CD",       "m-su_H")
    cycle.link_components("CD",       "m-ex_H", "Pump",     "m-su")

    cycle.add_source("Hot_source",  HSource, cycle.components["SG"], "m-su_H")
    cycle.add_source("Cold_source", CSource, cycle.components["CD"], "m-su_C")

    cycle.set_cycle_guess(target='Pump:su',     m_dot=1.0, T=T_pump_su_guess, p=P_low)
    cycle.set_fixed_properties(target='Pump:ex',     p=P_high)
    cycle.set_cycle_guess(target='Expander:su', m_dot=1.0, T=T_turb_su_guess, p=P_high)
    cycle.set_fixed_properties(target='Expander:ex', p=P_low)

    cycle.mute_print()

    try:
        cycle.solve()
    except Exception:
        return None

    if not all(cycle.components[c].model.solved for c in cycle.components):
        return None

    W_pump = pump.ex.h     - pump.su.h
    W_exp  = expander.su.h - expander.ex.h
    W_net  = W_exp - W_pump
    Q_sg   = sg.ex_C.h    - sg.su_C.h

    if Q_sg <= 0 or W_net <= 0:
        return None

    return {
        'W_net_sp':  W_net,
        'Q_sg_sp':   Q_sg,
        'eta':       W_net / Q_sg,
        'T_exp_su':  expander.su.T,
        'T_exp_ex':  expander.ex.T,
    }

#%%
def make_bounds(T_ex_turb_K):
    return [
        (150.0, T_ex_turb_K - 273.15 - PINCH_MIN),   # T_oil_su_C
        (40.0,  200.0),                                # T_oil_retour_C
        (2.0,   150.0),                                # P_high_bar
    ]

def evaluate(x, brayton, bounds):
    T_oil_su_C, T_oil_retour_C, P_high_bar = x[0], x[1], x[2]

    if not (bounds[0][0] <= T_oil_su_C     <= bounds[0][1]): return 0.0, None
    if not (bounds[1][0] <= T_oil_retour_C <= bounds[1][1]): return 0.0, None
    if not (bounds[2][0] <= P_high_bar     <= bounds[2][1]): return 0.0, None
    if T_oil_retour_C >= T_oil_su_C: return 0.0, None

    T_oil_su_K     = T_oil_su_C     + 273.15
    T_oil_retour_K = T_oil_retour_C + 273.15
    P_high         = P_high_bar * 1e5
    T_ex_turb_K    = brayton['T_ex_turb']

    r_hex1 = run_hex1(T_oil_retour_K, T_oil_su_K, T_ex_turb_K)
    if r_hex1 is None: return 0.0, None
    m_dot_oil, Q_HEX1, T_air_ex, pinch_hex1 = r_hex1
    if pinch_hex1 < PINCH_MIN: return 0.0, None

    r_rankine = run_rankine_unit(T_oil_su_K, m_dot_oil, P_high)
    if r_rankine is None: return 0.0, None

    m_dot_wf_scaled = Q_HEX1 / r_rankine['Q_sg_sp']
    W_net_rankine   = m_dot_wf_scaled * r_rankine['W_net_sp']
    eta_combined    = (brayton['W_net'] + W_net_rankine) / brayton['Q_in']

    details = {
        'm_dot_oil':     m_dot_oil,
        'm_dot_wf':      m_dot_wf_scaled,
        'Q_HEX1':        Q_HEX1,
        'T_air_ex':      T_air_ex,
        'pinch_hex1':    pinch_hex1,
        'W_net_rankine': W_net_rankine,
        'eta_rankine':   r_rankine['eta'],
        'eta_combined':  eta_combined,
        'T_sat':         PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water'),
        'T_exp_su':      r_rankine['T_exp_su'],
        'T_exp_ex':      r_rankine['T_exp_ex'],
    }
    return eta_combined, details

def objective(x, brayton, bounds):
    eta, _ = evaluate(x, brayton, bounds)
    return -eta

#%%
if __name__ == "__main__":

    # --- Step 1: Run Brayton ---
    print("=" * 45)
    print("  Step 1 — Brayton cycle")
    print("=" * 45)
    brayton = run_brayton(
        PR, T_hot_su, eta_comp, eta_turb_br, eta_HX,
        T_amb, P_atm, P_salt, T_salt_cold, T_salt_limit,
        print_results=True
    )

    # --- Step 2: Optimize Rankine ---
    print("\n" + "=" * 45)
    print("  Step 2 — Rankine optimizer")
    print("=" * 45)

    bounds = make_bounds(brayton['T_ex_turb'])
    x0     = [min(250.0, bounds[0][1] - 5), 80.0, 20.0]

    print(f"  T_ex_turb    = {brayton['T_ex_turb'] - 273.15:.1f} °C")
    print(f"  Bounds : T_oil_su     ∈ [{bounds[0][0]:.0f}, {bounds[0][1]:.1f}] °C")
    print(f"           T_oil_retour ∈ [{bounds[1][0]:.0f}, {bounds[1][1]:.0f}] °C")
    print(f"           P_high       ∈ [{bounds[2][0]:.0f}, {bounds[2][1]:.0f}] bar")
    print()

    result = minimize(
        objective,
        x0,
        args=(brayton, bounds),
        method='Nelder-Mead',
        bounds=bounds,
        options={'xatol': 0.1, 'fatol': 1e-5, 'maxiter': 1000, 'disp': True},
    )

    x_opt = result.x
    eta_opt, details = evaluate(x_opt, brayton, bounds)

    if details is None:
        print("Warning: optimum point is infeasible — try different x0.")
    else:
        T_sat_opt = details['T_sat'] - 273.15

        print(f"\n{'═'*52}")
        print(f"  Optimum — Combined Brayton + Rankine")
        print(f"{'═'*52}")
        print(f"  PR                = {PR:.1f}  [-]")
        print(f"  T_ex_turb         = {brayton['T_ex_turb'] - 273.15:.1f}  °C")
        print(f"  T_oil_su          = {x_opt[0]:.1f}  °C")
        print(f"  T_oil_retour      = {x_opt[1]:.1f}  °C")
        print(f"  P_high            = {x_opt[2]:.1f}  bar")
        print(f"  T_sat(P_high)     = {T_sat_opt:.1f}  °C")
        print(f"  m_dot_oil         = {details['m_dot_oil']:.3f}  kg/s  (per kg/s air)")
        print(f"  m_dot_wf          = {details['m_dot_wf']:.3f}  kg/s  (per kg/s air)")
        print(f"  T_air_ex_HEX1     = {details['T_air_ex'] - 273.15:.1f}  °C")
        print(f"  Pinch HEX1        = {details['pinch_hex1']:.1f}  K")
        print(f"  Pinch SG (T_sat)  = {x_opt[0] - T_sat_opt:.1f}  K")
        print(f"  T_turb_su         = {details['T_exp_su'] - 273.15:.1f}  °C")
        print(f"  T_turb_ex         = {details['T_exp_ex'] - 273.15:.1f}  °C")
        print(f"{'─'*52}")
        print(f"  W_net_Brayton     = {brayton['W_net']/1e3:.2f}  kW  (per kg/s air)")
        print(f"  W_net_Rankine     = {details['W_net_rankine']/1e3:.2f}  kW  (per kg/s air)")
        print(f"  Q_HEX1 = Q_sg     = {details['Q_HEX1']/1e3:.2f}  kW  (energy consistent ✓)")
        print(f"  eta_Rankine       = {details['eta_rankine']*100:.2f}  %")
        print(f"  eta_Brayton       = {brayton['eta']*100:.2f}  %")
        print(f"  eta_combined      = {eta_opt*100:.2f}  %")
        print(f"  Gain vs Brayton   = +{(eta_opt - brayton['eta'])*100:.2f}  pp")
        print(f"{'═'*52}")