# -*- coding: utf-8 -*-
"""
Combined Air Brayton + Steam Rankine cycle — scipy optimizer

Topology:
  Air -> Compressor -> Heater -> Turbine -> HEX1 (Air/Oil)
  Oil -> SH (Oil/Steam) [parallel]
      -> Evap (Oil/Steam) [parallel]
  Steam -> Expander -> Condenser -> Pump -> Evap

Free variables : [T_oil_su_C, T_oil_ret_C, P_high_bar, dT_sh]
Derived        : m_dot_oil from HEX1 balance, m_dot_wf scaled so Q_sg = Q_HEX1
Objective      : maximize eta = (W_Brayton + W_Rankine) / Q_in_Brayton
"""

#%%
import numpy as np
from scipy.optimize import minimize
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_it import IterativeCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

#%%
# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
W_net_target  = 100e6         # W — target net power output

T_amb         = 13.2 + 273.15 # K
P_atm         = 1.01325e5     # Pa
PR            = 11
eta_comp      = 0.8
eta_turb_br   = 0.9
eta_HX        = 1.0

P_salt        = 1e5           # Pa
T_salt_cold   = 290 + 273.15  # K
T_salt_limit  = 565 + 273.15  # K
T_hot_su      = 1000 + 273.15 # K

Cp_oil        = 2200.0        # J/kg/K
P_oil         = 2e5           # Pa
eta_HEX1      = 0.9
eta_evap      = 0.9
eta_sh        = 0.9
eta_pump      = 0.75
eta_turb_st   = 0.85

P_cond        = 0.10e5        # Pa — condenser pressure
PINCH_MIN     = 10.0          # K
X_MIN         = 0.88          # minimum steam quality at turbine exit

#%%
def run_brayton(print_results=False):
    if T_hot_su <= T_salt_limit:
        hot_source = SolarSaltConnector()
        hot_source.set_properties(T=T_hot_su, p=P_salt, m_dot=10.0)
        hot_fluid = 'SolarSalt'
    else:
        hot_source = MassConnector()
        hot_source.set_properties(fluid='Air', T=T_hot_su, p=P_salt, m_dot=10.0)
        hot_fluid = 'Air'

    cycle      = IterativeCircuit(fluid='Air')
    compressor = CompressorCstEff()
    heater     = HexCstEff()
    turbine    = ExpanderCstEff()

    compressor.set_parameters(eta_is=eta_comp)
    turbine.set_parameters(eta_is=eta_turb_br)
    heater.set_parameters(eta=eta_HX)

    cycle.add_component(compressor, "Compressor")
    cycle.add_component(heater,     "Heater")
    cycle.add_component(turbine,    "Turbine")

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

    W_comp = compressor.W.W_dot
    W_turb = turbine.W.W_dot
    W_net  = W_turb - W_comp
    Q_in   = heater.Q.Q_dot
    eta    = W_net / Q_in

    if T_hot_su <= T_salt_limit:
        Cp_avg    = SolarSaltConnector._Cp((T_hot_su + T_salt_cold) / 2)
        m_dot_hot = Q_in / (Cp_avg * (T_hot_su - T_salt_cold))
        hot_label = "m_dot_salt"
    else:
        h_hot_in  = PropsSI('H', 'T', T_hot_su,       'P', P_salt, 'Air')
        h_hot_ex  = PropsSI('H', 'T', heater.ex_H.T,  'P', P_salt, 'Air')
        m_dot_hot = Q_in / (h_hot_in - h_hot_ex)
        hot_label = "m_dot_air_hot"

    if print_results:
        print(f"{'─'*45}")
        print(f"  Open Air Brayton Cycle — hot source: {hot_fluid} at {T_hot_su-273.15:.0f} °C")
        print(f"{'─'*45}")
        print(f"  PR             = {PR:.1f}  [-]")
        print(f"  TIT            = {turbine.su.T - 273.15:.1f}  [°C]")
        print(f"  T_ex_comp      = {compressor.ex.T - 273.15:.1f}  [°C]")
        print(f"  T_ex_turb      = {turbine.ex.T - 273.15:.1f}  [°C]")
        print(f"  W_comp         = {W_comp/1e3:.2f}  [kW per kg/s air]")
        print(f"  W_turb         = {W_turb/1e3:.2f}  [kW per kg/s air]")
        print(f"  W_net          = {W_net/1e3:.2f}  [kW per kg/s air]")
        print(f"  Q_in           = {Q_in/1e3:.2f}  [kW per kg/s air]")
        print(f"  eta            = {eta*100:.2f}  [%]")
        print(f"  {hot_label:<15}= {m_dot_hot:.4f}  [kg/s per kg/s air]")
        print(f"{'─'*45}")

    return {
        'T_ex_turb': turbine.ex.T,
        'W_net':     W_net,
        'Q_in':      Q_in,
        'eta':       eta,
        'T_ex_comp': compressor.ex.T,
        'TIT':       turbine.su.T,
    }


#%%
def run_hex1(T_oil_ret_K, T_oil_su_K, T_ex_turb_K):
    """HEX1: Brayton exhaust air -> thermal oil. Returns (m_dot_oil, Q, T_air_ex, pinch)."""
    h_air_su = PropsSI('H', 'T', T_ex_turb_K, 'P', P_atm, 'Air')
    h_air_ex = PropsSI('H', 'T', T_oil_ret_K, 'P', P_atm, 'Air')

    Q_HEX1 = eta_HEX1 * (h_air_su - h_air_ex)
    if Q_HEX1 <= 0 or T_oil_su_K <= T_oil_ret_K:
        return None

    m_dot_oil       = Q_HEX1 / (Cp_oil * (T_oil_su_K - T_oil_ret_K))
    h_air_ex_actual = h_air_su - Q_HEX1
    T_air_ex        = PropsSI('T', 'H', h_air_ex_actual, 'P', P_atm, 'Air')

    pinch = min(T_ex_turb_K - T_oil_su_K, T_air_ex - T_oil_ret_K)
    return m_dot_oil, Q_HEX1, T_air_ex, pinch


#%%
def run_rankine(T_oil_su_K, P_high, T_sh_ex_K):
    """
    Analytical Rankine cycle (ref m_dot_wf = 1 kg/s).
    SH and Evap both fed by oil at T_oil_su independently.
    Returns perf dict or None if infeasible.
    """
    try:
        T_sat = PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water')
    except Exception:
        return None

    if T_sat  >= T_oil_su_K - PINCH_MIN: return None
    if T_sh_ex_K <= T_sat:               return None
    if T_sh_ex_K >= T_oil_su_K - PINCH_MIN: return None

    try:
        # pump
        h_pump_su    = PropsSI('H', 'P', P_cond,  'Q', 0, 'Water')
        s_pump_su    = PropsSI('S', 'P', P_cond,  'Q', 0, 'Water')
        h_pump_ex_is = PropsSI('H', 'P', P_high, 'S', s_pump_su, 'Water')
        h_pump_ex    = h_pump_su + (h_pump_ex_is - h_pump_su) / eta_pump
        W_pump       = h_pump_ex - h_pump_su

        # evaporator
        h_sat_liq = PropsSI('H', 'P', P_high, 'Q', 0, 'Water')
        h_sat_vap = PropsSI('H', 'P', P_high, 'Q', 1, 'Water')
        Q_evap    = eta_evap * (h_sat_vap - h_pump_ex)
        h_evap_ex = h_pump_ex + Q_evap
        if h_evap_ex < h_sat_liq:
            return None

        # superheater
        h_sh_ex_target = PropsSI('H', 'T', T_sh_ex_K, 'P', P_high, 'Water')
        Q_sh           = eta_sh * (h_sh_ex_target - h_evap_ex)
        if Q_sh <= 0:
            return None
        h_sh_ex        = h_evap_ex + Q_sh
        T_sh_ex_actual = PropsSI('T', 'H', h_sh_ex, 'P', P_high, 'Water')

        Q_sg = Q_evap + Q_sh
        if Q_sg <= 0:
            return None

        # expander
        s_exp_su    = PropsSI('S', 'H', h_sh_ex, 'P', P_high, 'Water')
        h_exp_ex_is = PropsSI('H', 'P', P_cond, 'S', s_exp_su, 'Water')
        h_exp_ex    = h_sh_ex - eta_turb_st * (h_sh_ex - h_exp_ex_is)
        W_exp       = h_sh_ex - h_exp_ex
        W_net       = W_exp - W_pump
        if W_net <= 0:
            return None

        T_exp_ex = PropsSI('T', 'H', h_exp_ex, 'P', P_cond, 'Water')

    except Exception:
        return None

    return {
        'W_net_sp': W_net,
        'Q_sg_sp':  Q_sg,
        'eta':      W_net / Q_sg,
        'T_exp_su': T_sh_ex_actual,
        'T_exp_ex': T_exp_ex,
        'h_exp_ex': h_exp_ex,
        'dT_sh':    T_sh_ex_actual - T_sat,
        'T_sat':    T_sat,
    }


#%%
def steam_quality(h_exp_ex):
    h_liq = PropsSI('H', 'P', P_cond, 'Q', 0, 'Water')
    h_vap = PropsSI('H', 'P', P_cond, 'Q', 1, 'Water')
    if h_exp_ex >= h_vap:
        return 1.0
    return (h_exp_ex - h_liq) / (h_vap - h_liq)


#%%
def make_bounds(T_ex_turb_K):
    return [
        (150.0, T_ex_turb_K - 273.15 - PINCH_MIN),  # T_oil_su_C
        (40.0,  200.0),                               # T_oil_ret_C
        (2.0,   150.0),                               # P_high_bar
        (0.0,   100.0),                               # dT_sh [K]
    ]


def evaluate(x, brayton, bounds):
    T_oil_su_C, T_oil_ret_C, P_high_bar, dT_sh = x

    if not (bounds[0][0] <= T_oil_su_C  <= bounds[0][1]): return 0.0, None
    if not (bounds[1][0] <= T_oil_ret_C <= bounds[1][1]): return 0.0, None
    if not (bounds[2][0] <= P_high_bar  <= bounds[2][1]): return 0.0, None
    if not (bounds[3][0] <= dT_sh       <= bounds[3][1]): return 0.0, None
    if T_oil_ret_C >= T_oil_su_C: return 0.0, None

    T_oil_su_K  = T_oil_su_C  + 273.15
    T_oil_ret_K = T_oil_ret_C + 273.15
    P_high      = P_high_bar * 1e5

    try:
        T_sat = PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water')
    except Exception:
        return 0.0, None

    T_sh_ex_K = T_sat + dT_sh
    if T_sh_ex_K >= T_oil_su_K - PINCH_MIN:
        return 0.0, None

    r_hex1 = run_hex1(T_oil_ret_K, T_oil_su_K, brayton['T_ex_turb'])
    if r_hex1 is None: return 0.0, None
    m_dot_oil, Q_HEX1, T_air_ex, pinch_hex1 = r_hex1
    if pinch_hex1 < PINCH_MIN: return 0.0, None

    r_st = run_rankine(T_oil_su_K, P_high, T_sh_ex_K)
    if r_st is None: return 0.0, None

    if steam_quality(r_st['h_exp_ex']) < X_MIN:
        return 0.0, None

    m_dot_wf  = Q_HEX1 / r_st['Q_sg_sp']
    W_net_st  = m_dot_wf * r_st['W_net_sp']
    eta       = (brayton['W_net'] + W_net_st) / brayton['Q_in']

    details = {
        'm_dot_oil':  m_dot_oil,
        'm_dot_wf':   m_dot_wf,
        'Q_HEX1':     Q_HEX1,
        'T_air_ex':   T_air_ex,
        'pinch_hex1': pinch_hex1,
        'W_net_st':   W_net_st,
        'eta_st':     r_st['eta'],
        'eta':        eta,
        'T_sat':      T_sat,
        'T_sh_ex':    r_st['T_exp_su'],
        'T_exp_su':   r_st['T_exp_su'],
        'T_exp_ex':   r_st['T_exp_ex'],
        'x_exp_ex':   steam_quality(r_st['h_exp_ex']),
        'dT_sh':      r_st['dT_sh'],
    }
    return eta, details


def objective(x, brayton, bounds):
    eta, _ = evaluate(x, brayton, bounds)
    return -eta


#%%
def diagnose_feasibility(brayton, bounds):
    print("\n  === Feasibility diagnostic ===")
    for T_oil in [200.0, 250.0, 300.0]:
        for T_ret in [50.0, 80.0]:
            for P in [2.0, 5.0, 10.0]:
                for dT in [10.0, 30.0]:
                    T_oil_su_K  = T_oil + 273.15
                    T_oil_ret_K = T_ret + 273.15
                    P_high      = P * 1e5

                    try:
                        T_sat = PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water')
                    except Exception:
                        print(f"  T_oil={T_oil:.0f} T_ret={T_ret:.0f} P={P:.0f} dT={dT:.0f} → FAIL: PropsSI")
                        continue

                    T_sh_ex_K = T_sat + dT
                    if T_sh_ex_K >= T_oil_su_K - PINCH_MIN:
                        print(f"  T_oil={T_oil:.0f} T_ret={T_ret:.0f} P={P:.0f} dT={dT:.0f} "
                              f"→ FAIL: T_sh_ex={T_sh_ex_K-273.15:.1f} >= T_oil_su-pinch={T_oil-PINCH_MIN:.1f}")
                        continue

                    r_hex1 = run_hex1(T_oil_ret_K, T_oil_su_K, brayton['T_ex_turb'])
                    if r_hex1 is None:
                        print(f"  T_oil={T_oil:.0f} T_ret={T_ret:.0f} P={P:.0f} dT={dT:.0f} → FAIL: HEX1")
                        continue
                    m_dot_oil, Q_HEX1, T_air_ex, pinch = r_hex1
                    if pinch < PINCH_MIN:
                        print(f"  T_oil={T_oil:.0f} T_ret={T_ret:.0f} P={P:.0f} dT={dT:.0f} "
                              f"→ FAIL: pinch={pinch:.1f} K")
                        continue

                    r_st = run_rankine(T_oil_su_K, P_high, T_sh_ex_K)
                    if r_st is None:
                        print(f"  T_oil={T_oil:.0f} T_ret={T_ret:.0f} P={P:.0f} dT={dT:.0f} → FAIL: Rankine")
                        continue

                    x_ex   = steam_quality(r_st['h_exp_ex'])
                    status = 'OK' if x_ex >= X_MIN else f'FAIL: x={x_ex:.3f}<{X_MIN}'
                    print(f"  T_oil={T_oil:.0f} T_ret={T_ret:.0f} P={P:.0f} dT={dT:.0f} → "
                          f"x_ex={x_ex:.3f}  W_net={r_st['W_net_sp']/1e3:.1f} kW  "
                          f"eta={r_st['eta']*100:.1f}%  {status}")
    print()


def find_x0(brayton, bounds):
    best_eta, best_x = 0.0, None
    for T_oil in [200.0, 250.0, min(300.0, bounds[0][1] - 5), bounds[0][1] - 5]:
        for T_ret in [50.0, 80.0, 120.0]:
            for P in [2.0, 3.0, 5.0, 10.0, 20.0, 40.0]:
                for dT in [5.0, 10.0, 20.0, 30.0, 50.0]:
                    eta, det = evaluate([T_oil, T_ret, P, dT], brayton, bounds)
                    if det is not None and eta > best_eta:
                        best_eta, best_x = eta, [T_oil, T_ret, P, dT]
    return best_x, best_eta


#%%
if __name__ == "__main__":

    print("=" * 50)
    print("  Step 1 — Brayton cycle")
    print("=" * 50)
    brayton = run_brayton(print_results=True)

    print("\n" + "=" * 50)
    print("  Step 2 — Combined cycle optimizer")
    print("=" * 50)

    bounds = make_bounds(brayton['T_ex_turb'])
    print(f"  T_ex_turb    = {brayton['T_ex_turb'] - 273.15:.1f} °C")
    print(f"  Bounds : T_oil_su  ∈ [{bounds[0][0]:.0f}, {bounds[0][1]:.1f}] °C")
    print(f"           T_oil_ret ∈ [{bounds[1][0]:.0f}, {bounds[1][1]:.0f}] °C")
    print(f"           P_high    ∈ [{bounds[2][0]:.0f}, {bounds[2][1]:.0f}] bar")
    print(f"           dT_sh     ∈ [{bounds[3][0]:.0f}, {bounds[3][1]:.0f}] K")
    print()

    diagnose_feasibility(brayton, bounds)

    x0, eta_x0 = find_x0(brayton, bounds)
    if x0 is None:
        print("ERROR: no feasible starting point found.")
        raise SystemExit
    print(f"  Best x0 : T_oil_su={x0[0]:.0f}°C  T_ret={x0[1]:.0f}°C  "
          f"P={x0[2]:.0f}bar  dT_sh={x0[3]:.0f}K  → eta={eta_x0*100:.2f}%")
    print()

    result = minimize(
        objective, x0,
        args=(brayton, bounds),
        method='Nelder-Mead',
        options={'xatol': 0.1, 'fatol': 1e-5, 'maxiter': 2000, 'disp': True},
    )

    x_opt = result.x
    eta_opt, det = evaluate(x_opt, brayton, bounds)
    if det is None:
        print("Warning: optimizer returned infeasible point — using best x0.")
        eta_opt, det = evaluate(x0, brayton, bounds)
        x_opt = x0

    T_sat_opt = det['T_sat'] - 273.15
    W_net_sp  = brayton['W_net']/1e3 + det['W_net_st']/1e3   # kW per kg/s air
    m_dot_air = (W_net_target / 1e3) / W_net_sp               # kg/s

    print(f"\n{'═'*52}")
    print(f"  Combined Brayton + Rankine — optimal point")
    print(f"{'═'*52}")
    print(f"  PR                = {PR:.1f}  [-]")
    print(f"  T_ex_turb         = {brayton['T_ex_turb'] - 273.15:.1f}  °C")
    print(f"  T_oil_su          = {x_opt[0]:.1f}  °C")
    print(f"  T_oil_ret         = {x_opt[1]:.1f}  °C")
    print(f"  P_high            = {x_opt[2]:.1f}  bar")
    print(f"  T_sat             = {T_sat_opt:.1f}  °C")
    print(f"  dT_sh             = {x_opt[3]:.1f}  K")
    print(f"  TIT_ST            = {det['T_sh_ex'] - 273.15:.1f}  °C")
    print(f"  T_air_ex_HEX1     = {det['T_air_ex'] - 273.15:.1f}  °C")
    print(f"  T_turb_su         = {det['T_exp_su'] - 273.15:.1f}  °C")
    print(f"  T_turb_ex         = {det['T_exp_ex'] - 273.15:.1f}  °C")
    print(f"  x_turb_ex         = {det['x_exp_ex']:.4f}  [-]  (min={X_MIN})")
    print(f"  Pinch HEX1        = {det['pinch_hex1']:.1f}  K")
    print(f"  Pinch SH          = {x_opt[0] - (det['T_sh_ex'] - 273.15):.1f}  K")
    print(f"  Pinch Evap        = {x_opt[0] - T_sat_opt:.1f}  K")
    print(f"{'─'*52}")
    print(f"  W_net_Brayton  = {brayton['W_net']/1e3:7.2f} kW/kg/s  →  {brayton['W_net']/1e6 * m_dot_air:6.2f} MW")
    print(f"  W_net_Rankine  = {det['W_net_st']/1e3:7.2f} kW/kg/s  →  {det['W_net_st']/1e6 * m_dot_air:6.2f} MW")
    print(f"  W_net_total    = {W_net_sp:7.2f} kW/kg/s  →  {W_net_sp/1e3 * m_dot_air:6.2f} MW")
    print(f"  Q_HEX1         = {det['Q_HEX1']/1e3:7.2f} kW/kg/s  →  {det['Q_HEX1']/1e6 * m_dot_air:6.2f} MW")
    print(f"  Q_in_Brayton   = {brayton['Q_in']/1e3:7.2f} kW/kg/s  →  {brayton['Q_in']/1e6 * m_dot_air:6.2f} MW")
    print(f"  m_dot_air      = {m_dot_air:.2f}  kg/s")
    print(f"  m_dot_oil      = {det['m_dot_oil']:.3f} kg/s/kg/s  →  {det['m_dot_oil'] * m_dot_air:.2f}  kg/s")
    print(f"  m_dot_wf       = {det['m_dot_wf']:.3f} kg/s/kg/s  →  {det['m_dot_wf'] * m_dot_air:.2f}  kg/s")
    print(f"{'─'*52}")
    print(f"  eta_Rankine    = {det['eta_st']*100:.2f}  %")
    print(f"  eta_Brayton    = {brayton['eta']*100:.2f}  %")
    print(f"  eta_combined   = {eta_opt*100:.2f}  %")
    print(f"  Gain           = +{(eta_opt - brayton['eta'])*100:.2f}  pp")
    print(f"{'═'*52}")