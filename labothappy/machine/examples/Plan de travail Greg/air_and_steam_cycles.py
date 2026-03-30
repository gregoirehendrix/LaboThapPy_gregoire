# -*- coding: utf-8 -*-
"""
Combined Air Brayton + Steam Rankine cycle — scipy optimizer

Topology:
  Air -> Compressor -> Heater -> Turbine -> HEX1 (Air/Oil)
  Oil -> [pipe losses] -> SH  (Oil/Steam) [parallel]
                       -> Evap (Oil/Steam) [parallel]
                       -> RH  (Oil/Steam) [parallel, optional]
  Steam -> Expander HP -> [RH] -> Expander LP -> Condenser -> Pump -> Evap

Free variables (no reheat) : [T_oil_su_C, T_oil_ret_C, P_high_bar, dT_sh]
Free variables (reheat)    : [T_oil_su_C, T_oil_ret_C, P_high_bar, dT_sh,
                               P_rh_bar, dT_rh]
Derived : m_dot_oil from HEX1 balance, m_dot_wf scaled so Q_sg = Q_HEX1
Objective : maximize eta = (W_Brayton + W_Rankine) / Q_in_Brayton

Pipe losses integrated in optimisation:
  T_oil_sg_in = T_amb + (T_oil_su - T_amb) * exp(-U*pi*D*L / (m_dot_oil * Cp_oil))
  Rankine and all downstream pinch checks use T_oil_sg_in, not T_oil_su.
"""

# ===========================================================================
# USER OPTIONS
# ===========================================================================
USE_REHEAT   = True     # True = Rankine with reheat, False = simple Rankine
W_net_target = 100e6   # W — target net power output

PIPE_LOSSES  = True     # True = include oil pipe thermal losses in optimisation
D_pipe       = 0.14     # m       — pipe outer diameter
L_pipe       = 350    # m       — root min square to be more precise
U_loss       = 1.5      # W/m²/K — overall heat loss coefficient
# ===========================================================================

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
T_amb         = 13.2 + 273.15  # K
P_atm         = 1.01325e5      # Pa
PR            = 11
eta_comp      = 0.8
eta_turb_br   = 0.9
eta_HX        = 1.0

P_salt        = 1e5            # Pa
T_salt_cold   = 290 + 273.15   # K
T_salt_limit  = 565 + 273.15   # K
T_hot_su      = 850 + 273.15   # K

Cp_oil        = 2200.0         # J/kg/K
P_oil         = 2e5            # Pa
eta_HEX1      = 0.9
eta_evap      = 0.9
eta_sh        = 0.9
eta_rh        = 0.9
eta_pump      = 0.75
eta_turb_st   = 0.85

P_cond        = 0.10e5         # Pa
PINCH_MIN     = 10.0           # K
X_MIN         = 0.88           # minimum steam quality at LP turbine exit

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
    cycle.set_source_properties(target="HotSource", fluid=hot_fluid,
                                T=T_hot_su, P=P_salt, m_dot=10.0)

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
        h_hot_in  = PropsSI('H', 'T', T_hot_su,      'P', P_salt, 'Air')
        h_hot_ex  = PropsSI('H', 'T', heater.ex_H.T, 'P', P_salt, 'Air')
        m_dot_hot = Q_in / (h_hot_in - h_hot_ex)
        hot_label = "m_dot_air_hot"

    if print_results:
        print(f"{'─'*48}")
        print(f"  Open Air Brayton — hot source: {hot_fluid} at {T_hot_su-273.15:.0f} °C")
        print(f"{'─'*48}")
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
        print(f"{'─'*48}")

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
    h_air_su = PropsSI('H', 'T', T_ex_turb_K, 'P', P_atm, 'Air')
    h_air_ex = PropsSI('H', 'T', T_oil_ret_K, 'P', P_atm, 'Air')

    Q_HEX1 = eta_HEX1 * (h_air_su - h_air_ex)
    if Q_HEX1 <= 0 or T_oil_su_K <= T_oil_ret_K:
        return None

    m_dot_oil       = Q_HEX1 / (Cp_oil * (T_oil_su_K - T_oil_ret_K))
    h_air_ex_actual = h_air_su - Q_HEX1
    T_air_ex        = PropsSI('T', 'H', h_air_ex_actual, 'P', P_atm, 'Air')
    pinch           = min(T_ex_turb_K - T_oil_su_K, T_air_ex - T_oil_ret_K)
    return m_dot_oil, Q_HEX1, T_air_ex, pinch


#%%
def oil_pipe_losses(T_oil_in_K, m_dot_oil_sp):
    """
    Exact exponential model for pipe thermal losses (constant Cp, uniform T_amb).
    Operates on specific flow rate [kg/s per kg/s air] — consistent with the
    rest of the specific model. NTU is correctly defined per unit air flow.
    Returns T_out [K], Q_loss [W/kg/s_air], dT [K], A [m²], NTU [-].
    """
    A_pipe  = np.pi * D_pipe * L_pipe
    NTU     = U_loss * A_pipe / (m_dot_oil_sp * Cp_oil)
    T_out_K = T_amb + (T_oil_in_K - T_amb) * np.exp(-NTU)
    Q_loss  = m_dot_oil_sp * Cp_oil * (T_oil_in_K - T_out_K)
    dT      = T_oil_in_K - T_out_K
    return T_out_K, Q_loss, dT, A_pipe, NTU


#%%
def run_rankine(T_oil_su_K, P_high, T_sh_ex_K):
    """Simple Rankine without reheat (ref m_dot_wf = 1 kg/s)."""
    try:
        T_sat = PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water')
    except Exception:
        return None

    if T_sat     >= T_oil_su_K - PINCH_MIN: return None
    if T_sh_ex_K <= T_sat:                  return None
    if T_sh_ex_K >= T_oil_su_K - PINCH_MIN: return None

    try:
        h_pump_su    = PropsSI('H', 'P', P_cond, 'Q', 0,   'Water')
        s_pump_su    = PropsSI('S', 'P', P_cond, 'Q', 0,   'Water')
        h_pump_ex_is = PropsSI('H', 'P', P_high, 'S', s_pump_su, 'Water')
        h_pump_ex    = h_pump_su + (h_pump_ex_is - h_pump_su) / eta_pump
        W_pump       = h_pump_ex - h_pump_su

        h_sat_liq = PropsSI('H', 'P', P_high, 'Q', 0, 'Water')
        h_sat_vap = PropsSI('H', 'P', P_high, 'Q', 1, 'Water')
        Q_evap    = eta_evap * (h_sat_vap - h_pump_ex)
        h_evap_ex = h_pump_ex + Q_evap
        if h_evap_ex < h_sat_liq: return None

        h_sh_target    = PropsSI('H', 'T', T_sh_ex_K, 'P', P_high, 'Water')
        Q_sh           = eta_sh * (h_sh_target - h_evap_ex)
        if Q_sh <= 0: return None
        h_sh_ex        = h_evap_ex + Q_sh
        T_sh_ex_actual = PropsSI('T', 'H', h_sh_ex, 'P', P_high, 'Water')

        Q_sg = Q_evap + Q_sh
        if Q_sg <= 0: return None

        s_exp_su    = PropsSI('S', 'H', h_sh_ex, 'P', P_high, 'Water')
        h_exp_ex_is = PropsSI('H', 'P', P_cond, 'S', s_exp_su, 'Water')
        h_exp_ex    = h_sh_ex - eta_turb_st * (h_sh_ex - h_exp_ex_is)
        W_exp       = h_sh_ex - h_exp_ex
        W_net       = W_exp - W_pump
        if W_net <= 0: return None

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
def run_rankine_rh(T_oil_su_K, P_high, T_sh_ex_K, P_rh, T_rh_ex_K):
    """Rankine with reheat (ref m_dot_wf = 1 kg/s)."""
    try:
        T_sat_hp = PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water')
        T_sat_rh = PropsSI('T', 'P', P_rh,   'Q', 0.5, 'Water')
    except Exception:
        return None

    if P_rh      >= P_high:                 return None
    if P_rh      <= P_cond:                 return None
    if T_sat_hp  >= T_oil_su_K - PINCH_MIN: return None
    if T_sh_ex_K <= T_sat_hp:               return None
    if T_sh_ex_K >= T_oil_su_K - PINCH_MIN: return None
    if T_rh_ex_K <= T_sat_rh:               return None
    if T_rh_ex_K >= T_oil_su_K - PINCH_MIN: return None

    try:
        h_pump_su    = PropsSI('H', 'P', P_cond, 'Q', 0,   'Water')
        s_pump_su    = PropsSI('S', 'P', P_cond, 'Q', 0,   'Water')
        h_pump_ex_is = PropsSI('H', 'P', P_high, 'S', s_pump_su, 'Water')
        h_pump_ex    = h_pump_su + (h_pump_ex_is - h_pump_su) / eta_pump
        W_pump       = h_pump_ex - h_pump_su

        h_sat_liq = PropsSI('H', 'P', P_high, 'Q', 0, 'Water')
        h_sat_vap = PropsSI('H', 'P', P_high, 'Q', 1, 'Water')
        Q_evap    = eta_evap * (h_sat_vap - h_pump_ex)
        h_evap_ex = h_pump_ex + Q_evap
        if h_evap_ex < h_sat_liq: return None

        h_sh_target    = PropsSI('H', 'T', T_sh_ex_K, 'P', P_high, 'Water')
        Q_sh           = eta_sh * (h_sh_target - h_evap_ex)
        if Q_sh <= 0: return None
        h_sh_ex        = h_evap_ex + Q_sh
        T_sh_ex_actual = PropsSI('T', 'H', h_sh_ex, 'P', P_high, 'Water')

        s_hp_su    = PropsSI('S', 'H', h_sh_ex, 'P', P_high, 'Water')
        h_hp_ex_is = PropsSI('H', 'P', P_rh,   'S', s_hp_su, 'Water')
        h_hp_ex    = h_sh_ex - eta_turb_st * (h_sh_ex - h_hp_ex_is)
        W_hp       = h_sh_ex - h_hp_ex

        h_rh_target    = PropsSI('H', 'T', T_rh_ex_K, 'P', P_rh, 'Water')
        Q_rh           = eta_rh * (h_rh_target - h_hp_ex)
        if Q_rh <= 0: return None
        h_rh_ex        = h_hp_ex + Q_rh
        T_rh_ex_actual = PropsSI('T', 'H', h_rh_ex, 'P', P_rh, 'Water')

        s_lp_su    = PropsSI('S', 'H', h_rh_ex, 'P', P_rh,  'Water')
        h_lp_ex_is = PropsSI('H', 'P', P_cond,  'S', s_lp_su, 'Water')
        h_lp_ex    = h_rh_ex - eta_turb_st * (h_rh_ex - h_lp_ex_is)
        W_lp       = h_rh_ex - h_lp_ex

        W_net = W_hp + W_lp - W_pump
        if W_net <= 0: return None

        Q_sg = Q_evap + Q_sh + Q_rh
        if Q_sg <= 0: return None

        T_lp_ex = PropsSI('T', 'H', h_lp_ex, 'P', P_cond, 'Water')

    except Exception:
        return None

    return {
        'W_net_sp': W_net,
        'Q_sg_sp':  Q_sg,
        'eta':      W_net / Q_sg,
        'T_sh_ex':  T_sh_ex_actual,
        'T_rh_ex':  T_rh_ex_actual,
        'T_lp_ex':  T_lp_ex,
        'h_lp_ex':  h_lp_ex,
        'dT_sh':    T_sh_ex_actual - T_sat_hp,
        'dT_rh':    T_rh_ex_actual - T_sat_rh,
        'T_sat_hp': T_sat_hp,
        'T_sat_rh': T_sat_rh,
        'W_hp':     W_hp,
        'W_lp':     W_lp,
        'Q_evap':   Q_evap,
        'Q_sh':     Q_sh,
        'Q_rh':     Q_rh,
    }


#%%
def steam_quality(h_ex):
    h_liq = PropsSI('H', 'P', P_cond, 'Q', 0, 'Water')
    h_vap = PropsSI('H', 'P', P_cond, 'Q', 1, 'Water')
    if h_ex >= h_vap:
        return 1.0
    return (h_ex - h_liq) / (h_vap - h_liq)


#%%
def make_bounds(T_ex_turb_K):
    T_oil_su_max = T_ex_turb_K - 273.15 - PINCH_MIN
    base = [
        (150.0, T_oil_su_max),
        (40.0,  200.0),
        (2.0,   150.0),
        (0.0,   100.0),
    ]
    if USE_REHEAT:
        base += [(1.0, 100.0), (0.0, 100.0)]
    return base


def evaluate(x, brayton, bounds):
    if USE_REHEAT:
        T_oil_su_C, T_oil_ret_C, P_high_bar, dT_sh, P_rh_bar, dT_rh = x
    else:
        T_oil_su_C, T_oil_ret_C, P_high_bar, dT_sh = x

    for i, xi in enumerate(x):
        if not (bounds[i][0] <= xi <= bounds[i][1]):
            return 0.0, None
    if T_oil_ret_C >= T_oil_su_C:
        return 0.0, None

    T_oil_su_K  = T_oil_su_C  + 273.15
    T_oil_ret_K = T_oil_ret_C + 273.15
    P_high      = P_high_bar * 1e5

    try:
        T_sat = PropsSI('T', 'P', P_high, 'Q', 0.5, 'Water')
    except Exception:
        return 0.0, None

    r_hex1 = run_hex1(T_oil_ret_K, T_oil_su_K, brayton['T_ex_turb'])
    if r_hex1 is None: return 0.0, None
    m_dot_oil, Q_HEX1, T_air_ex, pinch_hex1 = r_hex1
    if pinch_hex1 < PINCH_MIN: return 0.0, None

    # pipe losses between HEX1 and steam generator — integrated in optimisation
    if PIPE_LOSSES:
        T_oil_sg_K, Q_loss_sp, dT_pipe, A_pipe, NTU_pipe = oil_pipe_losses(
            T_oil_in_K   = T_oil_su_K,
            m_dot_oil_sp = m_dot_oil,
        )
        if T_oil_sg_K <= T_oil_ret_K + PINCH_MIN:
            return 0.0, None
    else:
        T_oil_sg_K = T_oil_su_K
        Q_loss_sp  = 0.0
        dT_pipe    = 0.0
        A_pipe     = 0.0
        NTU_pipe   = 0.0

    # all downstream constraints use T_oil_sg_K (temperature actually available at SG)
    T_sh_ex_K = T_sat + dT_sh
    if T_sh_ex_K >= T_oil_sg_K - PINCH_MIN:
        return 0.0, None

    if USE_REHEAT:
        P_rh = P_rh_bar * 1e5
        if P_rh >= P_high: return 0.0, None
        try:
            T_sat_rh = PropsSI('T', 'P', P_rh, 'Q', 0.5, 'Water')
        except Exception:
            return 0.0, None
        T_rh_ex_K = T_sat_rh + dT_rh
        if T_rh_ex_K >= T_oil_sg_K - PINCH_MIN: return 0.0, None

        r_st = run_rankine_rh(T_oil_sg_K, P_high, T_sh_ex_K, P_rh, T_rh_ex_K)
        if r_st is None: return 0.0, None
        x_lp_ex = steam_quality(r_st['h_lp_ex'])
        if x_lp_ex < X_MIN: return 0.0, None
    else:
        r_st = run_rankine(T_oil_sg_K, P_high, T_sh_ex_K)
        if r_st is None: return 0.0, None
        x_lp_ex = steam_quality(r_st['h_exp_ex'])
        if x_lp_ex < X_MIN: return 0.0, None

    m_dot_wf = Q_HEX1 / r_st['Q_sg_sp']
    W_net_st = m_dot_wf * r_st['W_net_sp']
    eta      = (brayton['W_net'] + W_net_st) / brayton['Q_in']

    details = {
        'm_dot_oil':  m_dot_oil,
        'm_dot_wf':   m_dot_wf,
        'Q_HEX1':     Q_HEX1,
        'T_air_ex':   T_air_ex,
        'pinch_hex1': pinch_hex1,
        'T_oil_sg':   T_oil_sg_K,
        'dT_pipe':    dT_pipe,
        'Q_loss_sp':  Q_loss_sp,
        'A_pipe':     A_pipe,
        'NTU_pipe':   NTU_pipe,
        'W_net_st':   W_net_st,
        'eta_st':     r_st['eta'],
        'eta':        eta,
        'T_sat_hp':   r_st.get('T_sat_hp', r_st.get('T_sat')),
        'T_sh_ex':    r_st.get('T_sh_ex', r_st.get('T_exp_su')),
        'T_exp_ex':   r_st.get('T_lp_ex', r_st.get('T_exp_ex')),
        'x_exp_ex':   x_lp_ex,
        'dT_sh':      r_st['dT_sh'],
        'r_st':       r_st,
    }
    if USE_REHEAT:
        details.update({
            'P_rh':     P_rh,
            'T_sat_rh': r_st['T_sat_rh'],
            'T_rh_ex':  r_st['T_rh_ex'],
            'dT_rh':    r_st['dT_rh'],
            'W_hp':     r_st['W_hp'],
            'W_lp':     r_st['W_lp'],
            'Q_rh':     r_st['Q_rh'],
        })
    return eta, details


def objective(x, brayton, bounds):
    eta, _ = evaluate(x, brayton, bounds)
    return -eta


#%%
def diagnose_feasibility(brayton, bounds):
    print("\n  === Feasibility diagnostic ===")
    P_rh_list  = [2.0, 5.0] if USE_REHEAT else [None]
    dT_rh_list = [20.0, 40.0] if USE_REHEAT else [None]

    for T_oil in [200.0, 300.0, bounds[0][1] - 5]:
        for T_ret in [50.0, 80.0]:
            for P in [5.0, 20.0, 60.0]:
                for dT in [10.0, 30.0]:
                    for P_rh in P_rh_list:
                        for dT_rh in dT_rh_list:
                            x_try = [T_oil, T_ret, P, dT, P_rh, dT_rh] if USE_REHEAT \
                                    else [T_oil, T_ret, P, dT]
                            eta, det = evaluate(x_try, brayton, bounds)
                            label = f"T_oil={T_oil:.0f} T_ret={T_ret:.0f} P={P:.0f} dT={dT:.0f}"
                            if USE_REHEAT:
                                label += f" P_rh={P_rh:.0f} dT_rh={dT_rh:.0f}"
                            if det is not None:
                                pipe_str = f"  dT_pipe={det['dT_pipe']:.2f}K  T_sg={det['T_oil_sg']-273.15:.1f}°C" \
                                           if PIPE_LOSSES else ""
                                print(f"  {label} → x_ex={det['x_exp_ex']:.3f}  "
                                      f"eta={eta*100:.2f}%{pipe_str}  OK")
                            else:
                                print(f"  {label} → FAIL")
    print()


def find_x0(brayton, bounds):
    best_eta, best_x = 0.0, None
    T_oils  = [200.0, 250.0, min(300.0, bounds[0][1]-5), bounds[0][1]-5]
    T_rets  = [50.0, 80.0, 120.0]
    P_highs = [5.0, 10.0, 20.0, 40.0, 80.0]
    dT_shs  = [10.0, 20.0, 30.0, 50.0]

    if USE_REHEAT:
        for T_oil in T_oils:
            for T_ret in T_rets:
                for P in P_highs:
                    for dT in dT_shs:
                        for P_rh in [2.0, 5.0, 10.0, 20.0]:
                            for dT_rh in [10.0, 20.0, 40.0]:
                                x_try = [T_oil, T_ret, P, dT, P_rh, dT_rh]
                                eta, det = evaluate(x_try, brayton, bounds)
                                if det is not None and eta > best_eta:
                                    best_eta, best_x = eta, x_try
    else:
        for T_oil in T_oils:
            for T_ret in T_rets:
                for P in [2.0, 3.0, 5.0, 10.0, 20.0, 40.0]:
                    for dT in [5.0, 10.0, 20.0, 30.0, 50.0]:
                        x_try = [T_oil, T_ret, P, dT]
                        eta, det = evaluate(x_try, brayton, bounds)
                        if det is not None and eta > best_eta:
                            best_eta, best_x = eta, x_try
    return best_x, best_eta


#%%
def print_results(x_opt, eta_opt, det, brayton, m_dot_air):
    T_sat_opt = det['T_sat_hp'] - 273.15
    W_net_sp  = brayton['W_net']/1e3 + det['W_net_st']/1e3
    rh_tag    = " + Reheat" if USE_REHEAT else ""

    print(f"\n{'═'*54}")
    print(f"  Combined Brayton + Rankine{rh_tag} — optimal point")
    print(f"{'═'*54}")
    print(f"  PR                = {PR:.1f}  [-]")
    print(f"  T_ex_turb         = {brayton['T_ex_turb'] - 273.15:.1f}  °C")
    print(f"  T_oil_su (HEX1)   = {x_opt[0]:.1f}  °C")
    print(f"  T_oil_ret         = {x_opt[1]:.1f}  °C")
    if PIPE_LOSSES:
        print(f"  ── Oil pipe losses ──────────────────────────")
        print(f"  D / L / U         = {D_pipe:.2f} m / {L_pipe:.0f} m / {U_loss:.2f} W/m²/K")
        print(f"  A_pipe            = {det['A_pipe']:.2f}  m²")
        print(f"  NTU_pipe          = {det['NTU_pipe']:.5f}  [-]")
        print(f"  dT_pipe           = {det['dT_pipe']:.3f}  K")
        print(f"  Q_loss            = {det['Q_loss_sp']/1e3:.4f}  kW per kg/s air")
        print(f"  T_oil_SG_in       = {det['T_oil_sg'] - 273.15:.2f}  °C")
        print(f"  ─────────────────────────────────────────────")
    print(f"  P_high            = {x_opt[2]:.1f}  bar")
    print(f"  T_sat_HP          = {T_sat_opt:.1f}  °C")
    print(f"  dT_sh             = {x_opt[3]:.1f}  K")
    print(f"  TIT_ST            = {det['T_sh_ex'] - 273.15:.1f}  °C")
    if USE_REHEAT:
        print(f"  P_reheat          = {det['P_rh']/1e5:.1f}  bar")
        print(f"  T_sat_RH          = {det['T_sat_rh'] - 273.15:.1f}  °C")
        print(f"  dT_rh             = {x_opt[5]:.1f}  K")
        print(f"  T_rh_ex           = {det['T_rh_ex'] - 273.15:.1f}  °C")
    print(f"  T_turb_ex         = {det['T_exp_ex'] - 273.15:.1f}  °C")
    print(f"  x_turb_ex         = {det['x_exp_ex']:.4f}  [-]  (min={X_MIN})")
    print(f"  T_air_ex_HEX1     = {det['T_air_ex'] - 273.15:.1f}  °C")
    print(f"  Pinch HEX1        = {det['pinch_hex1']:.1f}  K")
    print(f"  Pinch SH          = {det['T_oil_sg']-273.15 - (det['T_sh_ex']-273.15):.1f}  K")
    print(f"  Pinch Evap        = {det['T_oil_sg']-273.15 - T_sat_opt:.1f}  K")
    if USE_REHEAT:
        print(f"  Pinch RH          = {det['T_oil_sg']-273.15 - (det['T_rh_ex']-273.15):.1f}  K")
    print(f"{'─'*54}")
    print(f"  W_net_Brayton  = {brayton['W_net']/1e3:7.2f} kW/kg/s  →  {brayton['W_net']/1e6*m_dot_air:6.2f} MW")
    if USE_REHEAT:
        r_st = det['r_st']
        print(f"  W_ST_HP        = {r_st['W_hp']/1e3:7.2f} kW/kg/s_wf")
        print(f"  W_ST_LP        = {r_st['W_lp']/1e3:7.2f} kW/kg/s_wf")
    print(f"  W_net_Rankine  = {det['W_net_st']/1e3:7.2f} kW/kg/s  →  {det['W_net_st']/1e6*m_dot_air:6.2f} MW")
    print(f"  W_net_total    = {W_net_sp:7.2f} kW/kg/s  →  {W_net_sp/1e3*m_dot_air:6.2f} MW  ✓")
    print(f"  Q_HEX1         = {det['Q_HEX1']/1e3:7.2f} kW/kg/s  →  {det['Q_HEX1']/1e6*m_dot_air:6.2f} MW")
    print(f"  Q_in_Brayton   = {brayton['Q_in']/1e3:7.2f} kW/kg/s  →  {brayton['Q_in']/1e6*m_dot_air:6.2f} MW")
    if PIPE_LOSSES:
        print(f"  Q_loss_pipe    = {det['Q_loss_sp']/1e3:7.4f} kW/kg/s  →  {det['Q_loss_sp']/1e6*m_dot_air*1e3:6.2f} kW  (total)")
    print(f"  m_dot_air      = {m_dot_air:.2f}  kg/s")
    print(f"  m_dot_oil      = {det['m_dot_oil']:.3f} kg/s per kg/s_air  →  {det['m_dot_oil']*m_dot_air:.2f}  kg/s")
    print(f"  m_dot_wf       = {det['m_dot_wf']:.3f} kg/s per kg/s_air  →  {det['m_dot_wf']*m_dot_air:.2f}  kg/s")
    print(f"{'─'*54}")
    print(f"  eta_Rankine    = {det['eta_st']*100:.2f}  %")
    print(f"  eta_Brayton    = {brayton['eta']*100:.2f}  %")
    print(f"  eta_combined   = {eta_opt*100:.2f}  %")
    print(f"  Gain           = +{(eta_opt - brayton['eta'])*100:.2f}  pp")
    print(f"{'═'*54}")


#%%
if __name__ == "__main__":

    print("=" * 50)
    print("  Step 1 — Brayton cycle")
    print("=" * 50)
    brayton = run_brayton(print_results=True)

    rh_tag = " + Reheat" if USE_REHEAT else ""
    print("\n" + "=" * 50)
    print(f"  Step 2 — Combined cycle optimizer{rh_tag}")
    print("=" * 50)

    bounds = make_bounds(brayton['T_ex_turb'])
    print(f"  T_ex_turb    = {brayton['T_ex_turb'] - 273.15:.1f} °C")
    print(f"  Bounds : T_oil_su  ∈ [{bounds[0][0]:.0f}, {bounds[0][1]:.1f}] °C")
    print(f"           T_oil_ret ∈ [{bounds[1][0]:.0f}, {bounds[1][1]:.0f}] °C")
    print(f"           P_high    ∈ [{bounds[2][0]:.0f}, {bounds[2][1]:.0f}] bar")
    print(f"           dT_sh     ∈ [{bounds[3][0]:.0f}, {bounds[3][1]:.0f}] K")
    if USE_REHEAT:
        print(f"           P_rh      ∈ [{bounds[4][0]:.0f}, {bounds[4][1]:.0f}] bar")
        print(f"           dT_rh     ∈ [{bounds[5][0]:.0f}, {bounds[5][1]:.0f}] K")
    print()

    diagnose_feasibility(brayton, bounds)

    x0, eta_x0 = find_x0(brayton, bounds)
    if x0 is None:
        print("ERROR: no feasible starting point found.")
        raise SystemExit
    print(f"  Best x0 : {x0}  → eta={eta_x0*100:.2f}%")
    print()

    result = minimize(
        objective, x0,
        args=(brayton, bounds),
        method='Nelder-Mead',
        options={'xatol': 0.1, 'fatol': 1e-5, 'maxiter': 3000, 'disp': True},
    )

    x_opt = result.x
    eta_opt, det = evaluate(x_opt, brayton, bounds)
    if det is None:
        print("Warning: optimizer returned infeasible point — using best x0.")
        eta_opt, det = evaluate(x0, brayton, bounds)
        x_opt = x0

    W_net_sp  = brayton['W_net']/1e3 + det['W_net_st']/1e3
    m_dot_air = (W_net_target / 1e3) / W_net_sp

    print_results(x_opt, eta_opt, det, brayton, m_dot_air)