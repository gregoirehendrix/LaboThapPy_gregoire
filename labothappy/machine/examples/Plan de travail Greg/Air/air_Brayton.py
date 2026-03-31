# ============================================================
#   OPEN AIR BRAYTON — sCO2 style
#   Simple | Recuperated | Recuperated + Reheat + Intercooling
# ============================================================

#%%
import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_rec  import RecursiveCircuit
from labothappy.machine.circuit_it   import IterativeCircuit
from labothappy.connector.mass_connector        import MassConnector
from labothappy.connector.solar_salt_connector  import SolarSaltConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff     import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff    import HexCstEff


# ─────────────────────────────────────────────
#   CYCLE BUILDERS
# ─────────────────────────────────────────────

#%%
def brayton_simple(eta_cp, eta_tb, eta_heater,
                   HSource, AirInlet, P_low, P_high):

    cycle      = IterativeCircuit(fluid='Air')
    compressor = CompressorCstEff()
    heater     = HexCstEff()
    turbine    = ExpanderCstEff()

    compressor.set_parameters(eta_is=eta_cp)
    turbine.set_parameters(eta_is=eta_tb)
    heater.set_parameters(eta=eta_heater)

    cycle.add_component(compressor, "Compressor")
    cycle.add_component(heater,     "Heater")
    cycle.add_component(turbine,    "Turbine")

    cycle.link_components("Compressor", "m-ex",   "Heater",  "m-su_C")
    cycle.link_components("Heater",     "m-ex_C", "Turbine", "m-su")

    cycle.add_source("AirInlet",  AirInlet, cycle.components["Compressor"], "m-su")
    cycle.add_source("HotSource", HSource,  cycle.components["Heater"],     "m-su_H")

    cycle.set_cycle_input(target="Compressor:ex", p=P_high)
    cycle.set_cycle_input(target="Turbine:ex",    p=P_low)

    cycle._build_solve_order()
    for name in cycle.solve_start_components:
        cycle.components[name].solve()

    return cycle, compressor, turbine, heater


#%%
def brayton_recuperated(eta_cp, eta_tb, eta_heater, eta_recup,
                        HSource, AirInlet, P_low, P_high, T_tb_su_guess):

    cycle       = RecursiveCircuit('Air')
    compressor  = CompressorCstEff()
    recuperator = HexCstEff()
    heater      = HexCstEff()
    turbine     = ExpanderCstEff()

    compressor.set_parameters(eta_is=eta_cp)
    turbine.set_parameters(eta_is=eta_tb)
    heater.set_parameters(eta=eta_heater)
    recuperator.set_parameters(eta=eta_recup)

    cycle.add_component(compressor,  "Compressor")
    cycle.add_component(recuperator, "Recuperator")
    cycle.add_component(heater,      "Heater")
    cycle.add_component(turbine,     "Turbine")

    cycle.link_components("Compressor",  "m-ex",   "Recuperator", "m-su_C")
    cycle.link_components("Recuperator", "m-ex_C", "Heater",      "m-su_C")
    cycle.link_components("Heater",      "m-ex_C", "Turbine",     "m-su")
    cycle.link_components("Turbine",     "m-ex",   "Recuperator", "m-su_H")

    cycle.add_source("AirInlet",  AirInlet, cycle.components["Compressor"], "m-su")
    cycle.add_source("HotSource", HSource,  cycle.components["Heater"],     "m-su_H")

    cycle.set_fixed_properties(target="Compressor:ex", p=P_high)
    cycle.set_fixed_properties(target="Turbine:ex",    p=P_low)

    T_amb = AirInlet.T
    m_dot = AirInlet.m_dot

    h_cs    = PropsSI('H', 'T', T_amb,        'P', P_low,  'Air')
    s_cs    = PropsSI('S', 'T', T_amb,        'P', P_low,  'Air')
    h_cx_is = PropsSI('H', 'P', P_high,       'S', s_cs,   'Air')
    h_cx    = h_cs + (h_cx_is - h_cs) / eta_cp
    T_cx    = PropsSI('T', 'H', h_cx,         'P', P_high, 'Air')

    h_ts    = PropsSI('H', 'T', T_tb_su_guess,'P', P_high, 'Air')
    s_ts    = PropsSI('S', 'T', T_tb_su_guess,'P', P_high, 'Air')
    h_tx_is = PropsSI('H', 'P', P_low,        'S', s_ts,   'Air')
    h_tx    = h_ts - eta_tb * (h_ts - h_tx_is)
    T_tx    = PropsSI('T', 'H', h_tx,         'P', P_low,  'Air')

    T_rec_xC = T_cx + eta_recup * (T_tx - T_cx)
    T_rec_xH = T_tx - eta_recup * (T_tx - T_cx)

    cycle.set_cycle_guess(target='Recuperator:su_C', T=T_cx,     p=P_high, m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:su_H', T=T_tx,     p=P_low,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:ex_C', T=T_rec_xC, p=P_high)
    cycle.set_cycle_guess(target='Recuperator:ex_H', T=T_rec_xH, p=P_low)

    cycle.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance=1e-3)

    cycle.mute_print()
    cycle.solve(max_iter=50)

    return cycle, compressor, turbine, heater, recuperator


#%%
def brayton_recuperated_rh_ic(eta_cp, eta_tb, eta_heater, eta_recup, eta_cooler,
                               HSource_heater, HSource_reheater, AirInlet, CSource,
                               P_low, P_high, P_mid, T_tb_su_guess):
    """
    Recuperated Brayton with 1 reheat + 1 intercooling.

    Layout:
      AirInlet → Comp1 → Intercooler → Comp2 → Recuperator(C) → Heater
               → Turb1 → Reheater   → Turb2 → Recuperator(H) → exhaust

    Pressures:
      P_low  : Comp1 inlet  / Turb2 outlet
      P_mid  : Comp1 outlet / Turb1 outlet  (= sqrt(P_low × P_high))
      P_high : Comp2 outlet / Turb1 inlet
    """
    cycle       = RecursiveCircuit('Air')
    comp1       = CompressorCstEff()
    comp2       = CompressorCstEff()
    intercooler = HexCstEff()
    recuperator = HexCstEff()
    heater      = HexCstEff()
    turb1       = ExpanderCstEff()
    reheater    = HexCstEff()
    turb2       = ExpanderCstEff()

    comp1.set_parameters(eta_is=eta_cp)
    comp2.set_parameters(eta_is=eta_cp)
    turb1.set_parameters(eta_is=eta_tb)
    turb2.set_parameters(eta_is=eta_tb)
    heater.set_parameters(eta=eta_heater)
    reheater.set_parameters(eta=eta_heater)
    recuperator.set_parameters(eta=eta_recup)
    intercooler.set_parameters(eta=eta_cooler)

    cycle.add_component(comp1,       "Comp1")
    cycle.add_component(intercooler, "Intercooler")
    cycle.add_component(comp2,       "Comp2")
    cycle.add_component(recuperator, "Recuperator")
    cycle.add_component(heater,      "Heater")
    cycle.add_component(turb1,       "Turb1")
    cycle.add_component(reheater,    "Reheater")
    cycle.add_component(turb2,       "Turb2")

    cycle.link_components("Comp1",       "m-ex",   "Intercooler", "m-su_H")
    cycle.link_components("Intercooler", "m-ex_H", "Comp2",       "m-su")
    cycle.link_components("Comp2",       "m-ex",   "Recuperator", "m-su_C")
    cycle.link_components("Recuperator", "m-ex_C", "Heater",      "m-su_C")
    cycle.link_components("Heater",      "m-ex_C", "Turb1",       "m-su")
    cycle.link_components("Turb1",       "m-ex",   "Reheater",    "m-su_C")
    cycle.link_components("Reheater",    "m-ex_C", "Turb2",       "m-su")
    cycle.link_components("Turb2",       "m-ex",   "Recuperator", "m-su_H")

    cycle.add_source("AirInlet",   AirInlet,         cycle.components["Comp1"],       "m-su")
    cycle.add_source("HotSource",  HSource_heater,   cycle.components["Heater"],      "m-su_H")
    cycle.add_source("ReheatSrc",  HSource_reheater, cycle.components["Reheater"],    "m-su_H")
    cycle.add_source("ColdSource", CSource,           cycle.components["Intercooler"], "m-su_C")

    cycle.set_fixed_properties(target="Comp1:ex", p=P_mid)
    cycle.set_fixed_properties(target="Comp2:ex", p=P_high)
    cycle.set_fixed_properties(target="Turb1:ex", p=P_mid)
    cycle.set_fixed_properties(target="Turb2:ex", p=P_low)

    T_amb = AirInlet.T
    m_dot = AirInlet.m_dot

    h_c1s    = PropsSI('H', 'T', T_amb,        'P', P_low,  'Air')
    s_c1s    = PropsSI('S', 'T', T_amb,        'P', P_low,  'Air')
    h_c1x_is = PropsSI('H', 'P', P_mid,        'S', s_c1s,  'Air')
    h_c1x    = h_c1s + (h_c1x_is - h_c1s) / eta_cp
    T_c1x    = PropsSI('T', 'H', h_c1x,        'P', P_mid,  'Air')

    T_ic_ex  = T_amb + (1 - eta_cooler) * (T_c1x - T_amb)

    h_c2s    = PropsSI('H', 'T', T_ic_ex,      'P', P_mid,  'Air')
    s_c2s    = PropsSI('S', 'T', T_ic_ex,      'P', P_mid,  'Air')
    h_c2x_is = PropsSI('H', 'P', P_high,       'S', s_c2s,  'Air')
    h_c2x    = h_c2s + (h_c2x_is - h_c2s) / eta_cp
    T_c2x    = PropsSI('T', 'H', h_c2x,        'P', P_high, 'Air')

    h_t1s    = PropsSI('H', 'T', T_tb_su_guess,'P', P_high, 'Air')
    s_t1s    = PropsSI('S', 'T', T_tb_su_guess,'P', P_high, 'Air')
    h_t1x_is = PropsSI('H', 'P', P_mid,        'S', s_t1s,  'Air')
    h_t1x    = h_t1s - eta_tb * (h_t1s - h_t1x_is)
    T_t1x    = PropsSI('T', 'H', h_t1x,        'P', P_mid,  'Air')

    h_t2s    = PropsSI('H', 'T', T_tb_su_guess,'P', P_mid,  'Air')
    s_t2s    = PropsSI('S', 'T', T_tb_su_guess,'P', P_mid,  'Air')
    h_t2x_is = PropsSI('H', 'P', P_low,        'S', s_t2s,  'Air')
    h_t2x    = h_t2s - eta_tb * (h_t2s - h_t2x_is)
    T_t2x    = PropsSI('T', 'H', h_t2x,        'P', P_low,  'Air')

    T_rec_xC = T_c2x + eta_recup * (T_t2x - T_c2x)
    T_rec_xH = T_t2x - eta_recup * (T_t2x - T_c2x)

    cycle.set_cycle_guess(target='Comp1:su',         T=T_amb,        p=P_low,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Comp1:ex',         T=T_c1x,        p=P_mid)
    cycle.set_cycle_guess(target='Comp2:su',         T=T_ic_ex,      p=P_mid,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Comp2:ex',         T=T_c2x,        p=P_high)
    cycle.set_cycle_guess(target='Recuperator:su_C', T=T_c2x,        p=P_high, m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:ex_C', T=T_rec_xC,     p=P_high)
    cycle.set_cycle_guess(target='Recuperator:su_H', T=T_t2x,        p=P_low,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Recuperator:ex_H', T=T_rec_xH,     p=P_low)
    cycle.set_cycle_guess(target='Turb1:su',         T=T_tb_su_guess,p=P_high, m_dot=m_dot)
    cycle.set_cycle_guess(target='Turb1:ex',         T=T_t1x,        p=P_mid)
    cycle.set_cycle_guess(target='Turb2:su',         T=T_tb_su_guess,p=P_mid,  m_dot=m_dot)
    cycle.set_cycle_guess(target='Turb2:ex',         T=T_t2x,        p=P_low)

    cycle.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Turb1:ex',         variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Turb1:ex',         variable='p', tolerance=1e-3)
    cycle.set_residual_variable(target='Comp2:su',         variable='h', tolerance=1e-3)
    cycle.set_residual_variable(target='Comp2:su',         variable='p', tolerance=1e-3)

    cycle.mute_print()
    cycle.solve(max_iter=50)

    return cycle, comp1, comp2, turb1, turb2, heater, reheater, recuperator, intercooler


# ─────────────────────────────────────────────
#   PERFORMANCE
# ─────────────────────────────────────────────

#%%
def compute_cycle_performance(compressor, turbine, heater,
                               hot_fluid, T_hot_su, P_hot,
                               recuperator=None,
                               comp2=None, turb2=None,
                               reheater=None, intercooler=None):
    m_dot = compressor.su.m_dot

    W_comp  = m_dot * (compressor.ex.h - compressor.su.h) / 1000
    W_turb  = m_dot * (turbine.su.h    - turbine.ex.h)    / 1000
    W_comp2 = m_dot * (comp2.ex.h  - comp2.su.h)  / 1000 if comp2  is not None else 0.0
    W_turb2 = m_dot * (turb2.su.h  - turb2.ex.h)  / 1000 if turb2  is not None else 0.0

    W_net   = (W_turb + W_turb2) - (W_comp + W_comp2)
    Q_in    =  m_dot * (heater.ex_C.h    - heater.su_C.h)    / 1000
    Q_rh    = (m_dot * (reheater.ex_C.h  - reheater.su_C.h)  / 1000
               if reheater    is not None else 0.0)
    Q_ic    = (m_dot * (intercooler.su_H.h - intercooler.ex_H.h) / 1000
               if intercooler is not None else 0.0)
    Q_total = Q_in + Q_rh
    Q_cooler = Q_total - W_net   # first-law derived

    # hot-side mass flow from energy balance (HSource m_dot is a dummy)
    if hot_fluid == 'SolarSalt':
        Cp_avg    = SolarSaltConnector._Cp((T_hot_su + heater.ex_H.T) / 2)
        m_dot_hot = (Q_in * 1000) / (Cp_avg * (T_hot_su - heater.ex_H.T))
        if reheater is not None:
            Cp_rh      = SolarSaltConnector._Cp((T_hot_su + reheater.ex_H.T) / 2)
            m_dot_hot += (Q_rh * 1000) / (Cp_rh * (T_hot_su - reheater.ex_H.T))
        hot_label = 'm_dot_salt'
    else:
        h_hin     = PropsSI('H', 'T', T_hot_su,       'P', P_hot, 'Air')
        h_hex     = PropsSI('H', 'T', heater.ex_H.T,  'P', P_hot, 'Air')
        m_dot_hot = (Q_in * 1000) / (h_hin - h_hex)
        if reheater is not None:
            h_hex_rh   = PropsSI('H', 'T', reheater.ex_H.T, 'P', P_hot, 'Air')
            m_dot_hot += (Q_rh * 1000) / (h_hin - h_hex_rh)
        hot_label = 'm_dot_air_hot'

    # exhaust temperature (after recuperator if present)
    if recuperator is not None:
        T_exhaust = recuperator.ex_H.T
    elif turb2 is not None:
        T_exhaust = turb2.ex.T
    else:
        T_exhaust = turbine.ex.T

    result = {
        'W_comp':    W_comp,   'W_comp2':   W_comp2,
        'W_turb':    W_turb,   'W_turb2':   W_turb2,
        'W_net':     W_net,
        'Q_in':      Q_in,     'Q_rh':      Q_rh,     'Q_ic': Q_ic,
        'Q_total':   Q_total,  'Q_cooler':  Q_cooler,
        'eta':       W_net / Q_total,
        'T_exhaust': T_exhaust,
        'm_dot_hot': m_dot_hot, 'hot_label': hot_label,
    }
    if recuperator is not None:
        result['Q_recup'] = m_dot * (recuperator.ex_C.h - recuperator.su_C.h) / 1000

    return result


# ─────────────────────────────────────────────
#   PRINT HELPERS
# ─────────────────────────────────────────────

#%%
def print_states(compressor, turbine,
                 recuperator=None,
                 comp2=None, turb2=None,
                 intercooler=None, reheater=None):
    print("=== Comp1:su ===")
    print(f"  T     = {compressor.su.T - 273.15:.2f} °C")
    print(f"  p     = {compressor.su.p / 1e5:.3f} bar")
    print(f"  m_dot = {compressor.su.m_dot:.4f} kg/s")
    print("=== Comp1:ex ===")
    print(f"  T     = {compressor.ex.T - 273.15:.2f} °C")
    print(f"  p     = {compressor.ex.p / 1e5:.3f} bar")
    if intercooler is not None:
        print("=== Intercooler ===")
        print(f"  su_H T = {intercooler.su_H.T - 273.15:.2f} °C")
        print(f"  ex_H T = {intercooler.ex_H.T - 273.15:.2f} °C")
    if comp2 is not None:
        print("=== Comp2:su ===")
        print(f"  T     = {comp2.su.T - 273.15:.2f} °C")
        print(f"  p     = {comp2.su.p / 1e5:.3f} bar")
        print("=== Comp2:ex ===")
        print(f"  T     = {comp2.ex.T - 273.15:.2f} °C")
        print(f"  p     = {comp2.ex.p / 1e5:.3f} bar")
    if recuperator is not None:
        print("=== Recuperator ===")
        print(f"  su_C T = {recuperator.su_C.T - 273.15:.2f} °C")
        print(f"  ex_C T = {recuperator.ex_C.T - 273.15:.2f} °C")
        print(f"  su_H T = {recuperator.su_H.T - 273.15:.2f} °C")
        print(f"  ex_H T = {recuperator.ex_H.T - 273.15:.2f} °C")
    print("=== Turb1:su ===")
    print(f"  T     = {turbine.su.T - 273.15:.2f} °C")
    print(f"  p     = {turbine.su.p / 1e5:.3f} bar")
    print("=== Turb1:ex ===")
    print(f"  T     = {turbine.ex.T - 273.15:.2f} °C")
    print(f"  p     = {turbine.ex.p / 1e5:.3f} bar")
    if reheater is not None:
        print("=== Reheater ===")
        print(f"  su_C T = {reheater.su_C.T - 273.15:.2f} °C")
        print(f"  ex_C T = {reheater.ex_C.T - 273.15:.2f} °C")
        print(f"  su_H T = {reheater.su_H.T - 273.15:.2f} °C")
        print(f"  ex_H T = {reheater.ex_H.T - 273.15:.2f} °C")
    if turb2 is not None:
        print("=== Turb2:su ===")
        print(f"  T     = {turb2.su.T - 273.15:.2f} °C")
        print(f"  p     = {turb2.su.p / 1e5:.3f} bar")
        print("=== Turb2:ex ===")
        print(f"  T     = {turb2.ex.T - 273.15:.2f} °C")
        print(f"  p     = {turb2.ex.p / 1e5:.3f} bar")


#%%
def print_efficiency(perf, First_Law_Check=False):
    print("=== Cycle Efficiency ===")
    print(f"  W_comp   : {perf['W_comp']:.2f} kW")
    if perf['W_comp2'] > 0:
        print(f"  W_comp2  : {perf['W_comp2']:.2f} kW")
    print(f"  W_turb   : {perf['W_turb']:.2f} kW")
    if perf['W_turb2'] > 0:
        print(f"  W_turb2  : {perf['W_turb2']:.2f} kW")
    print(f"  W_net    : {perf['W_net']:.2f} kW")
    print(f"  Q_in     : {perf['Q_in']:.2f} kW")
    if perf['Q_rh'] > 0:
        print(f"  Q_reheat : {perf['Q_rh']:.2f} kW")
    if perf['Q_ic'] > 0:
        print(f"  Q_ic     : {perf['Q_ic']:.2f} kW")
    if 'Q_recup' in perf:
        print(f"  Q_recup  : {perf['Q_recup']:.2f} kW")
    print(f"  Q_total  : {perf['Q_total']:.2f} kW")
    print(f"  Q_cooler : {perf['Q_cooler']:.2f} kW")
    print(f"  T_exhaust: {perf['T_exhaust'] - 273.15:.2f} °C")
    print(f"  eta      : {perf['eta'] * 100:.2f} %")
    print(f"  {perf['hot_label']} : {perf['m_dot_hot']:.4f} kg/s / kg_air/s")
    if First_Law_Check:
        balance = perf['Q_total'] - perf['Q_cooler']
        print("  --- First law check ---")
        print(f"  Q_total - Q_cooler  : {balance:.4f} kW")
        print(f"  W_net               : {perf['W_net']:.4f} kW")
        print(f"  Discrepancy         : {abs(perf['W_net'] - balance):.6f} kW")
    print("========================")


#%%
def scale_to_power(W_net_target_MW, perf, PRINT_SCALE=True, First_Law_Check=False):
    """Scale all quantities to target net power. Reference run at m_dot_air = 1 kg/s."""
    scale = (W_net_target_MW * 1000) / perf['W_net']
    if PRINT_SCALE:
        print("=== Scaling to target net power ===")
        print(f"  W_net_target       : {W_net_target_MW:.2f} MW")
        print(f"  W_net ref (1 kg/s) : {perf['W_net']:.2f} kW")
        print(f"  Scale factor       : {scale:.4f}")
        print(f"  m_dot_air required : {scale:.4f} kg/s")
        print(f"  W_comp             : {perf['W_comp']   * scale / 1000:.3f} MW")
        if perf['W_comp2'] > 0:
            print(f"  W_comp2            : {perf['W_comp2'] * scale / 1000:.3f} MW")
        print(f"  W_turb             : {perf['W_turb']   * scale / 1000:.3f} MW")
        if perf['W_turb2'] > 0:
            print(f"  W_turb2            : {perf['W_turb2'] * scale / 1000:.3f} MW")
        print(f"  W_net              : {perf['W_net']    * scale / 1000:.3f} MW")
        print(f"  Q_total            : {perf['Q_total']  * scale / 1000:.3f} MW")
        print(f"  Q_cooler           : {perf['Q_cooler'] * scale / 1000:.3f} MW")
        if 'Q_recup' in perf:
            print(f"  Q_recup            : {perf['Q_recup']  * scale / 1000:.3f} MW")
        print(f"  {perf['hot_label']:<18} : {perf['m_dot_hot'] * scale:.3f} kg/s")
        print(f"  eta (thermo)       : {perf['eta'] * 100:.2f} %")
        if First_Law_Check:
            Q_h = perf['Q_total']  * scale / 1000
            Q_c = perf['Q_cooler'] * scale / 1000
            W   = perf['W_net']    * scale / 1000
            print("  --- First law check ---")
            print(f"  Q_total - Q_cooler  : {Q_h - Q_c:.6f} MW")
            print(f"  W_net               : {W:.6f} MW")
            print(f"  Discrepancy         : {abs(W - (Q_h - Q_c)):.8f} MW")
        print("====================================")
    return scale


# ─────────────────────────────────────────────
#   MAIN
# ─────────────────────────────────────────────

#%%
if __name__ == "__main__":

    study_case      = "Recuperated_RH_IC"  # "Simple" | "Recuperated" | "Recuperated_RH_IC"
    PRINT           = True
    DETAIL          = True
    PRINT_SCALE     = True
    First_Law_Check = False
    W_net_target    = 5                    # [MW]

    # ── Operating conditions ─────────────────
    T_amb     = 13.2 + 273.15   # [K]
    P_low     = 1.01325e5       # [Pa]
    PR        = 4
    P_high    = P_low * PR
    P_mid     = np.sqrt(P_low * P_high)   # optimal intermediate pressure [Pa]
    m_dot_air = 1.0             # [kg/s] reference

    # ── Efficiencies ─────────────────────────
    eta_cp     = 0.80
    eta_tb     = 0.90
    eta_heater = 1.00
    eta_recup  = 0.85
    eta_cooler = 0.99

    # ── Hot source ───────────────────────────
    T_hot_su     = 850 + 273.15   # [K]
    T_salt_limit = 565 + 273.15   # [K]
    P_hot        = 1e5            # [Pa]

    if T_hot_su <= T_salt_limit:
        HSource          = SolarSaltConnector()
        HSource.set_properties(T=T_hot_su, p=P_hot, m_dot=1.0)
        HSource_heater   = SolarSaltConnector()
        HSource_heater.set_properties(T=T_hot_su, p=P_hot, m_dot=1.0)
        HSource_reheater = SolarSaltConnector()
        HSource_reheater.set_properties(T=T_hot_su, p=P_hot, m_dot=1.0)
        hot_fluid = 'SolarSalt'
        print(f"Hot source : Solar Salt at {T_hot_su - 273.15:.1f} °C")
    else:
        HSource          = MassConnector()
        HSource.set_properties(fluid='Air', T=T_hot_su, p=P_hot, m_dot=1.0)
        HSource_heater   = MassConnector()
        HSource_heater.set_properties(fluid='Air', T=T_hot_su, p=P_hot, m_dot=1.0)
        HSource_reheater = MassConnector()
        HSource_reheater.set_properties(fluid='Air', T=T_hot_su, p=P_hot, m_dot=1.0)
        hot_fluid = 'Air'
        print(f"Hot source : Air at {T_hot_su - 273.15:.1f} °C")

    # ── Air inlet & cold source ───────────────
    AirInlet = MassConnector()
    AirInlet.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=m_dot_air)

    CSource = MassConnector()
    CSource.set_properties(fluid='Air', T=T_amb, P=P_low, m_dot=10.0)

    # ── Run ──────────────────────────────────
    if study_case == "Simple":
        cycle, comp, turb, heater = brayton_simple(
            eta_cp, eta_tb, eta_heater,
            HSource, AirInlet, P_low, P_high
        )
        perf = compute_cycle_performance(
            comp, turb, heater, hot_fluid, T_hot_su, P_hot
        )
        if PRINT:
            if DETAIL: print_states(comp, turb)
            print_efficiency(perf, First_Law_Check)
        scale_to_power(W_net_target, perf, PRINT_SCALE, First_Law_Check)

    elif study_case == "Recuperated":
        cycle, comp, turb, heater, recup = brayton_recuperated(
            eta_cp, eta_tb, eta_heater, eta_recup,
            HSource, AirInlet, P_low, P_high,
            T_tb_su_guess=T_hot_su
        )
        perf = compute_cycle_performance(
            comp, turb, heater, hot_fluid, T_hot_su, P_hot,
            recuperator=recup
        )
        if PRINT:
            if DETAIL: print_states(comp, turb, recuperator=recup)
            print_efficiency(perf, First_Law_Check)
        scale_to_power(W_net_target, perf, PRINT_SCALE, First_Law_Check)

    elif study_case == "Recuperated_RH_IC":
        cycle, comp1, comp2, turb1, turb2, heater, reheater, recup, intercooler = \
            brayton_recuperated_rh_ic(
                eta_cp, eta_tb, eta_heater, eta_recup, eta_cooler,
                HSource_heater, HSource_reheater, AirInlet, CSource,
                P_low, P_high, P_mid,
                T_tb_su_guess=T_hot_su
            )
        perf = compute_cycle_performance(
            comp1, turb1, heater, hot_fluid, T_hot_su, P_hot,
            recuperator=recup,
            comp2=comp2, turb2=turb2,
            reheater=reheater, intercooler=intercooler
        )
        if PRINT:
            if DETAIL:
                print_states(comp1, turb1, recuperator=recup,
                             comp2=comp2, turb2=turb2,
                             intercooler=intercooler, reheater=reheater)
            print_efficiency(perf, First_Law_Check)
        scale_to_power(W_net_target, perf, PRINT_SCALE, First_Law_Check)