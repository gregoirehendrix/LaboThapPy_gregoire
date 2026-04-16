# -*- coding: utf-8 -*-
"""
Created on Apr 2026
@author: gregoire.hendrix

Steam Rankine — Reheat + 2 cFWH LP (drain cascaded back) + DEA open FWH

Architecture :
    Passe 1 — RecursiveCircuit :
        Pump → Eco → Eva → SH → Turb_HP → RH → Turb_LP → Condenser
        f_rh=0.10, P_rh=25 bar (optimum)

    Passe 2 — Analytique :
        Condenser → Pump_LP(1kg/s) → cFWH_LP1 → cFWH_LP2 → DEA → Pump_HP → Eco

        Turb_LP : P_rh → P_fw → P_s2 → P_s1 → P_low
            LP1 : 1 kg/s          soutirage DEA   @ P_fw
            LP2 : (1-y_dea)       soutirage LP2   @ P_s2
            LP3 : (1-y_dea-y_lp2) soutirage LP1   @ P_s1
            LP4 : m_lp4           → condenseur

        Drains cascaded back vers condenseur.

    Pressions :
        P_fw  = 1.00 bar  (DEA)
        P_s2  = sqrt(P_low * P_fw)  ≈ 0.067 bar
        P_s1  = P_low^0.75 * P_fw^0.25 ≈ 0.030 bar
"""

import warnings
warnings.filterwarnings('ignore')

from CoolProp.CoolProp import PropsSI
from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.expander.expander_csteff       import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc
from labothappy.component.heat_exchanger.hex_cstpinch    import HexCstPinch
from labothappy.component.pump.pump_csteff               import PumpCstEff


# =============================================================================
# HELPERS
# =============================================================================

def _salt_connector(T_K, p_Pa, m_dot):
    c = SolarSaltConnector()
    c.set_properties(T=T_K, p=p_Pa, m_dot=m_dot)
    return c

def _make_disc_hx(eta_max=0.99, pinch_min=5.0, n_disc=10):
    hx = HexCstEffDisc()
    hx.set_parameters(eta_max=eta_max, Pinch_min=pinch_min, n_disc=n_disc)
    return hx

def _make_pinch_cond(pinch=5.0, dT_sc=2.0):
    hx = HexCstPinch()
    hx.set_parameters(Pinch=pinch, Delta_T_sh_sc=dT_sc, HX_type='condenser')
    return hx

def _expand(h_su, p_su, p_ex, eta):
    s_su    = PropsSI('S', 'H', h_su, 'P', p_su, 'Water')
    h_ex_is = PropsSI('H', 'P', p_ex, 'S', s_su, 'Water')
    return h_su - eta * (h_su - h_ex_is)

def _pump(h_su, p_su, p_ex, eta):
    s_su    = PropsSI('S', 'H', h_su, 'P', p_su, 'Water')
    h_ex_is = PropsSI('H', 'P', p_ex, 'S', s_su, 'Water')
    return h_su + (h_ex_is - h_su) / eta


# =============================================================================
# PASSE 1 — RecursiveCircuit reheat
# =============================================================================

def _build_reheat_circuit(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                           HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh,
                           CSource, P_low, P_high, P_reheat, m_dot_st, n_disc=10):

    T_sat_hi = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_salt   = HSource_sh.T

    cycle     = RecursiveCircuit('Water')
    pump      = PumpCstEff();      pump.set_parameters(eta_is=eta_pump)
    eco       = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    eva       = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    sh        = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    rh        = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    turb_hp   = ExpanderCstEff(); turb_hp.set_parameters(eta_is=eta_turb)
    turb_lp   = ExpanderCstEff(); turb_lp.set_parameters(eta_is=eta_turb)
    condenser = _make_pinch_cond(pinch=pinch_cond)

    for name, comp in [("Pump", pump), ("Economiser", eco), ("Evaporator", eva),
                        ("Superheater", sh), ("Turb_HP", turb_hp), ("Reheater", rh),
                        ("Turb_LP", turb_lp), ("Condenser", condenser)]:
        cycle.add_component(comp, name)

    cycle.link_components("Pump",        "m-ex",   "Economiser",  "m-su_C")
    cycle.link_components("Economiser",  "m-ex_C", "Evaporator",  "m-su_C")
    cycle.link_components("Evaporator",  "m-ex_C", "Superheater", "m-su_C")
    cycle.link_components("Superheater", "m-ex_C", "Turb_HP",     "m-su")
    cycle.link_components("Turb_HP",     "m-ex",   "Reheater",    "m-su_C")
    cycle.link_components("Reheater",    "m-ex_C", "Turb_LP",     "m-su")
    cycle.link_components("Turb_LP",     "m-ex",   "Condenser",   "m-su_H")
    cycle.link_components("Condenser",   "m-ex_H", "Pump",        "m-su")

    cycle.add_source("SaltSource_SH",  HSource_sh,  cycle.components["Superheater"], "m-su_H")
    cycle.add_source("SaltSource_EVA", Salt_eva_su, cycle.components["Evaporator"],  "m-su_H")
    cycle.add_source("SaltSource_ECO", Salt_eco_su, cycle.components["Economiser"],  "m-su_H")
    cycle.add_source("SaltSource_RH",  HSource_rh,  cycle.components["Reheater"],    "m-su_H")
    cycle.add_source("ColdSource",     CSource,     cycle.components["Condenser"],   "m-su_C")

    cycle.set_fixed_properties(target="Pump:ex",    p=P_high)
    cycle.set_fixed_properties(target="Turb_HP:ex", p=P_reheat)
    cycle.set_fixed_properties(target="Turb_LP:ex", p=P_low)

    h_sat_lo  = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    h_pump_ex = _pump(h_sat_lo, P_low, P_high, eta_pump)
    T_pump_ex = PropsSI('T', 'H', h_pump_ex, 'P', P_high, 'Water')
    T_HP_su_g = T_salt - pinch_hx
    h_HP_su   = PropsSI('H', 'T', T_HP_su_g, 'P', P_high,   'Water')
    h_HP_ex   = _expand(h_HP_su, P_high, P_reheat, eta_turb)
    T_HP_ex   = PropsSI('T', 'H', h_HP_ex, 'P', P_reheat, 'Water')
    h_LP_su   = PropsSI('H', 'T', T_HP_su_g, 'P', P_reheat, 'Water')
    h_LP_ex   = _expand(h_LP_su, P_reheat, P_low, eta_turb)
    T_LP_ex   = PropsSI('T', 'H', h_LP_ex, 'P', P_low, 'Water')

    cycle.set_cycle_guess(target='Pump:su',          p=P_low,    SC=3,          m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',          p=P_high,   T=T_pump_ex)
    cycle.set_cycle_guess(target='Economiser:ex_C',  p=P_high,   T=T_sat_hi-5)
    cycle.set_cycle_guess(target='Evaporator:ex_C',  p=P_high,   T=T_sat_hi,   x=1)
    cycle.set_cycle_guess(target='Superheater:ex_C', p=P_high,   T=T_HP_su_g)
    cycle.set_cycle_guess(target='Turb_HP:su',       p=P_high,   T=T_HP_su_g,  m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_HP:ex',       p=P_reheat, T=T_HP_ex)
    cycle.set_cycle_guess(target='Reheater:ex_C',    p=P_reheat, T=T_HP_su_g)
    cycle.set_cycle_guess(target='Turb_LP:su',       p=P_reheat, T=T_HP_su_g,  m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_LP:ex',       p=P_low,    T=T_LP_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',   p=P_low,    T=T_sat_lo-2)

    for target, var in [
        ('Turb_HP:ex','h'), ('Turb_HP:ex','p'),
        ('Turb_LP:ex','h'), ('Turb_LP:ex','p'),
        ('Pump:ex','h'),    ('Pump:ex','p'),
        ('Superheater:ex_C','h'), ('Superheater:ex_C','p'),
        ('Reheater:ex_C','h'),    ('Reheater:ex_C','p'),
        ('Evaporator:ex_C','h'),  ('Evaporator:ex_C','p'),
        ('Economiser:ex_C','h'),  ('Economiser:ex_C','p'),
        ('Condenser:ex_H','h'),   ('Condenser:ex_H','p'),
    ]:
        cycle.set_residual_variable(target=target, variable=var, tolerance=1e-6)

    return cycle, pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser


def solve_reheat(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                 T_salt_su, P_salt, m_dot_salt_total, f_rh,
                 CSource, P_low, P_high, P_reheat, m_dot_st,
                 tol_outer=0.1, max_outer=20, n_disc=10):

    T_sat_hi     = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    m_dot_boiler = m_dot_salt_total * (1.0 - f_rh)
    m_dot_rh     = m_dot_salt_total * f_rh
    HSource_rh   = _salt_connector(T_salt_su, P_salt, m_dot_rh)

    T_eva_su = T_salt_su - 10.0
    T_eco_su = T_sat_hi  + 5.0

    for _ in range(max_outer):
        HSource_sh  = _salt_connector(T_salt_su, P_salt, m_dot_boiler)
        Salt_eva_su = _salt_connector(T_eva_su,  P_salt, m_dot_boiler)
        Salt_eco_su = _salt_connector(T_eco_su,  P_salt, m_dot_boiler)

        cycle, pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser = \
            _build_reheat_circuit(
                eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh,
                CSource, P_low, P_high, P_reheat, m_dot_st, n_disc,
            )
        cycle.solve()

        T_eva_su_new = sh.ex_H.T
        T_eco_su_new = eva.ex_H.T
        err = max(abs(T_eva_su_new - T_eva_su), abs(T_eco_su_new - T_eco_su))
        T_eva_su = T_eva_su_new
        T_eco_su = T_eco_su_new
        if err < tol_outer:
            break

    return pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser


# =============================================================================
# PASSE 2 — Analytique : 2 cFWH LP + DEA
# =============================================================================

def compute_reheat_cfwh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                         eta_pump, eta_turb, P_high, P_low, P_fw, P_s1, P_s2, eta_g,
                         pinch_cfwh=5.0):

    h_lp_su   = rh.ex_C.h
    p_lp_su   = turb_lp.su.p
    h_cond_ex = condenser.ex_H.h
    h_sh_ex   = sh.ex_C.h

    # Turb_LP : P_rh → P_fw → P_s2 → P_s1 → P_low
    h_ext_fw  = _expand(h_lp_su,  p_lp_su, P_fw,  eta_turb)
    h_ext_s2  = _expand(h_ext_fw, P_fw,    P_s2,  eta_turb)
    h_ext_s1  = _expand(h_ext_s2, P_s2,    P_s1,  eta_turb)
    h_lp4_ex  = _expand(h_ext_s1, P_s1,    P_low, eta_turb)

    T_ext_fw  = PropsSI('T', 'H', h_ext_fw, 'P', P_fw,  'Water')
    T_ext_s2  = PropsSI('T', 'H', h_ext_s2, 'P', P_s2,  'Water')
    T_ext_s1  = PropsSI('T', 'H', h_ext_s1, 'P', P_s1,  'Water')

    # Pump_LP : 1 kg/s, P_low → P_fw
    h_plp_ex   = _pump(h_cond_ex, P_low, P_fw, eta_pump)

    # Drains sat liq
    h_drain_s1 = PropsSI('H', 'P', P_s1, 'Q', 0, 'Water')
    h_drain_s2 = PropsSI('H', 'P', P_s2, 'Q', 0, 'Water')
    h_fw_out   = PropsSI('H', 'P', P_fw, 'Q', 0, 'Water')
    T_fw_out   = PropsSI('T', 'P', P_fw, 'Q', 0, 'Water')

    # Sorties eau côté tube cFWH (T_sat - pinch, à pression P_fw)
    T_fw_out_lp1 = PropsSI('T', 'P', P_s1, 'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_lp2 = PropsSI('T', 'P', P_s2, 'Q', 0, 'Water') - pinch_cfwh
    h_fw_out_lp1 = PropsSI('H', 'T', T_fw_out_lp1, 'P', P_fw, 'Water')
    h_fw_out_lp2 = PropsSI('H', 'T', T_fw_out_lp2, 'P', P_fw, 'Water')

    # cFWH_LP1 : y_lp1*(h_ext_s1 - h_drain_s1) = 1*(h_fw_out_lp1 - h_plp_ex)
    y_lp1 = (h_fw_out_lp1 - h_plp_ex)     / (h_ext_s1 - h_drain_s1)

    # cFWH_LP2 : y_lp2*(h_ext_s2 - h_drain_s2) = 1*(h_fw_out_lp2 - h_fw_out_lp1)
    y_lp2 = (h_fw_out_lp2 - h_fw_out_lp1) / (h_ext_s2 - h_drain_s2)

    # DEA : 1*h_fw_out_lp2 + y_dea*h_ext_fw = (1+y_dea)*h_fw_out  [bilan masse+énergie]
    # → y_dea*(h_ext_fw - h_fw_out) = h_fw_out - h_fw_out_lp2
    y_dea = (h_fw_out - h_fw_out_lp2) / (h_ext_fw - h_fw_out)

    if y_dea < 0:
        print(f"  WARNING : y_dea={y_dea:.4f} < 0 — revoir P_fw")

    # Pump_HP : 1 kg/s, P_fw → P_high
    h_php_ex = _pump(h_fw_out, P_fw, P_high, eta_pump)

    m_lp4 = 1.0 - y_dea - y_lp2 - y_lp1
    if m_lp4 < 0:
        print(f"  WARNING : m_lp4={m_lp4:.4f} < 0")

    # Travaux [kJ/kg_total]
    W_pump_lp = 1.0                        * (h_plp_ex        - h_cond_ex)  / 1000
    W_pump_hp = 1.0                        * (h_php_ex        - h_fw_out)   / 1000
    W_HP      = 1.0                        * (turb_hp.su.h    - turb_hp.ex.h) / 1000
    W_LP1     = 1.0                        * (h_lp_su         - h_ext_fw)   / 1000
    W_LP2     = (1 - y_dea)                * (h_ext_fw        - h_ext_s2)   / 1000
    W_LP3     = (1 - y_dea - y_lp2)        * (h_ext_s2        - h_ext_s1)   / 1000
    W_LP4     = m_lp4                      * (h_ext_s1        - h_lp4_ex)   / 1000
    W_turb    = W_HP + W_LP1 + W_LP2 + W_LP3 + W_LP4
    W_net     = W_turb - W_pump_lp - W_pump_hp

    # Chaleurs [kJ/kg_total]
    Q_boiler  = (h_sh_ex   - h_php_ex)             / 1000
    Q_rh      = (rh.ex_C.h - turb_hp.ex.h)         / 1000
    Q_boil    = Q_boiler + Q_rh
    Q_cond = Q_boil + W_pump_lp + W_pump_hp - W_turb

    eta_th      = W_net / Q_boil
    discrepancy = abs(W_net - (Q_boil - Q_cond))
    eta_el = eta_g * eta_th


    return {
        'y_dea': y_dea, 'y_lp1': y_lp1, 'y_lp2': y_lp2, 'm_lp4': m_lp4,
        'P_fw': P_fw, 'P_s1': P_s1, 'P_s2': P_s2,
        'T_fw_out': T_fw_out, 'T_ext_fw': T_ext_fw,
        'T_ext_s1': T_ext_s1, 'T_ext_s2': T_ext_s2,
        'h_plp_ex': h_plp_ex, 'h_php_ex': h_php_ex,
        'h_fw_out': h_fw_out, 'h_fw_out_lp1': h_fw_out_lp1,
        'h_fw_out_lp2': h_fw_out_lp2, 'h_lp4_ex': h_lp4_ex,
        'h_drain_s1': h_drain_s1, 'h_drain_s2': h_drain_s2,
        'W_pump_lp': W_pump_lp, 'W_pump_hp': W_pump_hp,
        'W_HP': W_HP, 'W_LP1': W_LP1, 'W_LP2': W_LP2,
        'W_LP3': W_LP3, 'W_LP4': W_LP4,
        'W_turb': W_turb, 'W_net': W_net,
        'Q_boiler': Q_boiler, 'Q_rh': Q_rh, 'Q_boil': Q_boil,
        'Q_cond': Q_cond, 'eta_th': eta_th, 'discrepancy': discrepancy,
        'eta_el': eta_el,
    }


# =============================================================================
# PRINT
# =============================================================================

def print_results(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                  perf, f_rh, P_reheat):

    def _row(label, T_C, p_bar, h):
        print(f"  {label:<36} T={T_C:7.2f} °C   p={p_bar:7.3f} bar   h={h:10.1f} J/kg")

    def _T(h, p):
        return PropsSI('T', 'H', h, 'P', p, 'Water') - 273.15

    print("\n=== Steam states — Passe 1 (reheat) ===")
    _row("Pump:su",              pump.su.T-273.15,    pump.su.p/1e5,    pump.su.h)
    _row("Pump:ex / Eco:su_C",   pump.ex.T-273.15,    pump.ex.p/1e5,    pump.ex.h)
    _row("Eco:ex_C / Eva:su_C",  eco.ex_C.T-273.15,   eco.ex_C.p/1e5,   eco.ex_C.h)
    _row("Eva:ex_C / SH:su_C",   eva.ex_C.T-273.15,   eva.ex_C.p/1e5,   eva.ex_C.h)
    _row("SH:ex_C / Turb_HP:su", sh.ex_C.T-273.15,    sh.ex_C.p/1e5,    sh.ex_C.h)
    _row("Turb_HP:ex / RH:su_C", turb_hp.ex.T-273.15, turb_hp.ex.p/1e5, turb_hp.ex.h)
    _row("RH:ex_C / Turb_LP:su", rh.ex_C.T-273.15,    rh.ex_C.p/1e5,    rh.ex_C.h)
    _row("Turb_LP:ex",           turb_lp.ex.T-273.15, turb_lp.ex.p/1e5, turb_lp.ex.h)
    _row("Condenser:ex_H",       condenser.ex_H.T-273.15, condenser.ex_H.p/1e5, condenser.ex_H.h)

    P_fw = perf['P_fw']
    P_s1 = perf['P_s1']
    P_s2 = perf['P_s2']

    print(f"\n=== FWH states — Passe 2 (analytique) ===")
    _row("Condensat / Pump_LP:ex",
         _T(perf['h_plp_ex'], P_fw),    P_fw/1e5,  perf['h_plp_ex'])
    _row(f"FW out cFWH_LP1 / in cFWH_LP2",
         _T(perf['h_fw_out_lp1'], P_fw), P_fw/1e5, perf['h_fw_out_lp1'])
    _row(f"Soutirage cFWH_LP1 @ {P_s1/1e5:.3f} bar",
         perf['T_ext_s1']-273.15,        P_s1/1e5, perf['h_drain_s1'])
    _row(f"FW out cFWH_LP2 / in DEA",
         _T(perf['h_fw_out_lp2'], P_fw), P_fw/1e5, perf['h_fw_out_lp2'])
    _row(f"Soutirage cFWH_LP2 @ {P_s2/1e5:.3f} bar",
         perf['T_ext_s2']-273.15,        P_s2/1e5, perf['h_drain_s2'])
    _row(f"Soutirage DEA @ {P_fw/1e5:.2f} bar",
         perf['T_ext_fw']-273.15,        P_fw/1e5, perf['h_fw_out'])
    _row("DEA out / Pump_HP:su",
         perf['T_fw_out']-273.15,        P_fw/1e5, perf['h_fw_out'])
    _row("Pump_HP:ex / Eco:su_C",
         _T(perf['h_php_ex'], P_high),   P_high/1e5, perf['h_php_ex'])

    print(f"\n=== Performance — Reheat + 2 cFWH LP + DEA ===")
    print(f"  f_rh       : {f_rh:.2f}   P_reheat : {P_reheat/1e5:.0f} bar")
    print(f"  P_s1       : {P_s1/1e5:.3f} bar   T_sat={PropsSI('T','P',P_s1,'Q',0,'Water')-273.15:.1f} °C")
    print(f"  P_s2       : {P_s2/1e5:.3f} bar   T_sat={PropsSI('T','P',P_s2,'Q',0,'Water')-273.15:.1f} °C")
    print(f"  P_fw (DEA) : {P_fw/1e5:.2f} bar   T_sat={PropsSI('T','P',P_fw,'Q',0,'Water')-273.15:.1f} °C")
    print(f"  y_lp1      : {perf['y_lp1']:.4f}  [-]")
    print(f"  y_lp2      : {perf['y_lp2']:.4f}  [-]")
    print(f"  y_dea      : {perf['y_dea']:.4f}  [-]")
    print(f"  m_lp4      : {perf['m_lp4']:.4f}  [-]")
    print(f"  ---")
    print(f"  W_pump_lp  : {perf['W_pump_lp']:.4f} kW/(kg/s)")
    print(f"  W_pump_hp  : {perf['W_pump_hp']:.4f} kW/(kg/s)")
    print(f"  W_HP       : {perf['W_HP']:.2f} kW/(kg/s)")
    print(f"  W_LP1      : {perf['W_LP1']:.2f} kW/(kg/s)")
    print(f"  W_LP2      : {perf['W_LP2']:.2f} kW/(kg/s)")
    print(f"  W_LP3      : {perf['W_LP3']:.2f} kW/(kg/s)")
    print(f"  W_LP4      : {perf['W_LP4']:.2f} kW/(kg/s)")
    print(f"  W_turb     : {perf['W_turb']:.2f} kW/(kg/s)")
    print(f"  W_net      : {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Q_boiler   : {perf['Q_boiler']:.2f} kW/(kg/s)")
    print(f"  Q_rh       : {perf['Q_rh']:.2f} kW/(kg/s)")
    print(f"  Q_boil     : {perf['Q_boil']:.2f} kW/(kg/s)")
    print(f"  Q_cond     : {perf['Q_cond']:.2f} kW/(kg/s)")
    print(f"  eta_th     : {perf['eta_th']*100:.2f} %")
    print(f"  eta_el     : {perf['eta_el']*100:.2f} %")
    print(f"  --- First law check ---")
    print(f"  Q_boil - Q_cond : {perf['Q_boil'] - perf['Q_cond']:.4f} kW")
    print(f"  W_net           : {perf['W_net']:.4f} kW")
    print(f"  Discrepancy     : {perf['discrepancy']:.6f} kW")
    print(f"  ==========================================")


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":

    # --- Paramètres machine (Siemens SST-PAC 700+900 DCRH) ---
    eta_g = 0.986
    eta_pump   = 0.75
    eta_turb   = 0.85
    eta_hx     = 0.95
    pinch_hx   = 5.0
    pinch_cond = 5.0
    pinch_cfwh = 5.0

    # --- Pressions cycle (lues dans le PDF) ---
    P_high   = 160e5      # bar → Pa
    P_low    = 0.100e5    # 100 mbar (P_exh du PDF)
    P_reheat = 38.9e5     # bar reheat Siemens

    # --- f_rh : à garder optimisé sur notre sel ---
    f_rh = 0.10

    # --- Pressions FWH (lues dans le PDF) ---
    P_fw  = 11.64e5   # DEA @ 11.64 bar
    P_s2  =  3.17e5   # cFWH_LP2 @ 3.17 bar
    P_s1  =  1.38e5   # cFWH_LP1 @ 1.38 bar

    # --- Source sel ---
    m_dot_st       = 1.0
    T_salt_su      = 565 + 273.15
    P_salt         = 2e5
    m_dot_salt_tot = 100 * 5.99

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=20+273.15, P=3e5, m_dot=100.0)

    print("Passe 1 : résolution cycle reheat...")
    pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser = \
        solve_reheat(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                     T_salt_su, P_salt, m_dot_salt_tot, f_rh,
                     CSource, P_low, P_high, P_reheat, m_dot_st)
    print("  OK")

    print("Passe 2 : calcul analytique cFWH...")
    perf = compute_reheat_cfwh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                                eta_pump, eta_turb, P_high, P_low, P_fw, P_s1, P_s2, eta_g,
                                pinch_cfwh=pinch_cfwh)
    print("  OK")

    print_results(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                  perf, f_rh, P_reheat)

    h_sh_ex   = sh.ex_C.h
    h_pump_ex = pump.ex.h
    h_pump_su = pump.su.h
    W_net_base  = ((turb_hp.su.h - turb_hp.ex.h)
                 + (turb_lp.su.h - turb_lp.ex.h)
                 - (h_pump_ex - h_pump_su)) / 1000
    Q_boil_base = ((h_sh_ex - h_pump_ex) + (rh.ex_C.h - turb_hp.ex.h)) / 1000
    eta_base    = W_net_base / Q_boil_base
    eta_base_el = eta_base*eta_g

    print(f"\n=== Comparaison ===")
    print(f"  Reheat seul            : η = {eta_base_el*100:.2f} %  W_net = {W_net_base:.2f} kW/(kg/s)")
    print(f"  Reheat + 2cFWH + DEA   : η = {perf['eta_el']*100:.2f} %  W_net = {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Gain total             : {(perf['eta_el'] - eta_base_el)*100:+.2f} pp")