# -*- coding: utf-8 -*-
"""
Created on Apr 2026
@author: gregoire.hendrix

Steam Rankine — Reheat + 1 Open FWH (dégazeur)

Architecture :
    Passe 1 — RecursiveCircuit :
        Pump → Eco → Eva → SH → Turb_HP → RH → Turb_LP → Condenser
        Salt split : (1-f_rh) → SH→Eva→Eco series, f_rh → RH fresh
        Point fixé : f_rh=0.10, P_rh=25 bar (optimum sweep précédent)

    Passe 2 — Analytique :
        Soutirage y en cours Turb_LP à P_fw=1.0 bar
        Pump_LP (P_low→P_fw) + FWH open + Pump_HP (P_fw→P_high)
        y calculé par bilan énergie sur le FWH
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
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


# =============================================================================
# PASSE 1 — RecursiveCircuit reheat
# =============================================================================

def _build_reheat_circuit(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                           HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh, CSource,
                           P_low, P_high, P_reheat, m_dot_st, n_disc=10):

    T_sat_hi = PropsSI('T', 'P', P_high,   'Q', 1, 'Water')
    T_sat_lo = PropsSI('T', 'P', P_low,    'Q', 0, 'Water')
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

    # Guesses
    h_sat_lo     = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    s_pump_su    = PropsSI('S', 'P', P_low,  'Q', 0, 'Water')
    h_pump_ex    = h_sat_lo + (PropsSI('H','P',P_high,'S',s_pump_su,'Water') - h_sat_lo) / eta_pump
    T_pump_ex    = PropsSI('T', 'H', h_pump_ex, 'P', P_high, 'Water')

    T_HP_su_g  = T_salt - pinch_hx
    h_HP_su    = PropsSI('H', 'T', T_HP_su_g, 'P', P_high, 'Water')
    h_HP_ex    = h_HP_su - eta_turb * (h_HP_su - PropsSI('H','P',P_reheat,'S',PropsSI('S','T',T_HP_su_g,'P',P_high,'Water'),'Water'))
    T_HP_ex    = PropsSI('T', 'H', h_HP_ex, 'P', P_reheat, 'Water')

    h_LP_su    = PropsSI('H', 'T', T_HP_su_g, 'P', P_reheat, 'Water')
    h_LP_ex    = h_LP_su - eta_turb * (h_LP_su - PropsSI('H','P',P_low,'S',PropsSI('S','T',T_HP_su_g,'P',P_reheat,'Water'),'Water'))
    T_LP_ex    = PropsSI('T', 'H', h_LP_ex, 'P', P_low, 'Water')

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
    """Résout le cycle reheat avec outer loop salt series."""

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
                HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh, CSource,
                P_low, P_high, P_reheat, m_dot_st, n_disc,
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
# PASSE 2 — Analytique FWH
# =============================================================================

def compute_reheat_fwh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                        eta_pump, eta_turb, P_high, P_low, P_fw):
    """
    Calcul analytique du cycle reheat + open FWH.

    Soutirage y en cours Turb_LP à P_fw.
    Base : 1 kg/s en entrée Turb_HP.
    """

    # --- États depuis passe 1 ---
    h_lp_su   = rh.ex_C.h           # entrée Turb_LP (= sortie RH)
    p_lp_su   = turb_lp.su.p        # P_reheat
    h_cond_ex = condenser.ex_H.h    # sortie condensat
    h_sh_ex   = sh.ex_C.h           # sortie SH = entrée boiler côté eau

    # --- Pump_LP : P_low → P_fw ---
    s_cond    = PropsSI('S', 'H', h_cond_ex, 'P', P_low, 'Water')
    h_plp_ex  = h_cond_ex + (PropsSI('H','P',P_fw,'S',s_cond,'Water') - h_cond_ex) / eta_pump

    # --- Soutirage à P_fw : expansion isentropique depuis entrée Turb_LP ---
    s_lp_su   = PropsSI('S', 'H', h_lp_su, 'P', p_lp_su, 'Water')
    h_ext_is  = PropsSI('H', 'P', P_fw, 'S', s_lp_su, 'Water')
    h_extract = h_lp_su - eta_turb * (h_lp_su - h_ext_is)
    T_extract = PropsSI('T', 'H', h_extract, 'P', P_fw, 'Water')

    # --- FWH : sortie eau saturée liquide à P_fw ---
    h_fw_out  = PropsSI('H', 'P', P_fw, 'Q', 0, 'Water')
    T_fw_out  = PropsSI('T', 'P', P_fw, 'Q', 0, 'Water')

    # --- Fraction soutirée y ---
    # bilan : y*h_extract + (1-y)*h_plp_ex = h_fw_out
    y = (h_fw_out - h_plp_ex) / (h_extract - h_plp_ex)

    # --- Pump_HP : P_fw → P_high ---
    s_fw      = PropsSI('S', 'P', P_fw,   'Q', 0, 'Water')
    h_php_ex  = h_fw_out + (PropsSI('H','P',P_high,'S',s_fw,'Water') - h_fw_out) / eta_pump

    # --- Turb_LP2 : expansion P_fw → P_low sur (1-y) ---
    s_ext     = PropsSI('S', 'H', h_extract, 'P', P_fw,  'Water')
    h_lp2_ex  = h_extract - eta_turb * (h_extract - PropsSI('H','P',P_low,'S',s_ext,'Water'))

    # --- Travaux [kJ/kg_total] ---
    W_pump_lp = (1-y) * (h_plp_ex  - h_cond_ex)              / 1000
    W_pump_hp =         (h_php_ex  - h_fw_out)                / 1000
    W_HP      =         (turb_hp.su.h - turb_hp.ex.h)         / 1000
    W_LP1     =         (h_lp_su   - h_extract)               / 1000
    W_LP2     = (1-y) * (h_extract - h_lp2_ex)                / 1000
    W_turb    = W_HP + W_LP1 + W_LP2
    W_net     = W_turb - W_pump_lp - W_pump_hp

    # --- Chaleurs [kJ/kg_total] ---
    Q_boiler  =         (h_sh_ex        - h_php_ex)           / 1000
    Q_rh = 1.0 * (rh.ex_C.h - turb_hp.ex.h) / 1000
    Q_boil    = Q_boiler + Q_rh
    Q_cond    = (1-y) * (h_lp2_ex - h_cond_ex)                / 1000

    eta_th    = W_net / Q_boil

    return {
        'y': y, 'P_fw': P_fw, 'T_fw_out': T_fw_out,
        'T_extract': T_extract, 'h_extract': h_extract,
        'h_plp_ex': h_plp_ex, 'h_php_ex': h_php_ex,
        'h_fw_out': h_fw_out, 'h_lp2_ex': h_lp2_ex,
        'W_pump_lp': W_pump_lp, 'W_pump_hp': W_pump_hp,
        'W_HP': W_HP, 'W_LP1': W_LP1, 'W_LP2': W_LP2,
        'W_turb': W_turb, 'W_net': W_net,
        'Q_boiler': Q_boiler, 'Q_rh': Q_rh, 'Q_boil': Q_boil,
        'Q_cond': Q_cond, 'eta_th': eta_th,
    }


# =============================================================================
# PRINT
# =============================================================================

def print_results(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                  perf, f_rh, P_reheat):
    """Print états passe 1 + performances passe 2."""

    def _row(label, T_C, p_bar, h):
        print(f"  {label:<26} T={T_C:7.2f} °C   p={p_bar:7.3f} bar   h={h:10.1f} J/kg")

    print("\n=== Steam states (passe 1 — reheat) ===")
    _row("Pump:su",              pump.su.T-273.15,    pump.su.p/1e5,    pump.su.h)
    _row("Pump:ex / Eco:su_C",   pump.ex.T-273.15,    pump.ex.p/1e5,    pump.ex.h)
    _row("Eco:ex_C / Eva:su_C",  eco.ex_C.T-273.15,   eco.ex_C.p/1e5,   eco.ex_C.h)
    _row("Eva:ex_C / SH:su_C",   eva.ex_C.T-273.15,   eva.ex_C.p/1e5,   eva.ex_C.h)
    _row("SH:ex_C / Turb_HP:su", sh.ex_C.T-273.15,    sh.ex_C.p/1e5,    sh.ex_C.h)
    _row("Turb_HP:ex / RH:su_C", turb_hp.ex.T-273.15, turb_hp.ex.p/1e5, turb_hp.ex.h)
    _row("RH:ex_C / Turb_LP:su", rh.ex_C.T-273.15,    rh.ex_C.p/1e5,    rh.ex_C.h)
    _row("Turb_LP:ex",           turb_lp.ex.T-273.15, turb_lp.ex.p/1e5, turb_lp.ex.h)
    _row("Condenser:ex_H",       condenser.ex_H.T-273.15, condenser.ex_H.p/1e5, condenser.ex_H.h)

    print(f"\n=== FWH states (passe 2 — analytique, P_fw={perf['P_fw']/1e5:.2f} bar) ===")
    _row("Condensat → Pump_LP:ex",  perf['h_plp_ex']/1e3*0+
         (PropsSI('T','H',perf['h_plp_ex'],'P',perf['P_fw'],'Water')-273.15),
         perf['P_fw']/1e5, perf['h_plp_ex'])
    _row("Soutirage @ P_fw",        perf['T_extract']-273.15, perf['P_fw']/1e5, perf['h_extract'])
    _row("FWH out / Pump_HP:su",    perf['T_fw_out']-273.15,  perf['P_fw']/1e5, perf['h_fw_out'])
    _row("Pump_HP:ex / Eco:su_C",
         PropsSI('T','H',perf['h_php_ex'],'P',perf['P_fw']*0+160e5,'Water')-273.15,
         160.0, perf['h_php_ex'])

    print(f"\n=== Performance — Reheat + Open FWH ===")
    print(f"  f_rh            : {f_rh:.2f}  [-]")
    print(f"  P_reheat        : {P_reheat/1e5:.0f} bar")
    print(f"  P_fw            : {perf['P_fw']/1e5:.2f} bar  (T_sat={perf['T_fw_out']-273.15:.1f} °C)")
    print(f"  y (soutirage)   : {perf['y']:.4f}  [-]")
    print(f"  ---")
    print(f"  W_pump_lp  : {perf['W_pump_lp']:.4f} kW/(kg/s)")
    print(f"  W_pump_hp  : {perf['W_pump_hp']:.4f} kW/(kg/s)")
    print(f"  W_HP       : {perf['W_HP']:.2f} kW/(kg/s)")
    print(f"  W_LP1      : {perf['W_LP1']:.2f} kW/(kg/s)")
    print(f"  W_LP2      : {perf['W_LP2']:.2f} kW/(kg/s)")
    print(f"  W_turb     : {perf['W_turb']:.2f} kW/(kg/s)")
    print(f"  W_net      : {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Q_boiler   : {perf['Q_boiler']:.2f} kW/(kg/s)")
    print(f"  Q_rh       : {perf['Q_rh']:.2f} kW/(kg/s)")
    print(f"  Q_boil     : {perf['Q_boil']:.2f} kW/(kg/s)")
    print(f"  Q_cond     : {perf['Q_cond']:.2f} kW/(kg/s)")
    print(f"  eta_th     : {perf['eta_th']*100:.2f} %")
    print(f"  --- First law check ---")
    discrepancy = abs(perf['W_net'] - (perf['Q_boil'] - perf['Q_cond']))
    print(f"  Q_boil - Q_cond : {perf['Q_boil'] - perf['Q_cond']:.4f} kW")
    print(f"  W_net           : {perf['W_net']:.4f} kW")
    print(f"  Discrepancy     : {discrepancy:.6f} kW")
    print(f"  ========================================")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # --- Paramètres ---
    eta_pump = 0.75
    eta_turb = 0.85
    eta_hx   = 0.95
    pinch_hx   = 5.0    # [K]
    pinch_cond = 5.0    # [K]

    P_high   = 160e5    # [Pa]
    P_low    = 0.10e5   # [Pa]
    P_reheat =  25e5    # [Pa]  optimum du sweep précédent
    P_fw     =   1.0e5  # [Pa]  FWH open (T_sat ≈ 100°C)
    f_rh     = 0.10     # [-]   optimum du sweep précédent

    m_dot_st       = 1.0        # [kg/s] base spécifique
    T_salt_su      = 565+273.15 # [K]
    P_salt         = 2e5        # [Pa]
    m_dot_salt_tot = 100*5.99   # [kg/s]

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=20+273.15, P=3e5, m_dot=100.0)

    # --- Passe 1 : RecursiveCircuit reheat ---
    print("Passe 1 : résolution cycle reheat...")
    pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser = \
        solve_reheat(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                     T_salt_su, P_salt, m_dot_salt_tot, f_rh,
                     CSource, P_low, P_high, P_reheat, m_dot_st)
    print("  OK")

    # --- Passe 2 : analytique FWH ---
    print("Passe 2 : calcul analytique FWH...")
    perf = compute_reheat_fwh(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                               eta_pump, eta_turb, P_high, P_low, P_fw)
    print("  OK")

    # --- Résultats ---
    print_results(pump, eco, eva, sh, rh, turb_hp, turb_lp, condenser,
                  perf, f_rh, P_reheat)

    # --- Comparaison avec reheat sans FWH ---
    h_lp_ex_no_fwh = turb_lp.ex.h
    h_lp_su_no_fwh = turb_lp.su.h
    h_pump_ex_no_fwh = pump.ex.h
    h_pump_su_no_fwh = pump.su.h
    h_sh_ex = sh.ex_C.h
    Q_rh_no_fwh = rh.ex_C.h - turb_hp.ex.h
    W_net_no_fwh  = ((turb_hp.su.h - turb_hp.ex.h)
                   + (h_lp_su_no_fwh - h_lp_ex_no_fwh)
                   - (h_pump_ex_no_fwh - h_pump_su_no_fwh)) / 1000
    Q_boil_no_fwh = ((h_sh_ex - h_pump_ex_no_fwh) + Q_rh_no_fwh) / 1000
    eta_no_fwh    = W_net_no_fwh / Q_boil_no_fwh

    print(f"\n=== Comparaison ===")
    print(f"  Reheat seul      : η_th = {eta_no_fwh*100:.2f} %  "
          f"W_net = {W_net_no_fwh:.2f} kW/(kg/s)")
    print(f"  Reheat + FWH     : η_th = {perf['eta_th']*100:.2f} %  "
          f"W_net = {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Gain             : {(perf['eta_th'] - eta_no_fwh)*100:+.2f} pp")