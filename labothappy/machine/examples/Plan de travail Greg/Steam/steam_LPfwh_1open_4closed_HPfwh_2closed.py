# -*- coding: utf-8 -*-
"""
Created on Apr 2026
@author: gregoire.hendrix

Steam Rankine — Architecture Siemens SST-PAC 700+900 DCRH
Triple reheat : HP + IP_A + IP_B + IP_C + LP

Passe 1 — RecursiveCircuit :
    Pump → Eco → Eva → SH → Turb_HP → RH1 → Turb_IPA → RH2
         → Turb_IPB → Turb_IPC → RH3 → Turb_LP → Condenser

    Salt split :
        f_boiler → SH → Eva → Eco  (series)
        f_rh1    → RH1  (157 bar)
        f_rh2    → RH2  (41.75 bar)
        f_rh3    → RH3  (11.64 bar)

Passe 2 — Analytique :
    Turb_IPA : 157   → 41.75 bar  soutirage cFWH_HP2
    Turb_IPB : 38.9  → 19.70 bar  soutirage cFWH_HP1
    Turb_IPC : 19.70 → 11.64 bar  soutirage DEA
    Turb_LP  : 11.64 → 6.08 → 3.17 → 1.38 → 0.095 bar  4 cFWH LP

    Côté eau :
        Pump_LP → cFWH_LP1..4 → DEA → Pump_HP → cFWH_HP1..2 → Eco
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
# PASSE 1 — RecursiveCircuit triple reheat
# =============================================================================

def _build_circuit(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                   HSource_sh, Salt_eva_su, Salt_eco_su,
                   HSource_rh1, HSource_rh2, HSource_rh3, CSource,
                   P_low, P_high, P_rh1, P_rh2, P_rh3, m_dot_st, n_disc=10):

    T_sat_hi = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_salt   = HSource_sh.T

    cycle     = RecursiveCircuit('Water')
    pump      = PumpCstEff();   pump.set_parameters(eta_is=eta_pump)
    eco       = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    eva       = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    sh        = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    rh1       = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    rh2       = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    rh3       = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    turb_hp   = ExpanderCstEff(); turb_hp.set_parameters(eta_is=eta_turb)
    turb_ipa  = ExpanderCstEff(); turb_ipa.set_parameters(eta_is=eta_turb)
    turb_ipb  = ExpanderCstEff(); turb_ipb.set_parameters(eta_is=eta_turb)
    turb_ipc  = ExpanderCstEff(); turb_ipc.set_parameters(eta_is=eta_turb)
    turb_lp   = ExpanderCstEff(); turb_lp.set_parameters(eta_is=eta_turb)
    cond      = _make_pinch_cond(pinch=pinch_cond)

    for name, comp in [
        ("Pump", pump), ("Economiser", eco), ("Evaporator", eva),
        ("Superheater", sh),
        ("Turb_HP",  turb_hp),  ("Reheater1", rh1),
        ("Turb_IPA", turb_ipa), ("Reheater2", rh2),
        ("Turb_IPB", turb_ipb), ("Turb_IPC",  turb_ipc),
        ("Reheater3", rh3),     ("Turb_LP",   turb_lp),
        ("Condenser", cond),
    ]:
        cycle.add_component(comp, name)

    cycle.link_components("Pump",        "m-ex",   "Economiser",  "m-su_C")
    cycle.link_components("Economiser",  "m-ex_C", "Evaporator",  "m-su_C")
    cycle.link_components("Evaporator",  "m-ex_C", "Superheater", "m-su_C")
    cycle.link_components("Superheater", "m-ex_C", "Turb_HP",     "m-su")
    cycle.link_components("Turb_HP",     "m-ex",   "Reheater1",   "m-su_C")
    cycle.link_components("Reheater1",   "m-ex_C", "Turb_IPA",    "m-su")
    cycle.link_components("Turb_IPA",    "m-ex",   "Reheater2",   "m-su_C")
    cycle.link_components("Reheater2",   "m-ex_C", "Turb_IPB",    "m-su")
    cycle.link_components("Turb_IPB",    "m-ex",   "Turb_IPC",    "m-su")
    cycle.link_components("Turb_IPC",    "m-ex",   "Reheater3",   "m-su_C")
    cycle.link_components("Reheater3",   "m-ex_C", "Turb_LP",     "m-su")
    cycle.link_components("Turb_LP",     "m-ex",   "Condenser",   "m-su_H")
    cycle.link_components("Condenser",   "m-ex_H", "Pump",        "m-su")

    cycle.add_source("SaltSource_SH",  HSource_sh,  cycle.components["Superheater"], "m-su_H")
    cycle.add_source("SaltSource_EVA", Salt_eva_su, cycle.components["Evaporator"],  "m-su_H")
    cycle.add_source("SaltSource_ECO", Salt_eco_su, cycle.components["Economiser"],  "m-su_H")
    cycle.add_source("SaltSource_RH1", HSource_rh1, cycle.components["Reheater1"],   "m-su_H")
    cycle.add_source("SaltSource_RH2", HSource_rh2, cycle.components["Reheater2"],   "m-su_H")
    cycle.add_source("SaltSource_RH3", HSource_rh3, cycle.components["Reheater3"],   "m-su_H")
    cycle.add_source("ColdSource",     CSource,     cycle.components["Condenser"],   "m-su_C")

    cycle.set_fixed_properties(target="Pump:ex",     p=P_high)
    cycle.set_fixed_properties(target="Turb_HP:ex",  p=P_rh1)
    cycle.set_fixed_properties(target="Turb_IPA:ex", p=P_rh2)
    cycle.set_fixed_properties(target="Turb_IPB:ex", p=19.70e5)
    cycle.set_fixed_properties(target="Turb_IPC:ex", p=P_rh3)
    cycle.set_fixed_properties(target="Turb_LP:ex",  p=P_low)

    T_su_g   = T_salt - pinch_hx
    h_sat_lo = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    h_pu_ex  = _pump(h_sat_lo, P_low, P_high, eta_pump)
    T_pu_ex  = PropsSI('T', 'H', h_pu_ex, 'P', P_high, 'Water')

    h_HP_su  = PropsSI('H', 'T', T_su_g, 'P', P_high,  'Water')
    h_HP_ex  = _expand(h_HP_su,  P_high,  P_rh1,   eta_turb)
    T_HP_ex  = PropsSI('T', 'H', h_HP_ex,  'P', P_rh1,   'Water')

    h_IPA_su = PropsSI('H', 'T', T_su_g, 'P', P_rh1,   'Water')
    h_IPA_ex = _expand(h_IPA_su, P_rh1,   P_rh2,   eta_turb)
    T_IPA_ex = PropsSI('T', 'H', h_IPA_ex, 'P', P_rh2,   'Water')

    h_IPB_su = PropsSI('H', 'T', T_su_g, 'P', P_rh2,   'Water')
    h_IPB_ex = _expand(h_IPB_su, P_rh2,   19.70e5, eta_turb)
    T_IPB_ex = PropsSI('T', 'H', h_IPB_ex, 'P', 19.70e5, 'Water')

    h_IPC_su = PropsSI('H', 'T', T_su_g, 'P', 19.70e5, 'Water')
    h_IPC_ex = _expand(h_IPC_su, 19.70e5, P_rh3,   eta_turb)
    T_IPC_ex = PropsSI('T', 'H', h_IPC_ex, 'P', P_rh3,   'Water')

    h_LP_su  = PropsSI('H', 'T', T_su_g, 'P', P_rh3,   'Water')
    h_LP_ex  = _expand(h_LP_su,  P_rh3,   P_low,   eta_turb)
    T_LP_ex  = PropsSI('T', 'H', h_LP_ex,  'P', P_low,   'Water')

    cycle.set_cycle_guess(target='Pump:su',           p=P_low,    SC=3,        m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',           p=P_high,   T=T_pu_ex)
    cycle.set_cycle_guess(target='Economiser:ex_C',   p=P_high,   T=T_sat_hi-5)
    cycle.set_cycle_guess(target='Evaporator:ex_C',   p=P_high,   T=T_sat_hi,  x=1)
    cycle.set_cycle_guess(target='Superheater:ex_C',  p=P_high,   T=T_su_g)
    cycle.set_cycle_guess(target='Turb_HP:su',        p=P_high,   T=T_su_g,    m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_HP:ex',        p=P_rh1,    T=T_HP_ex)
    cycle.set_cycle_guess(target='Reheater1:ex_C',    p=P_rh1,    T=T_su_g)
    cycle.set_cycle_guess(target='Turb_IPA:su',       p=P_rh1,    T=T_su_g,    m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_IPA:ex',       p=P_rh2,    T=T_IPA_ex)
    cycle.set_cycle_guess(target='Reheater2:ex_C',    p=P_rh2,    T=T_su_g)
    cycle.set_cycle_guess(target='Turb_IPB:su',       p=P_rh2,    T=T_su_g,    m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_IPB:ex',       p=19.70e5,  T=T_IPB_ex)
    cycle.set_cycle_guess(target='Turb_IPC:su',       p=19.70e5,  T=T_su_g,    m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_IPC:ex',       p=P_rh3,    T=T_IPC_ex)
    cycle.set_cycle_guess(target='Reheater3:ex_C',    p=P_rh3,    T=T_su_g)
    cycle.set_cycle_guess(target='Turb_LP:su',        p=P_rh3,    T=T_su_g,    m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_LP:ex',        p=P_low,    T=T_LP_ex)
    cycle.set_cycle_guess(target='Condenser:ex_H',    p=P_low,    T=T_sat_lo-2)

    for target, var in [
        ('Turb_HP:ex','h'),  ('Turb_HP:ex','p'),
        ('Turb_IPA:ex','h'), ('Turb_IPA:ex','p'),
        ('Turb_IPB:ex','h'), ('Turb_IPB:ex','p'),
        ('Turb_IPC:ex','h'), ('Turb_IPC:ex','p'),
        ('Turb_LP:ex','h'),  ('Turb_LP:ex','p'),
        ('Pump:ex','h'),     ('Pump:ex','p'),
        ('Superheater:ex_C','h'),  ('Superheater:ex_C','p'),
        ('Reheater1:ex_C','h'),    ('Reheater1:ex_C','p'),
        ('Reheater2:ex_C','h'),    ('Reheater2:ex_C','p'),
        ('Reheater3:ex_C','h'),    ('Reheater3:ex_C','p'),
        ('Evaporator:ex_C','h'),   ('Evaporator:ex_C','p'),
        ('Economiser:ex_C','h'),   ('Economiser:ex_C','p'),
        ('Condenser:ex_H','h'),    ('Condenser:ex_H','p'),
    ]:
        cycle.set_residual_variable(target=target, variable=var, tolerance=1e-6)

    return cycle, pump, eco, eva, sh, rh1, rh2, rh3, \
           turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond


def solve_cycle(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                T_salt_su, P_salt, m_dot_salt_total,
                f_rh1, f_rh2, f_rh3,
                CSource, P_low, P_high, P_rh1, P_rh2, P_rh3, m_dot_st,
                tol_outer=0.1, max_outer=20, n_disc=10):

    T_sat_hi     = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    f_boiler     = 1.0 - f_rh1 - f_rh2 - f_rh3
    m_dot_boiler = m_dot_salt_total * f_boiler
    HSource_rh1  = _salt_connector(T_salt_su, P_salt, m_dot_salt_total * f_rh1)
    HSource_rh2  = _salt_connector(T_salt_su, P_salt, m_dot_salt_total * f_rh2)
    HSource_rh3  = _salt_connector(T_salt_su, P_salt, m_dot_salt_total * f_rh3)

    T_eva_su = T_salt_su - 10.0
    T_eco_su = T_sat_hi  + 5.0

    for _ in range(max_outer):
        HSource_sh  = _salt_connector(T_salt_su, P_salt, m_dot_boiler)
        Salt_eva_su = _salt_connector(T_eva_su,  P_salt, m_dot_boiler)
        Salt_eco_su = _salt_connector(T_eco_su,  P_salt, m_dot_boiler)

        cycle, pump, eco, eva, sh, rh1, rh2, rh3, \
        turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond = \
            _build_circuit(
                eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                HSource_sh, Salt_eva_su, Salt_eco_su,
                HSource_rh1, HSource_rh2, HSource_rh3, CSource,
                P_low, P_high, P_rh1, P_rh2, P_rh3, m_dot_st, n_disc,
            )
        cycle.solve()

        T_eva_su_new = sh.ex_H.T
        T_eco_su_new = eva.ex_H.T
        err = max(abs(T_eva_su_new - T_eva_su), abs(T_eco_su_new - T_eco_su))
        T_eva_su = T_eva_su_new
        T_eco_su = T_eco_su_new
        if err < tol_outer:
            break

    return pump, eco, eva, sh, rh1, rh2, rh3, \
           turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond


# =============================================================================
# PASSE 2 — Analytique
# =============================================================================

def compute_cfwh(pump, eco, eva, sh, rh1, rh2, rh3,
                 turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond,
                 eta_pump, eta_turb, P_high, P_low,
                 P_hp2, P_hp1, P_fw, P_s4, P_s3, P_s2, P_s1,
                 pinch_cfwh=5.0):
    """
    Soutirages :
        Turb_IPA : soutirage cFWH_HP2 @ P_hp2=41.75 bar
        Turb_IPB : soutirage cFWH_HP1 @ P_hp1=19.70 bar  (= sortie IPB)
        Turb_IPC : soutirage DEA       @ P_fw=11.64 bar   (= sortie IPC)
        Turb_LP  : 4 soutirages LP @ P_s4..P_s1
    """

    h_hp_su   = turb_hp.su.h
    p_hp_su   = turb_hp.su.p
    h_ipa_su  = rh1.ex_C.h
    p_ipa_su  = turb_ipa.su.p
    h_ipb_su  = rh2.ex_C.h
    p_ipb_su  = turb_ipb.su.p
    h_ipc_su  = turb_ipc.su.h
    p_ipc_su  = turb_ipc.su.p
    h_lp_su   = rh3.ex_C.h
    p_lp_su   = turb_lp.su.p
    h_cond_ex = cond.ex_H.h
    h_sh_ex   = sh.ex_C.h
    P_rh1     = turb_hp.ex.p
    P_rh2     = turb_ipa.ex.p
    P_ipb_ex  = turb_ipb.ex.p   # = P_hp1 = 19.70 bar
    P_rh3     = turb_ipc.ex.p   # = P_fw  = 11.64 bar

    # Turb_HP (pas de soutirage)
    h_hp_ex  = _expand(h_hp_su,  p_hp_su,  P_rh1,   eta_turb)

    # Turb_IPA : P_rh1 → P_hp2 → P_rh2
    h_ext_hp2 = _expand(h_ipa_su,  p_ipa_su, P_hp2,  eta_turb)
    h_ipa_ex  = _expand(h_ext_hp2, P_hp2,    P_rh2,  eta_turb)
    T_ext_hp2 = PropsSI('T', 'H', h_ext_hp2, 'P', P_hp2, 'Water')

    # Turb_IPB : P_rh2 → P_hp1  (soutirage = sortie)
    h_ipb_ex  = _expand(h_ipb_su, p_ipb_su, P_ipb_ex, eta_turb)
    T_ext_hp1 = PropsSI('T', 'H', h_ipb_ex, 'P', P_ipb_ex, 'Water')

    # Turb_IPC : P_hp1 → P_fw   (soutirage DEA = sortie)
    h_ipc_ex  = _expand(h_ipc_su, p_ipc_su, P_rh3, eta_turb)
    T_ext_fw  = PropsSI('T', 'H', h_ipc_ex, 'P', P_rh3, 'Water')

    # Turb_LP segmentée : P_rh3 → P_s4 → P_s3 → P_s2 → P_s1 → P_low
    h_ext_s4 = _expand(h_lp_su,  p_lp_su, P_s4,  eta_turb)
    h_ext_s3 = _expand(h_ext_s4, P_s4,    P_s3,  eta_turb)
    h_ext_s2 = _expand(h_ext_s3, P_s3,    P_s2,  eta_turb)
    h_ext_s1 = _expand(h_ext_s2, P_s2,    P_s1,  eta_turb)
    h_lp5_ex = _expand(h_ext_s1, P_s1,    P_low, eta_turb)
    T_ext_s4 = PropsSI('T', 'H', h_ext_s4, 'P', P_s4, 'Water')
    T_ext_s3 = PropsSI('T', 'H', h_ext_s3, 'P', P_s3, 'Water')
    T_ext_s2 = PropsSI('T', 'H', h_ext_s2, 'P', P_s2, 'Water')
    T_ext_s1 = PropsSI('T', 'H', h_ext_s1, 'P', P_s1, 'Water')

    # Pump_LP : 1 kg/s, P_low → P_fw
    h_plp_ex = _pump(h_cond_ex, P_low, P_fw, eta_pump)

    # Drains sat liq
    h_drain_hp2 = PropsSI('H', 'P', P_hp2,    'Q', 0, 'Water')
    h_drain_hp1 = PropsSI('H', 'P', P_ipb_ex, 'Q', 0, 'Water')
    h_fw_out    = PropsSI('H', 'P', P_fw,      'Q', 0, 'Water')
    T_fw_out    = PropsSI('T', 'P', P_fw,      'Q', 0, 'Water')
    h_drain_s4  = PropsSI('H', 'P', P_s4,      'Q', 0, 'Water')
    h_drain_s3  = PropsSI('H', 'P', P_s3,      'Q', 0, 'Water')
    h_drain_s2  = PropsSI('H', 'P', P_s2,      'Q', 0, 'Water')
    h_drain_s1  = PropsSI('H', 'P', P_s1,      'Q', 0, 'Water')

    # Sorties eau côté tube cFWH
    T_fw_out_lp1 = PropsSI('T', 'P', P_s1,      'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_lp2 = PropsSI('T', 'P', P_s2,      'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_lp3 = PropsSI('T', 'P', P_s3,      'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_lp4 = PropsSI('T', 'P', P_s4,      'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_hp1 = PropsSI('T', 'P', P_ipb_ex,  'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_hp2 = PropsSI('T', 'P', P_hp2,     'Q', 0, 'Water') - pinch_cfwh

    h_fw_out_lp1 = PropsSI('H', 'T', T_fw_out_lp1, 'P', P_fw,   'Water')
    h_fw_out_lp2 = PropsSI('H', 'T', T_fw_out_lp2, 'P', P_fw,   'Water')
    h_fw_out_lp3 = PropsSI('H', 'T', T_fw_out_lp3, 'P', P_fw,   'Water')
    h_fw_out_lp4 = PropsSI('H', 'T', T_fw_out_lp4, 'P', P_fw,   'Water')
    h_fw_out_hp1 = PropsSI('H', 'T', T_fw_out_hp1, 'P', P_high, 'Water')
    h_fw_out_hp2 = PropsSI('H', 'T', T_fw_out_hp2, 'P', P_high, 'Water')

    # cFWH LP
    y_lp1 = (h_fw_out_lp1 - h_plp_ex)     / (h_ext_s1 - h_drain_s1)
    y_lp2 = (h_fw_out_lp2 - h_fw_out_lp1) / (h_ext_s2 - h_drain_s2)
    y_lp3 = (h_fw_out_lp3 - h_fw_out_lp2) / (h_ext_s3 - h_drain_s3)
    y_lp4 = (h_fw_out_lp4 - h_fw_out_lp3) / (h_ext_s4 - h_drain_s4)

    # DEA + cFWH HP : boucle itérative (drains HP → DEA)
    y_hp1 = 0.0
    y_hp2 = 0.0
    for _ in range(20):
        y_dea = (h_fw_out
                 - h_fw_out_lp4
                 - y_hp1 * (h_drain_hp1 - h_fw_out)
                 - y_hp2 * (h_drain_hp2 - h_fw_out)) / (h_ipc_ex - h_fw_out)

        h_php_ex  = _pump(h_fw_out, P_fw, P_high, eta_pump)
        y_hp1_new = (h_fw_out_hp1 - h_php_ex)       / (h_ipb_ex  - h_drain_hp1)
        y_hp2_new = (h_fw_out_hp2 - h_fw_out_hp1)   / (h_ext_hp2 - h_drain_hp2)

        if abs(y_hp1_new - y_hp1) < 1e-8 and abs(y_hp2_new - y_hp2) < 1e-8:
            y_hp1, y_hp2 = y_hp1_new, y_hp2_new
            break
        y_hp1, y_hp2 = y_hp1_new, y_hp2_new

    y_hp  = y_hp1 + y_hp2
    y_lp  = y_lp1 + y_lp2 + y_lp3 + y_lp4
    m_lp5 = 1.0 - y_hp - y_dea - y_lp

    for label, val in [('y_dea', y_dea), ('y_hp1', y_hp1), ('y_hp2', y_hp2),
                       ('y_lp1', y_lp1), ('y_lp2', y_lp2), ('y_lp3', y_lp3),
                       ('y_lp4', y_lp4), ('m_lp5', m_lp5)]:
        if val < -1e-4:
            print(f"  WARNING : {label}={val:.4f} < 0")

    # Travaux [kJ/kg_total]
    W_pump_lp = (h_plp_ex - h_cond_ex)                              / 1000
    W_pump_hp = (h_php_ex - h_fw_out)                               / 1000
    W_HP      = (h_hp_su  - h_hp_ex)                                / 1000
    W_IPA1    = 1.0           * (h_ipa_su  - h_ext_hp2)             / 1000
    W_IPA2    = (1 - y_hp2)   * (h_ext_hp2 - h_ipa_ex)              / 1000
    W_IPB     = (1 - y_hp2 - y_hp1) * (h_ipb_su - h_ipb_ex)        / 1000
    W_IPC     = (1 - y_hp)   * (h_ipc_su  - h_ipc_ex)              / 1000
    W_LP1     = (1 - y_hp)   * (h_lp_su   - h_ext_s4)              / 1000
    W_LP2     = (1 - y_hp - y_dea)           * (h_ext_s4 - h_ext_s3) / 1000
    W_LP3     = (1 - y_hp - y_dea - y_lp4)   * (h_ext_s3 - h_ext_s2) / 1000
    W_LP4     = (1 - y_hp - y_dea - y_lp4
                   - y_lp3)                   * (h_ext_s2 - h_ext_s1) / 1000
    W_LP5     = m_lp5                         * (h_ext_s1 - h_lp5_ex) / 1000
    W_turb    = W_HP + W_IPA1 + W_IPA2 + W_IPB + W_IPC + W_LP1 + W_LP2 + W_LP3 + W_LP4 + W_LP5
    W_net     = W_turb - W_pump_lp - W_pump_hp

    # Chaleurs [kJ/kg_total]
    Q_boiler = (h_sh_ex    - h_fw_out_hp2)                    / 1000
    Q_rh1    = (rh1.ex_C.h - h_hp_ex)                         / 1000
    Q_rh2    = (1 - y_hp2) * (rh2.ex_C.h - h_ipa_ex)          / 1000
    Q_rh3    = (1 - y_hp)  * (rh3.ex_C.h - h_ipc_ex)          / 1000
    Q_boil   = Q_boiler + Q_rh1 + Q_rh2 + Q_rh3
    Q_cond   = Q_boil + W_pump_lp + W_pump_hp - W_turb

    eta_th      = W_net / Q_boil
    discrepancy = abs(W_net - (Q_boil - Q_cond))

    return {
        'y_hp1': y_hp1, 'y_hp2': y_hp2, 'y_hp': y_hp,
        'y_dea': y_dea,
        'y_lp1': y_lp1, 'y_lp2': y_lp2, 'y_lp3': y_lp3, 'y_lp4': y_lp4,
        'm_lp5': m_lp5,
        'P_hp2': P_hp2, 'P_hp1': P_hp1, 'P_fw': P_fw,
        'P_s4': P_s4, 'P_s3': P_s3, 'P_s2': P_s2, 'P_s1': P_s1,
        'T_fw_out': T_fw_out,
        'T_ext_hp2': T_ext_hp2, 'T_ext_hp1': T_ext_hp1,
        'T_ext_fw': T_ext_fw,
        'T_ext_s4': T_ext_s4, 'T_ext_s3': T_ext_s3,
        'T_ext_s2': T_ext_s2, 'T_ext_s1': T_ext_s1,
        'h_plp_ex': h_plp_ex, 'h_php_ex': h_php_ex,
        'h_fw_out': h_fw_out,
        'h_fw_out_lp1': h_fw_out_lp1, 'h_fw_out_lp2': h_fw_out_lp2,
        'h_fw_out_lp3': h_fw_out_lp3, 'h_fw_out_lp4': h_fw_out_lp4,
        'h_fw_out_hp1': h_fw_out_hp1, 'h_fw_out_hp2': h_fw_out_hp2,
        'h_lp5_ex': h_lp5_ex,
        'h_drain_hp1': h_drain_hp1, 'h_drain_hp2': h_drain_hp2,
        'h_drain_s1': h_drain_s1, 'h_drain_s2': h_drain_s2,
        'h_drain_s3': h_drain_s3, 'h_drain_s4': h_drain_s4,
        'W_pump_lp': W_pump_lp, 'W_pump_hp': W_pump_hp,
        'W_HP': W_HP, 'W_IPA1': W_IPA1, 'W_IPA2': W_IPA2,
        'W_IPB': W_IPB, 'W_IPC': W_IPC,
        'W_LP1': W_LP1, 'W_LP2': W_LP2, 'W_LP3': W_LP3,
        'W_LP4': W_LP4, 'W_LP5': W_LP5,
        'W_turb': W_turb, 'W_net': W_net,
        'Q_boiler': Q_boiler, 'Q_rh1': Q_rh1, 'Q_rh2': Q_rh2, 'Q_rh3': Q_rh3,
        'Q_boil': Q_boil, 'Q_cond': Q_cond,
        'eta_th': eta_th, 'discrepancy': discrepancy,
    }


# =============================================================================
# PRINT
# =============================================================================

def print_results(pump, eco, eva, sh, rh1, rh2, rh3,
                  turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond,
                  perf, f_rh1, f_rh2, f_rh3):

    def _row(label, T_C, p_bar, h):
        print(f"  {label:<46} T={T_C:7.2f} °C   p={p_bar:7.3f} bar   h={h:10.1f} J/kg")

    def _T(h, p):
        return PropsSI('T', 'H', h, 'P', p, 'Water') - 273.15

    P_high = turb_hp.su.p
    P_fw   = perf['P_fw']

    print("\n=== Steam states — Passe 1 ===")
    _row("Pump:su",                  pump.su.T-273.15,     pump.su.p/1e5,     pump.su.h)
    _row("Pump:ex / Eco:su_C",       pump.ex.T-273.15,     pump.ex.p/1e5,     pump.ex.h)
    _row("Eco:ex_C / Eva:su_C",      eco.ex_C.T-273.15,    eco.ex_C.p/1e5,    eco.ex_C.h)
    _row("Eva:ex_C / SH:su_C",       eva.ex_C.T-273.15,    eva.ex_C.p/1e5,    eva.ex_C.h)
    _row("SH:ex_C / Turb_HP:su",     sh.ex_C.T-273.15,     sh.ex_C.p/1e5,     sh.ex_C.h)
    _row("Turb_HP:ex / RH1:su_C",    turb_hp.ex.T-273.15,  turb_hp.ex.p/1e5,  turb_hp.ex.h)
    _row("RH1:ex_C / Turb_IPA:su",   rh1.ex_C.T-273.15,    rh1.ex_C.p/1e5,    rh1.ex_C.h)
    _row("Turb_IPA:ex / RH2:su_C",   turb_ipa.ex.T-273.15, turb_ipa.ex.p/1e5, turb_ipa.ex.h)
    _row("RH2:ex_C / Turb_IPB:su",   rh2.ex_C.T-273.15,    rh2.ex_C.p/1e5,    rh2.ex_C.h)
    _row("Turb_IPB:ex / Turb_IPC:su",turb_ipb.ex.T-273.15, turb_ipb.ex.p/1e5, turb_ipb.ex.h)
    _row("Turb_IPC:ex / RH3:su_C",   turb_ipc.ex.T-273.15, turb_ipc.ex.p/1e5, turb_ipc.ex.h)
    _row("RH3:ex_C / Turb_LP:su",    rh3.ex_C.T-273.15,    rh3.ex_C.p/1e5,    rh3.ex_C.h)
    _row("Turb_LP:ex",               turb_lp.ex.T-273.15,  turb_lp.ex.p/1e5,  turb_lp.ex.h)
    _row("Condenser:ex_H",           cond.ex_H.T-273.15,   cond.ex_H.p/1e5,   cond.ex_H.h)

    print(f"\n=== FWH states — Passe 2 ===")
    _row("Condensat / Pump_LP:ex",
         _T(perf['h_plp_ex'], P_fw),      P_fw/1e5,   perf['h_plp_ex'])
    _row("FW out cFWH_LP1 / in LP2",
         _T(perf['h_fw_out_lp1'], P_fw),  P_fw/1e5,   perf['h_fw_out_lp1'])
    _row("FW out cFWH_LP2 / in LP3",
         _T(perf['h_fw_out_lp2'], P_fw),  P_fw/1e5,   perf['h_fw_out_lp2'])
    _row("FW out cFWH_LP3 / in LP4",
         _T(perf['h_fw_out_lp3'], P_fw),  P_fw/1e5,   perf['h_fw_out_lp3'])
    _row("FW out cFWH_LP4 / in DEA",
         _T(perf['h_fw_out_lp4'], P_fw),  P_fw/1e5,   perf['h_fw_out_lp4'])
    _row("DEA out / Pump_HP:su",
         perf['T_fw_out']-273.15,         P_fw/1e5,   perf['h_fw_out'])
    _row("Pump_HP:ex / cFWH_HP1:su_C",
         _T(perf['h_php_ex'], P_high),    P_high/1e5, perf['h_php_ex'])
    _row("FW out cFWH_HP1 / in HP2",
         _T(perf['h_fw_out_hp1'], P_high), P_high/1e5, perf['h_fw_out_hp1'])
    _row("FW out cFWH_HP2 / Eco:su_C",
         _T(perf['h_fw_out_hp2'], P_high), P_high/1e5, perf['h_fw_out_hp2'])

    print(f"\n  Soutirages :")
    print(f"  cFWH_LP1 @ {perf['P_s1']/1e5:.3f} bar  T={perf['T_ext_s1']-273.15:.1f}°C  y={perf['y_lp1']:.4f}")
    print(f"  cFWH_LP2 @ {perf['P_s2']/1e5:.3f} bar  T={perf['T_ext_s2']-273.15:.1f}°C  y={perf['y_lp2']:.4f}")
    print(f"  cFWH_LP3 @ {perf['P_s3']/1e5:.3f} bar  T={perf['T_ext_s3']-273.15:.1f}°C  y={perf['y_lp3']:.4f}")
    print(f"  cFWH_LP4 @ {perf['P_s4']/1e5:.3f} bar  T={perf['T_ext_s4']-273.15:.1f}°C  y={perf['y_lp4']:.4f}")
    print(f"  DEA      @ {perf['P_fw']/1e5:.2f} bar   T={perf['T_ext_fw']-273.15:.1f}°C  y={perf['y_dea']:.4f}")
    print(f"  cFWH_HP1 @ {perf['P_hp1']/1e5:.2f} bar  T={perf['T_ext_hp1']-273.15:.1f}°C  y={perf['y_hp1']:.4f}")
    print(f"  cFWH_HP2 @ {perf['P_hp2']/1e5:.2f} bar  T={perf['T_ext_hp2']-273.15:.1f}°C  y={perf['y_hp2']:.4f}")
    print(f"  m_lp5 (→ condenseur) = {perf['m_lp5']:.4f}")

    print(f"\n=== Performance ===")
    print(f"  f_rh1={f_rh1:.3f}  f_rh2={f_rh2:.3f}  f_rh3={f_rh3:.3f}")
    print(f"  ---")
    print(f"  W_pump_lp : {perf['W_pump_lp']:.4f} kW/(kg/s)")
    print(f"  W_pump_hp : {perf['W_pump_hp']:.4f} kW/(kg/s)")
    print(f"  W_HP      : {perf['W_HP']:.2f} kW/(kg/s)")
    print(f"  W_IPA1    : {perf['W_IPA1']:.2f} kW/(kg/s)")
    print(f"  W_IPA2    : {perf['W_IPA2']:.2f} kW/(kg/s)")
    print(f"  W_IPB     : {perf['W_IPB']:.2f} kW/(kg/s)")
    print(f"  W_IPC     : {perf['W_IPC']:.2f} kW/(kg/s)")
    print(f"  W_LP1     : {perf['W_LP1']:.2f} kW/(kg/s)")
    print(f"  W_LP2     : {perf['W_LP2']:.2f} kW/(kg/s)")
    print(f"  W_LP3     : {perf['W_LP3']:.2f} kW/(kg/s)")
    print(f"  W_LP4     : {perf['W_LP4']:.2f} kW/(kg/s)")
    print(f"  W_LP5     : {perf['W_LP5']:.2f} kW/(kg/s)")
    print(f"  W_turb    : {perf['W_turb']:.2f} kW/(kg/s)")
    print(f"  W_net     : {perf['W_net']:.2f} kW/(kg/s)")
    print(f"  Q_boiler  : {perf['Q_boiler']:.2f} kW/(kg/s)")
    print(f"  Q_rh1     : {perf['Q_rh1']:.2f} kW/(kg/s)")
    print(f"  Q_rh2     : {perf['Q_rh2']:.2f} kW/(kg/s)")
    print(f"  Q_rh3     : {perf['Q_rh3']:.2f} kW/(kg/s)")
    print(f"  Q_boil    : {perf['Q_boil']:.2f} kW/(kg/s)")
    print(f"  Q_cond    : {perf['Q_cond']:.2f} kW/(kg/s)")
    print(f"  eta_th    : {perf['eta_th']*100:.2f} %")
    print(f"  --- First law check ---")
    print(f"  Q_boil - Q_cond : {perf['Q_boil'] - perf['Q_cond']:.4f} kW")
    print(f"  W_net           : {perf['W_net']:.4f} kW")
    print(f"  Discrepancy     : {perf['discrepancy']:.6f} kW")
    print(f"  ==========================================")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    eta_pump   = 0.75
    eta_turb   = 0.85
    eta_hx     = 0.95
    pinch_hx   = 5.0
    pinch_cond = 5.0
    pinch_cfwh = 5.0

    P_high = 160e5
    P_low  =  0.100e5
    P_rh1  = 157e5      # sortie Turb_HP
    P_rh2  =  41.75e5   # sortie Turb_IPA
    P_rh3  =  11.64e5   # sortie Turb_IPC = DEA

    f_rh1  = 0.05
    f_rh2  = 0.05
    f_rh3  = 0.05

    P_hp2 = 41.75e5    # = P_rh2, soutirage fin Turb_IPA
    P_hp1 = 19.70e5    # soutirage = sortie Turb_IPB
    P_fw  = 11.64e5    # DEA = sortie Turb_IPC
    P_s4  =  6.08e5
    P_s3  =  3.17e5
    P_s2  =  1.38e5
    P_s1  =  0.20e5

    m_dot_st       = 1.0
    T_salt_su      = 565 + 273.15
    P_salt         = 2e5
    m_dot_salt_tot = 100 * 5.99

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=20+273.15, P=3e5, m_dot=100.0)

    print("Passe 1 : résolution cycle triple reheat...")
    pump, eco, eva, sh, rh1, rh2, rh3, \
    turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond = \
        solve_cycle(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                    T_salt_su, P_salt, m_dot_salt_tot,
                    f_rh1, f_rh2, f_rh3,
                    CSource, P_low, P_high, P_rh1, P_rh2, P_rh3, m_dot_st)
    print("  OK")

    print("Passe 2 : calcul analytique FWH...")
    perf = compute_cfwh(pump, eco, eva, sh, rh1, rh2, rh3,
                        turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond,
                        eta_pump, eta_turb, P_high, P_low,
                        P_hp2, P_hp1, P_fw, P_s4, P_s3, P_s2, P_s1,
                        pinch_cfwh=pinch_cfwh)
    print("  OK")

    print_results(pump, eco, eva, sh, rh1, rh2, rh3,
                  turb_hp, turb_ipa, turb_ipb, turb_ipc, turb_lp, cond,
                  perf, f_rh1, f_rh2, f_rh3)

    h_sh_ex   = sh.ex_C.h
    h_pump_ex = pump.ex.h
    h_pump_su = pump.su.h
    W_net_base = ((turb_hp.su.h  - turb_hp.ex.h)
                + (turb_ipa.su.h - turb_ipa.ex.h)
                + (turb_ipb.su.h - turb_ipb.ex.h)
                + (turb_ipc.su.h - turb_ipc.ex.h)
                + (turb_lp.su.h  - turb_lp.ex.h)
                - (h_pump_ex - h_pump_su)) / 1000
    Q_boil_base = ((h_sh_ex      - h_pump_ex)
                 + (rh1.ex_C.h   - turb_hp.ex.h)
                 + (rh2.ex_C.h   - turb_ipa.ex.h)
                 + (rh3.ex_C.h   - turb_ipc.ex.h)) / 1000
    eta_base    = W_net_base / Q_boil_base

    print(f"\n=== Comparaison ===")
    print(f"  Triple reheat seul         : η = {eta_base*100:.2f} %")
    print(f"  + FWH complets             : η = {perf['eta_th']*100:.2f} %")
    print(f"  Gain FWH                   : {(perf['eta_th'] - eta_base)*100:+.2f} pp")
    print(f"  Siemens référence          : η = 46.28 %")
    print(f"  Gap restant                : {46.28 - perf['eta_th']*100:.2f} pp")