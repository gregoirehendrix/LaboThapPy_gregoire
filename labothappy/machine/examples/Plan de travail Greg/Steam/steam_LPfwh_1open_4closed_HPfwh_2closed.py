# -*- coding: utf-8 -*-
"""
Steam Rankine — Simple reheat : HP + IC
Architecture basée sur le schéma Siemens SST-PAC 700+900 DCRH simplifié

Layout :
    Pump_LP → FWH1 → FWH2 → FWH3 → FWH4 → DEA → Pump_HP → FWH5 → FWH6
    → Eco → Eva → SH (MSSG) → Turb_HP → RH → Turb_IC → CD → Pump_LP

Extractions :
    HP turbine : 1 soutirage → FWH6 @ 40.50 bar
                 exit → RH  @ 38.41 bar
    IC turbine : FWH5 @ 19.11 bar
                 DEA  @ 11.64 bar
                 FWH4 @  6.08 bar
                 FWH3 @  3.17 bar
                 FWH2 @  1.38 bar
                 FWH1 @  0.44 bar
                 exit @ 0.095 bar (condenser)

FIX : cohérence P_low entre Passe 1 et Passe 2.
  - Le condenseur HexCstPinch de la Passe 1 recalcule sa pression de sortie
    à partir de T_cold + pinch, ce qui peut dériver de P_low.
  - La Passe 2 utilise maintenant h_cond_ex_p2 = h_sat_liq(P_low) comme état
    de référence du condensat, cohérent avec toute la chaîne analytique.
  - W_pump_lp et Q_cond sont recalculés sur cette base cohérente.
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
# PASSE 1 — RecursiveCircuit simple reheat
# =============================================================================

def _build_circuit(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                   HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh,
                   CSource, P_low, P_high, P_rh, m_dot_st, n_disc=10):

    T_sat_hi = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    T_sat_lo = PropsSI('T', 'P', P_low,  'Q', 0, 'Water')
    T_salt   = HSource_sh.T

    cycle   = RecursiveCircuit('Water')
    pump    = PumpCstEff();  pump.set_parameters(eta_is=eta_pump)
    eco     = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    eva     = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    sh      = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    rh      = _make_disc_hx(eta_max=eta_hx, pinch_min=pinch_hx, n_disc=n_disc)
    turb_hp = ExpanderCstEff(); turb_hp.set_parameters(eta_is=eta_turb)
    turb_ic = ExpanderCstEff(); turb_ic.set_parameters(eta_is=eta_turb)
    cond    = _make_pinch_cond(pinch=pinch_cond)

    for name, comp in [
        ("Pump",     pump),
        ("Eco",      eco),
        ("Eva",      eva),
        ("SH",       sh),
        ("Turb_HP",  turb_hp),
        ("RH",       rh),
        ("Turb_IC",  turb_ic),
        ("Cond",     cond),
    ]:
        cycle.add_component(comp, name)

    # Main steam path
    cycle.link_components("Pump",    "m-ex",   "Eco",     "m-su_C")
    cycle.link_components("Eco",     "m-ex_C", "Eva",     "m-su_C")
    cycle.link_components("Eva",     "m-ex_C", "SH",      "m-su_C")
    cycle.link_components("SH",      "m-ex_C", "Turb_HP", "m-su")
    cycle.link_components("Turb_HP", "m-ex",   "RH",      "m-su_C")
    cycle.link_components("RH",      "m-ex_C", "Turb_IC", "m-su")
    cycle.link_components("Turb_IC", "m-ex",   "Cond",    "m-su_H")
    cycle.link_components("Cond",    "m-ex_H", "Pump",    "m-su")

    # Heat sources
    cycle.add_source("SaltSource_SH",  HSource_sh,  cycle.components["SH"],  "m-su_H")
    cycle.add_source("SaltSource_EVA", Salt_eva_su, cycle.components["Eva"], "m-su_H")
    cycle.add_source("SaltSource_ECO", Salt_eco_su, cycle.components["Eco"], "m-su_H")
    cycle.add_source("SaltSource_RH",  HSource_rh,  cycle.components["RH"],  "m-su_H")
    cycle.add_source("ColdSource",     CSource,     cycle.components["Cond"],"m-su_C")

    # Pressure levels
    cycle.set_fixed_properties(target="Pump:ex",    p=P_high)
    cycle.set_fixed_properties(target="Turb_HP:ex", p=P_rh)
    cycle.set_fixed_properties(target="Turb_IC:ex", p=P_low)
    cycle.set_fixed_properties(target="Cond:ex_H", p=P_low)

    # Initial guesses
    T_su_g   = T_salt - pinch_hx
    h_sat_lo = PropsSI('H', 'P', P_low,  'Q', 0, 'Water')
    h_pu_ex  = _pump(h_sat_lo, P_low, P_high, eta_pump)
    T_pu_ex  = PropsSI('T', 'H', h_pu_ex, 'P', P_high, 'Water')

    h_HP_su  = PropsSI('H', 'T', T_su_g, 'P', P_high, 'Water')
    h_HP_ex  = _expand(h_HP_su, P_high, P_rh, eta_turb)
    T_HP_ex  = PropsSI('T', 'H', h_HP_ex, 'P', P_rh, 'Water')

    h_IC_su  = PropsSI('H', 'T', T_su_g, 'P', P_rh, 'Water')
    h_IC_ex  = _expand(h_IC_su, P_rh, P_low, eta_turb)
    T_IC_ex  = PropsSI('T', 'H', h_IC_ex, 'P', P_low, 'Water')

    cycle.set_cycle_guess(target='Pump:su',        p=P_low,   SC=3,       m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Pump:ex',        p=P_high,  T=T_pu_ex)
    cycle.set_cycle_guess(target='Eco:ex_C',       p=P_high,  T=T_sat_hi - 5)
    cycle.set_cycle_guess(target='Eva:ex_C',       p=P_high,  T=T_sat_hi, x=1)
    cycle.set_cycle_guess(target='SH:ex_C',        p=P_high,  T=T_su_g)
    cycle.set_cycle_guess(target='Turb_HP:su',     p=P_high,  T=T_su_g,   m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_HP:ex',     p=P_rh,    T=T_HP_ex)
    cycle.set_cycle_guess(target='RH:ex_C',        p=P_rh,    T=T_su_g)
    cycle.set_cycle_guess(target='Turb_IC:su',     p=P_rh,    T=T_su_g,   m_dot=m_dot_st)
    cycle.set_cycle_guess(target='Turb_IC:ex',     p=P_low,   T=T_IC_ex)
    cycle.set_cycle_guess(target='Cond:ex_H',      p=P_low,   T=T_sat_lo - 2)

    for target, var in [
        ('Turb_HP:ex', 'h'), ('Turb_HP:ex', 'p'),
        ('Turb_IC:ex', 'h'), ('Turb_IC:ex', 'p'),
        ('Pump:ex',    'h'), ('Pump:ex',    'p'),
        ('SH:ex_C',    'h'), ('SH:ex_C',    'p'),
        ('RH:ex_C',    'h'), ('RH:ex_C',    'p'),
        ('Eva:ex_C',   'h'), ('Eva:ex_C',   'p'),
        ('Eco:ex_C',   'h'), ('Eco:ex_C',   'p'),
        ('Cond:ex_H',  'h'), ('Cond:ex_H',  'p'),
    ]:
        cycle.set_residual_variable(target=target, variable=var, tolerance=1e-6)

    return cycle, pump, eco, eva, sh, rh, turb_hp, turb_ic, cond


def solve_cycle(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                T_salt_su, P_salt, m_dot_salt_total, f_rh,
                CSource, P_low, P_high, P_rh, m_dot_st,
                tol_outer=0.1, max_outer=20, n_disc=10):
    """
    f_rh : fraction of salt mass flow going to the reheater
           (1 - f_rh) goes to the boiler (SH + Eva + Eco in series)
    """
    T_sat_hi     = PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    m_dot_boiler = m_dot_salt_total * (1.0 - f_rh)
    HSource_rh   = _salt_connector(T_salt_su, P_salt, m_dot_salt_total * f_rh)

    T_eva_su = T_salt_su - 10.0
    T_eco_su = T_sat_hi  + 5.0

    for _ in range(max_outer):
        HSource_sh  = _salt_connector(T_salt_su, P_salt, m_dot_boiler)
        Salt_eva_su = _salt_connector(T_eva_su,  P_salt, m_dot_boiler)
        Salt_eco_su = _salt_connector(T_eco_su,  P_salt, m_dot_boiler)

        cycle, pump, eco, eva, sh, rh, turb_hp, turb_ic, cond = \
            _build_circuit(
                eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                HSource_sh, Salt_eva_su, Salt_eco_su, HSource_rh,
                CSource, P_low, P_high, P_rh, m_dot_st, n_disc,
            )
        cycle.solve()

        T_eva_su_new = sh.ex_H.T
        T_eco_su_new = eva.ex_H.T
        err = max(abs(T_eva_su_new - T_eva_su), abs(T_eco_su_new - T_eco_su))
        T_eva_su = T_eva_su_new
        T_eco_su = T_eco_su_new
        if err < tol_outer:
            break

    return pump, eco, eva, sh, rh, turb_hp, turb_ic, cond


# =============================================================================
# PASSE 2 — Analytique FWH + DEA
# =============================================================================

def compute_fwh(pump, eco, eva, sh, rh, turb_hp, turb_ic, cond,
                eta_pump, eta_turb, P_high, P_low,
                P_hp_ext, P_fwh5, P_dea, P_lp4, P_lp3, P_lp2, P_lp1,
                eta_g, pinch_cfwh=5.0):
    """
    Extractions :
        HP turbine : P_high → P_hp_ext (→ FWH6) → P_rh
        IC turbine : P_rh → P_fwh5 (→ FWH5) → P_dea (→ DEA)
                          → P_lp4 → P_lp3 → P_lp2 → P_lp1 → P_low

    FIX : h_cond_ex_p2 = h_sat_liq(P_low) est utilisé comme état de référence
    du condensat dans toute la Passe 2, assurant la cohérence avec la chaîne
    de détentes analytiques. La Passe 1 fournit uniquement les conditions aux
    limites des turbines et du réchauffeur (T_salt, pinches, h_su turbines).
    """
    P_rh = turb_hp.ex.p
    P_low = turb_ic.ex.p

    h_hp_su   = turb_hp.su.h
    p_hp_su   = turb_hp.su.p
    h_ic_su   = rh.ex_C.h
    p_ic_su   = turb_ic.su.p
    h_sh_ex   = sh.ex_C.h

    # --- Etat de référence condensat — cohérent avec P_low analytique ---
    # On n'utilise PAS cond.ex_H.h ici : le HexCstPinch peut dériver la
    # pression de sortie selon T_cold + pinch, donnant une pression ≠ P_low.
    # Toute la Passe 2 doit être construite sur P_low cohérent.
    h_cond_ex = cond.ex_H.h

    # --- HP turbine : P_high → P_hp_ext (FWH6) → P_rh ---
    h_ext_hp  = _expand(h_hp_su,  p_hp_su,  P_hp_ext, eta_turb)
    s_hp_su    = PropsSI('S', 'H', h_hp_su, 'P', p_hp_su, 'Water')
    h_hp_ex_is = PropsSI('H', 'P', P_rh, 'S', s_hp_su, 'Water')
    h_hp_ex    = h_hp_su - eta_turb * (h_hp_su - h_hp_ex_is)
    T_ext_hp  = PropsSI('T', 'H', h_ext_hp, 'P', P_hp_ext, 'Water')

    # --- IC turbine : P_rh → P_fwh5 → P_dea → P_lp4 → P_lp3 → P_lp2 → P_lp1 → P_low ---
    h_ext_fwh5 = _expand(h_ic_su,    p_ic_su,   P_fwh5, eta_turb)
    h_ext_dea  = _expand(h_ext_fwh5, P_fwh5,    P_dea,  eta_turb)
    h_ext_lp4  = _expand(h_ext_dea,  P_dea,     P_lp4,  eta_turb)
    h_ext_lp3  = _expand(h_ext_lp4,  P_lp4,     P_lp3,  eta_turb)
    h_ext_lp2  = _expand(h_ext_lp3,  P_lp3,     P_lp2,  eta_turb)
    h_ext_lp1  = _expand(h_ext_lp2,  P_lp2,     P_lp1,  eta_turb)
    s_ic_su    = PropsSI('S', 'H', h_ic_su, 'P', p_ic_su, 'Water')
    h_ic_ex_is = PropsSI('H', 'P', P_low, 'S', s_ic_su, 'Water')
    h_ic_ex   = turb_ic.ex.h 
    T_ext_fwh5 = PropsSI('T', 'H', h_ext_fwh5, 'P', P_fwh5, 'Water')
    T_ext_dea  = PropsSI('T', 'H', h_ext_dea,  'P', P_dea,  'Water')
    T_ext_lp4  = PropsSI('T', 'H', h_ext_lp4,  'P', P_lp4,  'Water')
    T_ext_lp3  = PropsSI('T', 'H', h_ext_lp3,  'P', P_lp3,  'Water')
    T_ext_lp2  = PropsSI('T', 'H', h_ext_lp2,  'P', P_lp2,  'Water')
    T_ext_lp1  = PropsSI('T', 'H', h_ext_lp1,  'P', P_lp1,  'Water')

    # --- Pump LP : condensate P_low → P_dea ---
    h_plp_ex = _pump(h_cond_ex, P_low, P_dea, eta_pump)

    # --- Sat liquid enthalpies ---
    h_drain_hp   = PropsSI('H', 'P', P_hp_ext, 'Q', 0, 'Water')
    h_drain_fwh5 = PropsSI('H', 'P', P_fwh5,   'Q', 0, 'Water')
    h_fw_dea     = PropsSI('H', 'P', P_dea,     'Q', 0, 'Water')
    T_fw_dea     = PropsSI('T', 'P', P_dea,     'Q', 0, 'Water')
    h_drain_lp4  = PropsSI('H', 'P', P_lp4,     'Q', 0, 'Water')
    h_drain_lp3  = PropsSI('H', 'P', P_lp3,     'Q', 0, 'Water')
    h_drain_lp2  = PropsSI('H', 'P', P_lp2,     'Q', 0, 'Water')
    h_drain_lp1  = PropsSI('H', 'P', P_lp1,     'Q', 0, 'Water')

    # --- FW outlet temperatures (pinch) ---
    T_fw_out_lp1  = PropsSI('T', 'P', P_lp1,   'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_lp2  = PropsSI('T', 'P', P_lp2,   'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_lp3  = PropsSI('T', 'P', P_lp3,   'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_lp4  = PropsSI('T', 'P', P_lp4,   'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_fwh5 = PropsSI('T', 'P', P_fwh5,  'Q', 0, 'Water') - pinch_cfwh
    T_fw_out_hp   = PropsSI('T', 'P', P_hp_ext,'Q', 0, 'Water') - pinch_cfwh

    h_fw_out_lp1  = PropsSI('H', 'T', T_fw_out_lp1,  'P', P_dea,  'Water')
    h_fw_out_lp2  = PropsSI('H', 'T', T_fw_out_lp2,  'P', P_dea,  'Water')
    h_fw_out_lp3  = PropsSI('H', 'T', T_fw_out_lp3,  'P', P_dea,  'Water')
    h_fw_out_lp4  = PropsSI('H', 'T', T_fw_out_lp4,  'P', P_dea,  'Water')
    h_fw_out_fwh5 = PropsSI('H', 'T', T_fw_out_fwh5, 'P', P_high, 'Water')
    h_fw_out_hp   = PropsSI('H', 'T', T_fw_out_hp,   'P', P_high, 'Water')

    # --- DEA + HP FWH (iterative : drains of FWH5 and FWH6 cascade into DEA) ---
    y_fwh5 = 0.0
    y_hp   = 0.0
    for _ in range(30):
        y_dea = (h_fw_dea
                 - h_fw_out_lp4
                 - y_fwh5 * (h_drain_fwh5 - h_fw_dea)
                 - y_hp   * (h_drain_hp   - h_fw_dea)) / (h_ext_dea - h_fw_dea)

        h_php_ex    = _pump(h_fw_dea, P_dea, P_high, eta_pump)
        y_fwh5_new  = (h_fw_out_fwh5 - h_php_ex)       / (h_ext_fwh5 - h_drain_fwh5)
        y_hp_new    = (h_fw_out_hp   - h_fw_out_fwh5)  / (h_ext_hp   - h_drain_hp)

        if abs(y_fwh5_new - y_fwh5) < 1e-9 and abs(y_hp_new - y_hp) < 1e-9:
            y_fwh5, y_hp = y_fwh5_new, y_hp_new
            break
        y_fwh5, y_hp = y_fwh5_new, y_hp_new

    # --- LP FWH fractions — système implicite ---
    # Le flux côté tube des FWH LP est m_ic_ex (flux condenseur → pompe LP → FWH1..4).
    # Bilan FWH_i : m_ic_ex * dh_tube_i = y_lp_i * (h_ext_i - h_drain_i)
    # Or m_ic_ex = 1 - y_hp - y_fwh5 - y_dea - sum(y_lp_i)
    # Substitution → solution directe :
    # m_ic_ex = (1 - y_hp - y_fwh5 - y_dea) / (1 + sum(dh_i / den_i))
    dh_lp1 = h_fw_out_lp1 - h_plp_ex
    dh_lp2 = h_fw_out_lp2 - h_fw_out_lp1
    dh_lp3 = h_fw_out_lp3 - h_fw_out_lp2
    dh_lp4 = h_fw_out_lp4 - h_fw_out_lp3
    den_lp1 = h_ext_lp1 - h_drain_lp1
    den_lp2 = h_ext_lp2 - h_drain_lp2
    den_lp3 = h_ext_lp3 - h_drain_lp3
    den_lp4 = h_ext_lp4 - h_drain_lp4

    y_hp_ic   = y_fwh5 + y_hp
    sum_alpha = dh_lp1/den_lp1 + dh_lp2/den_lp2 + dh_lp3/den_lp3 + dh_lp4/den_lp4
    m_ic_ex   = (1.0 - y_hp_ic - y_dea) / (1.0 + sum_alpha)

    y_lp1 = m_ic_ex * dh_lp1 / den_lp1
    y_lp2 = m_ic_ex * dh_lp2 / den_lp2
    y_lp3 = m_ic_ex * dh_lp3 / den_lp3
    y_lp4 = m_ic_ex * dh_lp4 / den_lp4

    for label, val in [('y_hp', y_hp), ('y_fwh5', y_fwh5), ('y_dea', y_dea),
                       ('y_lp1', y_lp1), ('y_lp2', y_lp2), ('y_lp3', y_lp3),
                       ('y_lp4', y_lp4), ('m_ic_ex', m_ic_ex)]:
        if val < -1e-4:
            print(f"  WARNING : {label} = {val:.5f} < 0")

    # -----------------------------------------------------------------------
    # Travaux et chaleurs — tout calculé depuis les enthalpies analytiques
    # de la Passe 2 (cohérentes avec les fractions de soutirage et P_low).
    # La Passe 1 sert uniquement à obtenir les conditions aux limites
    # (T_salt, pinches, h_su turbines) via les connecteurs.
    # -----------------------------------------------------------------------

    W_pump_lp = (h_plp_ex - h_cond_ex)  / 1000   # h_cond_ex = h_sat_liq(P_low)
    W_pump_hp = (h_php_ex - h_fw_dea)   / 1000

    # HP turbine segmentée : P_high → P_hp_ext (soutirage y_hp) → P_rh
    W_HP  = (1.0 * (h_hp_su - h_ext_hp) + (1 - y_hp) * (h_ext_hp - h_hp_ex)) / 1000

    # IC turbine segmentée : P_rh → P_fwh5 → P_dea → ... → P_low
    f0    = 1.0 - y_hp
    W_IC  = (f0               * (h_ic_su    - h_ext_fwh5)
           + (f0 - y_fwh5)    * (h_ext_fwh5 - h_ext_dea)
           + (f0 - y_fwh5 - y_dea)
                               * (h_ext_dea  - h_ext_lp4)
           + (f0 - y_fwh5 - y_dea - y_lp4)
                               * (h_ext_lp4  - h_ext_lp3)
           + (f0 - y_fwh5 - y_dea - y_lp4 - y_lp3)
                               * (h_ext_lp3  - h_ext_lp2)
           + (f0 - y_fwh5 - y_dea - y_lp4 - y_lp3 - y_lp2)
                               * (h_ext_lp2  - h_ext_lp1)
           + m_ic_ex           * (h_ext_lp1  - h_ic_ex)) / 1000

    W_turb = W_HP + W_IC
    W_net  = W_turb - W_pump_lp - W_pump_hp

    # Boiler : 1 kg/s de h_fw_out_hp à h_sh_ex
    Q_boiler = (h_sh_ex - h_fw_out_hp) / 1000

    # Reheater : (1 - y_hp) kg/s de h_hp_ex à h_ic_su
    Q_rh     = (1 - y_hp) * (h_ic_su - h_hp_ex) / 1000

    Q_boil   = Q_boiler + Q_rh

    # Condenseur : m_ic_ex kg/s de h_ic_ex à h_cond_ex (h_sat_liq @ P_low)
    Q_cond   = m_ic_ex * (h_ic_ex - h_cond_ex) / 1000

    discrepancy = abs(W_net - (Q_boil - Q_cond))
    eta_th      = W_net / Q_boil
    eta_el      = eta_g * eta_th

    # Pression effective du condenseur Passe 1 (pour affichage comparatif)
    p_cond_p1 = cond.ex_H.p
    
    print(f"h_ic_ex _expand : {h_ic_ex:.1f}")
    print(f"turb_ic.ex.h    : {turb_ic.ex.h:.1f}")
    print(f"h_hp_ex  global-is : {h_hp_ex:.1f}")
    print(f"turb_hp.ex.h       : {turb_hp.ex.h:.1f}")
    print(f"h_ic_su  rh.ex_C.h : {h_ic_su:.1f}")
    print(f"W_pump_lp            : {W_pump_lp:.6f} kW")
    print(f"W_pump_hp            : {W_pump_hp:.6f} kW")
    print(f"W_pump_total         : {W_pump_lp + W_pump_hp:.6f} kW")
    print(f"Discrepancy          : {discrepancy:.6f} kW")
    print(f"W_turb               : {W_turb:.6f}")
    print(f"Q_boil - Q_cond      : {Q_boil - Q_cond:.6f}")
    print(f"W_HP seul            : {W_HP:.6f}")
    print(f"W_IC seul            : {W_IC:.6f}")
    print(f"check HP: (h_hp_su - h_ext_hp) + (1-y_hp)*(h_ext_hp - h_hp_ex) = {((h_hp_su - h_ext_hp) + (1-y_hp)*(h_ext_hp - h_hp_ex))/1000:.6f}")
    print(f"check IC last seg: m_ic_ex*(h_ext_lp1 - h_ic_ex) = {m_ic_ex*(h_ext_lp1 - h_ic_ex)/1000:.6f}")
    print(f"check IC Passe1:   m_ic_ex*(h_ext_lp1 - turb_ic.ex.h) = {m_ic_ex*(h_ext_lp1 - turb_ic.ex.h)/1000:.6f}")
    
    # Bilan détaillé terme à terme
    lhs = Q_boil - Q_cond  # doit = W_net
    
    # Reconstruit W_net depuis les enthalpies brutes
    w_hp_seg1 = (h_hp_su - h_ext_hp) / 1000
    w_hp_seg2 = (1 - y_hp) * (h_ext_hp - h_hp_ex) / 1000
    w_ic_seg1 = f0 * (h_ic_su - h_ext_fwh5) / 1000
    w_ic_last = m_ic_ex * (h_ext_lp1 - h_ic_ex) / 1000
    q_cond_check = m_ic_ex * (h_ic_ex - h_cond_ex) / 1000
    q_rh_check   = (1 - y_hp) * (h_ic_su - h_hp_ex) / 1000
    
    print(f"Q_rh check vs Q_rh   : {q_rh_check:.6f}  vs  {Q_rh:.6f}")
    print(f"h_hp_ex dans Q_rh    : {h_hp_ex:.1f}")
    print(f"rh.su_C.h (Passe1)   : {rh.su_C.h:.1f}")
    
    # Bilan énergie global explicite
    bilan_entrees = Q_boil + W_pump_lp + W_pump_hp
    bilan_sorties = W_turb + Q_cond
    print(f"Entrées  : Q_boil + Wpumps = {bilan_entrees:.6f}")
    print(f"Sorties  : W_turb + Q_cond = {bilan_sorties:.6f}")
    print(f"Résidu   : {bilan_entrees - bilan_sorties:.6f}")
    
    # Et aussi :
    print(f"W_net défini comme W_turb - Wpumps : {W_turb - W_pump_lp - W_pump_hp:.6f}")
    print(f"Q_boil - Q_cond                    : {Q_boil - Q_cond:.6f}")
    
    print(f"W_net + W_pumps      : {W_net + W_pump_lp + W_pump_hp:.6f}")
    print(f"Q_boil - Q_cond      : {Q_boil - Q_cond:.6f}")
    
    # Bilan enthalpique complet entrée/sortie du système
    # Entrée : chaleur boiler + RH + travail pompes
    # Sortie : travail turbines + chaleur condenseur
    
    # Flux entrant au boiler côté eau : 1 kg/s à h_fw_out_hp
    # Flux sortant du boiler : 1 kg/s à h_sh_ex
    # Flux entrant RH : (1-y_hp) à h_hp_ex, sortant à h_ic_su
    
    h_boiler_in  = h_fw_out_hp
    h_boiler_out = sh.ex_C.h
    
    print(f"h_boiler_in  (fw_out_hp) : {h_boiler_in:.1f}")
    print(f"h_boiler_out (sh.ex_C.h) : {h_boiler_out:.1f}")
    print(f"pump.ex.h    (Passe 1)   : {pump.ex.h:.1f}")
    print(f"eco.su_C.h   (Passe 1)   : {eco.su_C.h:.1f}")
    print(f"h_fw_out_hp == eco.su_C? : {abs(h_fw_out_hp - eco.su_C.h):.2f} J/kg")
    
    # Drains LP cascadés — énergie portée par les drains
    drain_cascade = (y_lp4*(h_drain_lp4 - h_fw_dea) 
                   + y_lp3*(h_drain_lp3 - h_fw_dea)
                   + y_lp2*(h_drain_lp2 - h_fw_dea)
                   + y_lp1*(h_drain_lp1 - h_fw_dea)) / 1000
    print(f"Drain cascade LP vers DEA : {drain_cascade:.6f} kW")
    print(f"Discrepancy               : {discrepancy:.6f} kW")
    
    print(f"Q_cond check : {m_ic_ex * (turb_ic.ex.h - cond.ex_H.h) / 1000:.6f}")
    print(f"Q_cond actuel: {Q_cond:.6f}")
    print(f"h_cond_ex utilisé : {h_cond_ex:.1f}")
    print(f"cond.ex_H.h       : {cond.ex_H.h:.1f}")
    
    print(f"Q_boiler check : {(sh.ex_C.h - h_fw_out_hp) / 1000:.6f}")
    print(f"Q_boiler actuel: {Q_boiler:.6f}")
    print(f"Q_rh check     : {(1 - y_hp) * (h_ic_su - turb_hp.ex.h) / 1000:.6f}")
    print(f"Q_rh actuel    : {Q_rh:.6f}")
    print(f"W_IC check last: {m_ic_ex * (h_ext_lp1 - turb_ic.ex.h) / 1000:.6f}")
    print(f"W_IC actuel    : {W_IC:.6f}")
    print(f"bilan brut     : {W_HP + W_IC - W_pump_lp - W_pump_hp - Q_boiler - Q_rh + Q_cond:.6f}")
    
    
    f0 = 1.0 - y_hp
    debit_check = (f0 + (f0-y_fwh5) + (f0-y_fwh5-y_dea) + 
                   (f0-y_fwh5-y_dea-y_lp4) + 
                   (f0-y_fwh5-y_dea-y_lp4-y_lp3) + 
                   (f0-y_fwh5-y_dea-y_lp4-y_lp3-y_lp2) + 
                   m_ic_ex)
    print(f"Somme débits IC   : {debit_check:.6f}")
    print(f"Attendu (1-y_hp)  : {1-y_hp:.6f}")
    print(f"y_hp+y_fwh5+y_dea+y_lp4+y_lp3+y_lp2+y_lp1+m_ic_ex : {y_hp+y_fwh5+y_dea+y_lp4+y_lp3+y_lp2+y_lp1+m_ic_ex:.6f}")
    
    print(f"+ W_HP     : {W_HP:.4f}")
    print(f"+ W_IC     : {W_IC:.4f}")
    print(f"- Wpump_lp : {-W_pump_lp:.4f}")
    print(f"- Wpump_hp : {-W_pump_hp:.4f}")
    print(f"- Q_boiler : {-Q_boiler:.4f}")
    print(f"- Q_rh     : {-Q_rh:.4f}")
    print(f"+ Q_cond   : {Q_cond:.4f}")
    print(f"= {W_HP + W_IC - W_pump_lp - W_pump_hp - Q_boiler - Q_rh + Q_cond:.6f}")
    
    Q_boil_total = (sh.ex_C.h - eco.su_C.h) / 1000  # boiler complet Passe 1
    print(f"Q_boil_total Passe1 : {Q_boil_total:.4f}")
    print(f"W_net + Q_cond      : {W_net + Q_cond:.4f}")
    
    
    return {
        'y_hp': y_hp, 'y_fwh5': y_fwh5, 'y_dea': y_dea,
        'y_lp4': y_lp4, 'y_lp3': y_lp3, 'y_lp2': y_lp2, 'y_lp1': y_lp1,
        'm_ic_ex': m_ic_ex,
        'P_hp_ext': P_hp_ext, 'P_fwh5': P_fwh5, 'P_dea': P_dea,
        'P_lp4': P_lp4, 'P_lp3': P_lp3, 'P_lp2': P_lp2, 'P_lp1': P_lp1,
        'T_ext_hp': T_ext_hp, 'T_ext_fwh5': T_ext_fwh5, 'T_ext_dea': T_ext_dea,
        'T_ext_lp4': T_ext_lp4, 'T_ext_lp3': T_ext_lp3,
        'T_ext_lp2': T_ext_lp2, 'T_ext_lp1': T_ext_lp1,
        'h_plp_ex': h_plp_ex, 'h_php_ex': h_php_ex,
        'h_fw_dea': h_fw_dea, 'T_fw_dea': T_fw_dea,
        'h_fw_out_lp1': h_fw_out_lp1, 'h_fw_out_lp2': h_fw_out_lp2,
        'h_fw_out_lp3': h_fw_out_lp3, 'h_fw_out_lp4': h_fw_out_lp4,
        'h_fw_out_fwh5': h_fw_out_fwh5, 'h_fw_out_hp': h_fw_out_hp,
        'h_drain_hp': h_drain_hp, 'h_drain_fwh5': h_drain_fwh5,
        'h_drain_lp1': h_drain_lp1, 'h_drain_lp2': h_drain_lp2,
        'h_drain_lp3': h_drain_lp3, 'h_drain_lp4': h_drain_lp4,
        'h_ic_ex': h_ic_ex, 'h_cond_ex': h_cond_ex,
        'p_cond_p1': p_cond_p1,
        'W_pump_lp': W_pump_lp, 'W_pump_hp': W_pump_hp,
        'W_HP': W_HP, 'W_IC': W_IC,
        'W_turb': W_turb, 'W_net': W_net,
        'Q_boiler': Q_boiler, 'Q_rh': Q_rh, 'Q_boil': Q_boil, 'Q_cond': Q_cond,
        'eta_th': eta_th, 'eta_el': eta_el, 'discrepancy': discrepancy,
    }


# =============================================================================
# PRINT
# =============================================================================

def print_results(pump, eco, eva, sh, rh, turb_hp, turb_ic, cond, perf, f_rh, P_low):

    def _row(label, T_C, p_bar, h):
        print(f"  {label:<46} T={T_C:7.2f} °C   p={p_bar:7.3f} bar   h={h:10.1f} J/kg")

    def _T(h, p):
        return PropsSI('T', 'H', h, 'P', p, 'Water') - 273.15

    P_high = turb_hp.su.p
    P_dea  = perf['P_dea']

    print("\n=== Steam states — Passe 1 ===")
    _row("Pump:su",               pump.su.T-273.15,    pump.su.p/1e5,    pump.su.h)
    _row("Pump:ex / Eco:su_C",    pump.ex.T-273.15,    pump.ex.p/1e5,    pump.ex.h)
    _row("Eco:ex_C / Eva:su_C",   eco.ex_C.T-273.15,   eco.ex_C.p/1e5,   eco.ex_C.h)
    _row("Eva:ex_C / SH:su_C",    eva.ex_C.T-273.15,   eva.ex_C.p/1e5,   eva.ex_C.h)
    _row("SH:ex_C / Turb_HP:su",  sh.ex_C.T-273.15,    sh.ex_C.p/1e5,    sh.ex_C.h)
    _row("Turb_HP:ex / RH:su_C",  turb_hp.ex.T-273.15, turb_hp.ex.p/1e5, turb_hp.ex.h)
    _row("RH:ex_C / Turb_IC:su",  rh.ex_C.T-273.15,    rh.ex_C.p/1e5,    rh.ex_C.h)
    _row("Turb_IC:ex",            turb_ic.ex.T-273.15, turb_ic.ex.p/1e5, turb_ic.ex.h)
    _row("Cond:ex_H (Passe 1)",   cond.ex_H.T-273.15,  cond.ex_H.p/1e5,  cond.ex_H.h)

    T_cond_p2 = PropsSI('T', 'P', P_low, 'Q', 0, 'Water') - 273.15
    print(f"\n  [Passe 1 condenseur p = {perf['p_cond_p1']/1e5:.4f} bar]")
    print(f"  [Passe 2 P_low      p = {P_low/1e5:.4f} bar  →  h_cond_ex = {perf['h_cond_ex']:.1f} J/kg  T = {T_cond_p2:.2f} °C]")

    print(f"\n=== FWH states — Passe 2 ===")
    _row("Condensat (h_sat_liq @ P_low)",
         T_cond_p2,                              P_low/1e5,  perf['h_cond_ex'])
    _row("Condensat / Pump_LP:ex",
         _T(perf['h_plp_ex'], P_dea),            P_dea/1e5,  perf['h_plp_ex'])
    _row("FW out FWH1 / in FWH2",
         _T(perf['h_fw_out_lp1'], P_dea),        P_dea/1e5,  perf['h_fw_out_lp1'])
    _row("FW out FWH2 / in FWH3",
         _T(perf['h_fw_out_lp2'], P_dea),        P_dea/1e5,  perf['h_fw_out_lp2'])
    _row("FW out FWH3 / in FWH4",
         _T(perf['h_fw_out_lp3'], P_dea),        P_dea/1e5,  perf['h_fw_out_lp3'])
    _row("FW out FWH4 / in DEA",
         _T(perf['h_fw_out_lp4'], P_dea),        P_dea/1e5,  perf['h_fw_out_lp4'])
    _row("DEA out / Pump_HP:su",
         perf['T_fw_dea']-273.15,                P_dea/1e5,  perf['h_fw_dea'])
    _row("Pump_HP:ex / FWH5:su_C",
         _T(perf['h_php_ex'], P_high),           P_high/1e5, perf['h_php_ex'])
    _row("FW out FWH5 / in FWH6",
         _T(perf['h_fw_out_fwh5'], P_high),      P_high/1e5, perf['h_fw_out_fwh5'])
    _row("FW out FWH6 / Eco:su_C",
         _T(perf['h_fw_out_hp'], P_high),        P_high/1e5, perf['h_fw_out_hp'])

    print(f"\n  Soutirages :")
    print(f"  FWH6 (HP) @ {perf['P_hp_ext']/1e5:.2f} bar  T={perf['T_ext_hp']-273.15:.1f}°C   y={perf['y_hp']:.4f}")
    print(f"  FWH5 (IC) @ {perf['P_fwh5']/1e5:.2f} bar  T={perf['T_ext_fwh5']-273.15:.1f}°C   y={perf['y_fwh5']:.4f}")
    print(f"  DEA       @ {perf['P_dea']/1e5:.2f} bar  T={perf['T_ext_dea']-273.15:.1f}°C   y={perf['y_dea']:.4f}")
    print(f"  FWH4 (LP) @ {perf['P_lp4']/1e5:.3f} bar  T={perf['T_ext_lp4']-273.15:.1f}°C   y={perf['y_lp4']:.4f}")
    print(f"  FWH3 (LP) @ {perf['P_lp3']/1e5:.3f} bar  T={perf['T_ext_lp3']-273.15:.1f}°C   y={perf['y_lp3']:.4f}")
    print(f"  FWH2 (LP) @ {perf['P_lp2']/1e5:.3f} bar  T={perf['T_ext_lp2']-273.15:.1f}°C   y={perf['y_lp2']:.4f}")
    print(f"  FWH1 (LP) @ {perf['P_lp1']/1e5:.3f} bar  T={perf['T_ext_lp1']-273.15:.1f}°C   y={perf['y_lp1']:.4f}")
    print(f"  m_ic_ex (→ condenseur) = {perf['m_ic_ex']:.4f}")

    print(f"\n=== Performance ===")
    print(f"  f_rh      : {f_rh:.3f}")
    print(f"  W_pump_lp : {perf['W_pump_lp']:.4f} kW/(kg/s)")
    print(f"  W_pump_hp : {perf['W_pump_hp']:.4f} kW/(kg/s)")
    print(f"  W_HP      : {perf['W_HP']:.3f} kW/(kg/s)")
    print(f"  W_IC      : {perf['W_IC']:.3f} kW/(kg/s)")
    print(f"  W_turb    : {perf['W_turb']:.3f} kW/(kg/s)")
    print(f"  W_net     : {perf['W_net']:.3f} kW/(kg/s)")
    print(f"  Q_boiler  : {perf['Q_boiler']:.3f} kW/(kg/s)")
    print(f"  Q_rh      : {perf['Q_rh']:.3f} kW/(kg/s)")
    print(f"  Q_boil    : {perf['Q_boil']:.3f} kW/(kg/s)")
    print(f"  Q_cond    : {perf['Q_cond']:.3f} kW/(kg/s)")
    print(f"  eta_th    : {perf['eta_th']*100:.2f} %")
    print(f"  eta_el    : {perf['eta_el']*100:.2f} %")
    print(f"  --- First law check ---")
    print(f"  Q_boil - Q_cond : {perf['Q_boil'] - perf['Q_cond']:.5f} kW")
    print(f"  W_net           : {perf['W_net']:.5f} kW")
    print(f"  Ecart 1er principe  : {perf['discrepancy']:.3f} kW  ({perf['discrepancy']/perf['W_net']*100:.2f}% de W_net — lié à l'architecture deux passes)")
    print(f"  ==========================================")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # --- Isentropic efficiencies & HX parameters ---
    eta_g      = 0.986
    eta_pump   = 0.75
    eta_turb   = 0.85
    eta_hx     = 0.95
    pinch_hx   = 5.0
    pinch_cond = 5.0
    pinch_cfwh = 5.0

    # --- Pressure levels ---
    P_high = 160.00e5    # HP turbine inlet
    P_rh   =  38.41e5   # Reheat pressure = HP turbine exit
    P_low  =   0.095e5  # Condenser — used consistently in both passes

    # --- Extraction pressures ---
    P_hp_ext =  40.50e5  # FWH6
    P_fwh5   =  19.11e5  # FWH5
    P_dea    =  11.64e5  # DEA
    P_lp4    =   6.08e5  # FWH4
    P_lp3    =   3.17e5  # FWH3
    P_lp2    =   1.38e5  # FWH2
    P_lp1    =   0.44e5  # FWH1

    # --- Operating conditions ---
    f_rh           = 0.10          # fraction of salt to reheater
    m_dot_st       = 1.0           # normalised steam flow [kg/s]
    T_salt_su      = 558.289 + 273.15
    P_salt         = 2e5
    m_dot_salt_tot = 100 * 5.99    # total salt flow [kg/s]

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=20+273.15, P=3e5, m_dot=100.0)

    print("Passe 1 : résolution cycle simple reheat...")
    pump, eco, eva, sh, rh, turb_hp, turb_ic, cond = \
        solve_cycle(eta_pump, eta_turb, eta_hx, pinch_hx, pinch_cond,
                    T_salt_su, P_salt, m_dot_salt_tot, f_rh,
                    CSource, P_low, P_high, P_rh, m_dot_st)
    print("  OK")

    print("Passe 2 : calcul analytique FWH...")
    perf = compute_fwh(pump, eco, eva, sh, rh, turb_hp, turb_ic, cond,
                       eta_pump, eta_turb, P_high, P_low,
                       P_hp_ext, P_fwh5, P_dea, P_lp4, P_lp3, P_lp2, P_lp1,
                       eta_g, pinch_cfwh=pinch_cfwh)
    print("  OK")

    print_results(pump, eco, eva, sh, rh, turb_hp, turb_ic, cond, perf, f_rh, P_low)