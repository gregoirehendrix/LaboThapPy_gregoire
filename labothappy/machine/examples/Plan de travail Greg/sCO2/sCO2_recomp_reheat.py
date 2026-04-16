# -*- coding: utf-8 -*-
"""
sCO2 Recompression + Reheat cycle
Derived from the validated recompression cycle (Hanwha HSC-90-RC-600).

Architecture:
  GasHeater → TurbHP → Reheater → TurbLP
            → RecupHT → RecupLT → Spliter
                 ├─ GasCooler → Compressor ─┐
                 └─ Recompressor ────────────┤
                                             ↓
                                   Mixer → RecupLT → RecupHT → GasHeater

The only structural change vs the base cycle:
  - Expander split into TurbHP (P_high → P_interm) + TurbLP (P_interm → P_low)
  - Reheater inserted between TurbHP and TurbLP (same hot source)
  - P_interm is a free parameter to optimise

@author: gregoire.hendrix
"""

import numpy as np
from scipy.optimize import minimize_scalar
import warnings
warnings.filterwarnings('ignore')

from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.solar_salt_connector import SolarSaltConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc
from labothappy.component.tank.tank_spliter import TankSpliter
from labothappy.component.tank.tank_mixer import TankMixer


# =============================================================================
# CYCLE BUILDER
# =============================================================================
def sCO2_recompression_reheat(
        eta_cp, eta_rc, eta_tbHP, eta_tbLP,
        eta_heater, eta_reheater, eta_cooler,
        eta_recup_lt, eta_recup_ht,
        split_ratio,
        HSource, RHSource, CSource,
        P_low, P_interm, P_high,
        T_c_su, m_dot_CO2,
        T_hot,
        pinch_heater=5.0, pinch_recup=2.0, n_disc=20):

    cycle = RecursiveCircuit('CO2')
    rep_spliter = [split_ratio, 1 - split_ratio]

    compressor   = CompressorCstEff()
    recompressor = CompressorCstEff()
    turbHP       = ExpanderCstEff()
    turbLP       = ExpanderCstEff()
    gas_heater   = HexCstEffDisc()
    reheater     = HexCstEffDisc()
    gas_cooler   = HexCstEff()
    recup_lt     = HexCstEffDisc()
    recup_ht     = HexCstEffDisc()
    spliter      = TankSpliter(outlet_repartition=rep_spliter)
    mixer        = TankMixer(n_inlets=2)

    compressor.set_parameters(eta_is=eta_cp)
    recompressor.set_parameters(eta_is=eta_rc)
    turbHP.set_parameters(eta_is=eta_tbHP)
    turbLP.set_parameters(eta_is=eta_tbLP)
    gas_heater.set_parameters(eta_max=eta_heater,   n_disc=n_disc, Pinch_min=pinch_heater)
    reheater.set_parameters(  eta_max=eta_reheater, n_disc=n_disc, Pinch_min=pinch_heater)
    gas_cooler.set_parameters(eta=eta_cooler)
    recup_lt.set_parameters(eta_max=eta_recup_lt, n_disc=n_disc, Pinch_min=pinch_recup)
    recup_ht.set_parameters(eta_max=eta_recup_ht, n_disc=n_disc, Pinch_min=pinch_recup)

    cycle.add_component(compressor,   "Compressor")
    cycle.add_component(recompressor, "Recompressor")
    cycle.add_component(turbHP,       "TurbHP")
    cycle.add_component(turbLP,       "TurbLP")
    cycle.add_component(gas_heater,   "GasHeater")
    cycle.add_component(reheater,     "Reheater")
    cycle.add_component(gas_cooler,   "GasCooler")
    cycle.add_component(recup_lt,     "RecupLT")
    cycle.add_component(recup_ht,     "RecupHT")
    cycle.add_component(spliter,      "Spliter")
    cycle.add_component(mixer,        "Mixer")

    # -------------------------------------------------------------------------
    # Topology — identical to base recompression except turbine split + reheater
    # -------------------------------------------------------------------------
    # HP expansion → reheat → LP expansion
    cycle.link_components("TurbHP",  "m-ex",   "Reheater", "m-su_C")
    cycle.link_components("Reheater","m-ex_C", "TurbLP",   "m-su")

    # Low pressure side (unchanged vs base)
    cycle.link_components("TurbLP",  "m-ex",   "RecupHT", "m-su_H")
    cycle.link_components("RecupHT", "m-ex_H", "RecupLT", "m-su_H")
    cycle.link_components("RecupLT", "m-ex_H", "Spliter", "m-su")

    # Split branches (unchanged)
    cycle.link_components("Spliter",  "m-ex_1", "GasCooler",   "m-su_H")
    cycle.link_components("GasCooler","m-ex_H", "Compressor",  "m-su")
    cycle.link_components("Compressor","m-ex",  "Mixer",       "m-su_1")

    cycle.link_components("Spliter",     "m-ex_2", "Recompressor", "m-su")
    cycle.link_components("Recompressor","m-ex",   "Mixer",        "m-su_2")

    # High pressure side (unchanged)
    cycle.link_components("Mixer",     "m-ex",   "RecupLT",   "m-su_C")
    cycle.link_components("RecupLT",   "m-ex_C", "RecupHT",   "m-su_C")
    cycle.link_components("RecupHT",   "m-ex_C", "GasHeater", "m-su_C")
    cycle.link_components("GasHeater", "m-ex_C", "TurbHP",    "m-su")

    # External sources
    cycle.add_source("Hot_source",    HSource,  cycle.components["GasHeater"], "m-su_H")
    cycle.add_source("Reheat_source", RHSource, cycle.components["Reheater"],  "m-su_H")
    cycle.add_source("Cold_source",   CSource,  cycle.components["GasCooler"], "m-su_C")

    # -------------------------------------------------------------------------
    # Initial guesses — physically derived from actual pressure ratios
    # Isentropic ideal gas approximation (gamma=1.30 ~ sCO2 supercritical)
    # Valid across full T_hot range (600 to 900 C).
    # -------------------------------------------------------------------------
    gamma    = 1.30
    T_TIT    = T_hot - pinch_heater          # turbine HP inlet [K]
    T_RH_out = T_TIT                         # reheat back to TIT

    # HP turbine exit: P_high -> P_interm
    PR_HP        = P_high / P_interm
    T_turbHP_ex  = T_TIT * (1 - 0.85 * (1 - (1/PR_HP)**((gamma-1)/gamma)))

    # LP turbine exit: P_interm -> P_low
    PR_LP        = P_interm / P_low
    T_turbLP_ex  = T_RH_out * (1 - 0.85 * (1 - (1/PR_LP)**((gamma-1)/gamma)))

    # Recuperator guesses — scale with available temperature span
    T_span         = T_turbLP_ex - T_c_su
    T_recupHT_ex_H = T_c_su + T_span * 0.35
    T_recupHT_ex_C = T_TIT   - (T_TIT - T_c_su) * 0.25

    # TurbHP: P_high → P_interm
    cycle.set_cycle_guess(target='TurbHP:su', m_dot=m_dot_CO2, T=T_TIT,       p=P_high)
    cycle.set_cycle_guess(target='TurbHP:ex', m_dot=m_dot_CO2, T=T_turbHP_ex, p=P_interm)

    # Reheater
    cycle.set_cycle_guess(target='Reheater:su_C', m_dot=m_dot_CO2, T=T_turbHP_ex, p=P_interm)
    cycle.set_cycle_guess(target='Reheater:ex_C', m_dot=m_dot_CO2, T=T_RH_out,    p=P_interm)

    # TurbLP: P_interm → P_low
    cycle.set_cycle_guess(target='TurbLP:su', m_dot=m_dot_CO2, T=T_RH_out,    p=P_interm)
    cycle.set_cycle_guess(target='TurbLP:ex', m_dot=m_dot_CO2, T=T_turbLP_ex, p=P_low)

    # RecupHT — same as base but using T_turbLP_ex
    cycle.set_cycle_guess(target='RecupHT:su_H', m_dot=m_dot_CO2, T=T_turbLP_ex,    p=P_low)
    cycle.set_cycle_guess(target='RecupHT:ex_H', m_dot=m_dot_CO2, T=T_recupHT_ex_H, p=P_low)
    cycle.set_cycle_guess(target='RecupHT:su_C', m_dot=m_dot_CO2, T=T_c_su + 100,   p=P_high)
    cycle.set_cycle_guess(target='RecupHT:ex_C', m_dot=m_dot_CO2, T=T_recupHT_ex_C, p=P_high)

    # RecupLT — identical to base
    cycle.set_cycle_guess(target='RecupLT:su_H', m_dot=m_dot_CO2, T=T_recupHT_ex_H, p=P_low)
    cycle.set_cycle_guess(target='RecupLT:ex_H', m_dot=m_dot_CO2, T=T_c_su + 50,    p=P_low)
    cycle.set_cycle_guess(target='RecupLT:su_C', m_dot=m_dot_CO2, T=T_c_su + 50,    p=P_high)
    cycle.set_cycle_guess(target='RecupLT:ex_C', m_dot=m_dot_CO2, T=T_c_su + 150,   p=P_high)

    # Compressor — identical to base
    cycle.set_cycle_guess(target='Compressor:su',
                          m_dot=m_dot_CO2 * split_ratio, T=T_c_su,      p=P_low)
    cycle.set_cycle_guess(target='Compressor:ex',
                          m_dot=m_dot_CO2 * split_ratio, T=T_c_su + 50, p=P_high)

    # Recompressor — identical to base
    cycle.set_cycle_guess(target='Recompressor:su',
                          m_dot=m_dot_CO2 * (1-split_ratio), T=T_c_su + 60,  p=P_low)
    cycle.set_cycle_guess(target='Recompressor:ex',
                          m_dot=m_dot_CO2 * (1-split_ratio), T=T_c_su + 110, p=P_high)

    # Residuals — same targets as base + new turbine/reheater nodes
    for tgt in ['TurbHP:ex', 'TurbLP:ex',
                'Compressor:ex', 'Recompressor:ex',
                'GasHeater:ex_C', 'Reheater:ex_C', 'GasCooler:ex_H',
                'RecupHT:su_H', 'RecupHT:ex_H', 'RecupLT:ex_H',
                'Spliter:ex_2']:
        cycle.set_residual_variable(target=tgt, variable='h', tolerance=1e-6)
        cycle.set_residual_variable(target=tgt, variable='p', tolerance=1e-6)

    return (cycle, compressor, recompressor, turbHP, turbLP,
            gas_heater, reheater, gas_cooler, recup_lt, recup_ht)


# =============================================================================
# PERFORMANCE
# =============================================================================
def compute_performance(compressor, recompressor, turbHP, turbLP,
                        gas_heater, reheater, gas_cooler, recup_lt, recup_ht):

    m_dot = turbHP.su.m_dot

    W_turbHP = m_dot * (turbHP.su.h - turbHP.ex.h) / 1000
    W_turbLP = m_dot * (turbLP.su.h - turbLP.ex.h) / 1000
    W_comp   = compressor.su.m_dot   * (compressor.ex.h   - compressor.su.h)   / 1000
    W_recomp = recompressor.su.m_dot * (recompressor.ex.h - recompressor.su.h) / 1000

    W_net    = W_turbHP + W_turbLP - W_comp - W_recomp
    Q_heater = m_dot * (gas_heater.ex_C.h - gas_heater.su_C.h) / 1000
    Q_reheat = m_dot * (reheater.ex_C.h   - reheater.su_C.h)   / 1000
    Q_total  = Q_heater + Q_reheat

    return {
        'W_turbHP':    W_turbHP,
        'W_turbLP':    W_turbLP,
        'W_comp':      W_comp,
        'W_recomp':    W_recomp,
        'W_net':       W_net,
        'Q_heater':    Q_heater,
        'Q_reheat':    Q_reheat,
        'Q_total':     Q_total,
        'Q_recup_HT':  m_dot * (recup_ht.ex_C.h - recup_ht.su_C.h) / 1000,
        'Q_recup_LT':  m_dot * (recup_lt.ex_C.h - recup_lt.su_C.h) / 1000,
        'eta':         W_net / Q_total,
        'TIT_HP':      turbHP.su.T - 273.15,
        'TIT_LP':      turbLP.su.T - 273.15,
        'T_turbLP_ex': turbLP.ex.T - 273.15,
        'T_comp_su':   compressor.su.T - 273.15,
        'T_recomp_su': recompressor.su.T - 273.15,
    }


def print_results(perf, P_interm_bar, split_ratio):
    print("\n" + "="*60)
    print("sCO2 RECOMPRESSION + REHEAT — RESULTS")
    print("="*60)
    print(f"  P_interm        = {P_interm_bar:.1f} bar")
    print(f"  split_ratio     = {split_ratio:.3f}")
    print(f"\n  TIT HP          = {perf['TIT_HP']:8.2f} °C")
    print(f"  TIT LP (reheat) = {perf['TIT_LP']:8.2f} °C")
    print(f"  TurbLP exit     = {perf['T_turbLP_ex']:8.2f} °C")
    print(f"  Comp inlet      = {perf['T_comp_su']:8.2f} °C")
    print(f"  Recomp inlet    = {perf['T_recomp_su']:8.2f} °C")
    print(f"\n  W_turbHP        = {perf['W_turbHP']:8.2f} kW")
    print(f"  W_turbLP        = {perf['W_turbLP']:8.2f} kW")
    print(f"  W_comp          = {perf['W_comp']:8.2f} kW")
    print(f"  W_recomp        = {perf['W_recomp']:8.2f} kW")
    print(f"  W_net           = {perf['W_net']:8.2f} kW")
    print(f"\n  Q_heater        = {perf['Q_heater']:8.2f} kW")
    print(f"  Q_reheat        = {perf['Q_reheat']:8.2f} kW")
    print(f"  Q_total         = {perf['Q_total']:8.2f} kW")
    print(f"  Q_recup_HT      = {perf['Q_recup_HT']:8.2f} kW")
    print(f"  Q_recup_LT      = {perf['Q_recup_LT']:8.2f} kW")
    print(f"\n{'─'*40}")
    print(f"  η CYCLE         = {perf['eta']*100:.2f} %")
    print(f"{'─'*40}")


# =============================================================================
# SINGLE-POINT RUNNER
# =============================================================================
def make_hot_sources(T_hot, T_cold):
    """
    Returns (HSource, RHSource, PINCH_HEATER, eta_heater, label)
    - T_hot <= 600°C : Solar Salt (SolarSaltConnector, pinch=1.2°C, eta=0.999)
    - T_hot >  600°C : Air        (MassConnector,       pinch=5.0°C, eta=0.999)
    Two independent sources (infinite reservoir assumption — valid for
    parametric cycle studies, see Crespi 2017).
    """
    T_salt_limit = 600.0 + 273.15

    if T_hot <= T_salt_limit:
        HSource = SolarSaltConnector()
        HSource.set_properties(T=T_hot, p=2e5, m_dot=100)
        RHSource = SolarSaltConnector()
        RHSource.set_properties(T=T_hot, p=2e5, m_dot=100)
        pinch   = 1.2
        eta_hx  = 0.999
        label   = f"Solar Salt at {T_hot-273.15:.1f} °C"
    else:
        HSource = MassConnector()
        HSource.set_properties(fluid='air', T=T_hot, p=1.2e5, m_dot=500)
        RHSource = MassConnector()
        RHSource.set_properties(fluid='air', T=T_hot, p=1.2e5, m_dot=500)
        pinch   = 5.0
        eta_hx  = 0.999
        label   = f"Air at {T_hot-273.15:.1f} °C"

    return HSource, RHSource, pinch, eta_hx, label


def run_one(P_interm_bar, split_ratio,
            P_low, P_high, T_hot, T_cold, m_dot=1.0,
            verbose=False):
    """Returns eta or -1 on failure."""
    try:
        P_interm = P_interm_bar * 1e5

        HSource, RHSource, PINCH_HEATER, eta_heater, _ = make_hot_sources(T_hot, T_cold)
        CSource = MassConnector()
        CSource.set_properties(fluid='air', T=T_cold, p=101325, m_dot=500)

        (cycle, compressor, recompressor, turbHP, turbLP,
         gas_heater, reheater, gas_cooler, recup_lt, recup_ht) = \
            sCO2_recompression_reheat(
                eta_cp=0.80, eta_rc=0.80,
                eta_tbHP=0.90, eta_tbLP=0.90,
                eta_heater=eta_heater, eta_reheater=eta_heater,
                eta_cooler=0.999,
                eta_recup_lt=0.85, eta_recup_ht=0.85,
                split_ratio=split_ratio,
                HSource=HSource, RHSource=RHSource, CSource=CSource,
                P_low=P_low, P_interm=P_interm, P_high=P_high,
                T_c_su=T_cold, m_dot_CO2=m_dot,
                T_hot=T_hot,
                pinch_heater=PINCH_HEATER, pinch_recup=2.0, n_disc=20
            )

        cycle.solve()

        perf = compute_performance(
            compressor, recompressor, turbHP, turbLP,
            gas_heater, reheater, gas_cooler, recup_lt, recup_ht)

        # Sanity checks
        if perf['eta'] < 0.10 or perf['eta'] > 0.80:
            return -1, None
        if perf['W_net'] < 0:
            return -1, None
        if perf['TIT_HP'] < 850:
            return -1, None
        if perf['T_comp_su'] > 60:        # compressor must be near T_cold
            return -1, None

        if verbose:
            print(f"  P_interm={P_interm_bar:.1f} bar | split={split_ratio:.3f} "
                  f"| TIT_LP={perf['TIT_LP']:.1f}°C "
                  f"| TurbLP_ex={perf['T_turbLP_ex']:.1f}°C "
                  f"| η={perf['eta']*100:.2f}%")

        return perf['eta'], perf

    except Exception:
        return -1, None


# =============================================================================
# OPTIMISATION over P_interm (split_ratio fixed from PFD)
# =============================================================================
def optimise(P_low, P_high, T_hot, T_cold,
             split_ratio=47.68/68.11):

    P_low_bar  = P_low  / 1e5
    P_high_bar = P_high / 1e5
    P_geo      = np.sqrt(P_low_bar * P_high_bar)

    _, _, PINCH_HEATER, _, source_label = make_hot_sources(T_hot, T_cold)

    print("="*60)
    print("OPTIMISATION — sCO2 Recompression + Reheat")
    print(f"T_hot={T_hot-273.15:.0f}°C | Hot source: {source_label}")
    print(f"P=[{P_low_bar:.1f}, {P_high_bar:.1f}] bar | pinch_heater={PINCH_HEATER}°C")
    print(f"split_ratio = {split_ratio:.3f} (from Hanwha PFD)")
    print(f"P_geo = {P_geo:.1f} bar (starting point)")
    print("="*60)

    best_eta    = -1
    best_Pinterm = P_geo
    best_perf   = None

    # ------------------------------------------------------------------
    # Step 1 : coarse scan over P_interm
    # ------------------------------------------------------------------
    print("\n[Step 1] Coarse scan over P_interm...")
    P_range = np.linspace(P_low_bar + 15, P_high_bar - 15, 12)

    for P_i in P_range:
        eta, perf = run_one(P_i, split_ratio, P_low, P_high, T_hot, T_cold, verbose=True)
        if eta > best_eta:
            best_eta     = eta
            best_Pinterm = P_i
            best_perf    = perf

    if best_perf is None:
        print("\n✘ No convergence on grid — check your LaboThapPy environment")
        return

    print(f"\n  → Best on grid: P_interm={best_Pinterm:.1f} bar | η={best_eta*100:.2f}%")

    # ------------------------------------------------------------------
    # Step 2 : fine scan in ±20 bar window around best
    # ------------------------------------------------------------------
    print("\n[Step 2] Fine scan ±20 bar around best point...")
    P_fine = np.linspace(max(P_low_bar+5, best_Pinterm-20),
                         min(P_high_bar-5, best_Pinterm+20), 15)

    for P_i in P_fine:
        eta, perf = run_one(P_i, split_ratio, P_low, P_high, T_hot, T_cold, verbose=True)
        if eta > best_eta:
            best_eta     = eta
            best_Pinterm = P_i
            best_perf    = perf

    print(f"\n  → Best after fine scan: P_interm={best_Pinterm:.1f} bar | η={best_eta*100:.2f}%")

    # ------------------------------------------------------------------
    # Step 3 : final run + comparison
    # ------------------------------------------------------------------
    print("\n[Step 3] Final diagnostic run...")
    _, best_perf = run_one(best_Pinterm, split_ratio, P_low, P_high,
                           T_hot, T_cold, verbose=True)

    print_results(best_perf, best_Pinterm, split_ratio)

    eta_base = 0.3766 if T_hot <= 600+273.15 else 0.3244  # recompression simple reference
    print(f"\n  Δη vs recompression simple @ 900°C : "
          f"+{(best_perf['eta'] - eta_base)*100:.2f} pts")
    print("="*60)


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":

    P_low  = 89.6e5
    #P_high = 215.6e5
    PR = 3
    P_high = P_low * PR
    T_hot  = 900 + 273.15
    T_cold = 37.1 + 273.15
    
    SR = 1                           #47.68/68.11

    optimise(P_low, P_high, T_hot, T_cold,
             split_ratio=SR)