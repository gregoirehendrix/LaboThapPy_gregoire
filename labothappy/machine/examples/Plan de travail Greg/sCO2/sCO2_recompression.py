# -*- coding: utf-8 -*-
"""
sCO2 Recompression cycle — with pressure drops from Hanwha PFD
Derived from validated base cycle.

Pressure drops added (from Hanwha Normal-Design PFD):
  RecupHT  hot side  : 1.3 bar  (P2→P4)
  RecupLT  hot side  : 0.8 bar  (P4→P5)
  GasCooler hot side : 1.0 bar  (P6→P7)
  RecupLT  cold side : 1.9 bar  (P9→P11)
  RecupHT  cold side : 2.0 bar  (P11→P14)
  GasHeater cold side: 4.7 bar  (P14→P1, includes piping)

NOTE: RecursiveCircuit propagates pressures automatically via connectors.
      DP_h / DP_c are passed as parameters to HexCstEffDisc and applied
      internally during discretization.
      The compressor P_high guess must account for the total HP-side drop
      so the turbine inlet pressure remains at P_high.
"""

import numpy as np

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
# PRESSURE DROP VALUES (from Hanwha PFD, bar → Pa)
# =============================================================================
DP = {
    'recupHT_hot':    1.3e5,   # BP side of RecupHT
    'recupLT_hot':    0.8e5,   # BP side of RecupLT
    'cooler_hot':     1.0e5,   # BP side of GasCooler
    'recupLT_cold':   1.9e5,   # HP side of RecupLT
    'recupHT_cold':   2.0e5,   # HP side of RecupHT
    'heater_cold':    4.7e5,   # HP side of GasHeater (incl. piping P14→P1)
}

# Total HP-side pressure drop — compressor must overcome this
DP_HP_total = DP['recupLT_cold'] + DP['recupHT_cold'] + DP['heater_cold']  # 8.6 bar
# Total BP-side pressure drop
DP_BP_total = DP['recupHT_hot'] + DP['recupLT_hot'] + DP['cooler_hot']     # 3.1 bar


# =============================================================================
# CYCLE BUILDER
# =============================================================================
def sCO2_recompression_dp(
        eta_cp, eta_rc, eta_tb, eta_heater, eta_cooler,
        eta_recup_lt, eta_recup_ht, split_ratio,
        HSource, CSource, P_low, P_high, T_c_su, m_dot_CO2,
        pinch_heater=1.2, pinch_recup=2, n_disc=20,
        T_exp_su_guess=None, T_exp_ex_guess=None,
        with_dp=True):
    """
    with_dp=True  : include Hanwha pressure drops
    with_dp=False : ideal cycle (no pressure drops) — for comparison
    """

    if T_exp_su_guess is None: T_exp_su_guess = 598.8 + 273.15
    if T_exp_ex_guess is None: T_exp_ex_guess = 488.4 + 273.15

    # Compressor must deliver P_high + HP losses so turbine sees P_high
    P_comp_ex = P_high + (DP_HP_total if with_dp else 0)

    cycle = RecursiveCircuit('CO2')
    rep_spliter = [split_ratio, 1 - split_ratio]

    compressor   = CompressorCstEff()
    recompressor = CompressorCstEff()
    expander     = ExpanderCstEff()
    gas_heater   = HexCstEffDisc()
    gas_cooler   = HexCstEff()
    recup_lt     = HexCstEffDisc()
    recup_ht     = HexCstEffDisc()
    spliter      = TankSpliter(outlet_repartition=rep_spliter)
    mixer        = TankMixer(n_inlets=2)

    compressor.set_parameters(eta_is=eta_cp)
    recompressor.set_parameters(eta_is=eta_rc)
    expander.set_parameters(eta_is=eta_tb)

    gas_cooler.set_parameters(eta=eta_cooler)

    if with_dp:
        gas_heater.set_parameters(
            eta_max=eta_heater, n_disc=n_disc, Pinch_min=pinch_heater,
            DP_c=DP['heater_cold'], DP_h=0)
        recup_lt.set_parameters(
            eta_max=eta_recup_lt, n_disc=n_disc, Pinch_min=pinch_recup,
            DP_h=DP['recupLT_hot'], DP_c=DP['recupLT_cold'])
        recup_ht.set_parameters(
            eta_max=eta_recup_ht, n_disc=n_disc, Pinch_min=pinch_recup,
            DP_h=DP['recupHT_hot'], DP_c=DP['recupHT_cold'])
    else:
        gas_heater.set_parameters(eta_max=eta_heater, n_disc=n_disc, Pinch_min=pinch_heater)
        recup_lt.set_parameters(eta_max=eta_recup_lt, n_disc=n_disc, Pinch_min=pinch_recup)
        recup_ht.set_parameters(eta_max=eta_recup_ht, n_disc=n_disc, Pinch_min=pinch_recup)

    cycle.add_component(compressor,   "Compressor")
    cycle.add_component(recompressor, "Recompressor")
    cycle.add_component(expander,     "Expander")
    cycle.add_component(gas_heater,   "Gas_Heater")
    cycle.add_component(gas_cooler,   "Gas_Cooler")
    cycle.add_component(recup_lt,     "RecupLT")
    cycle.add_component(recup_ht,     "RecupHT")
    cycle.add_component(spliter,      "Spliter")
    cycle.add_component(mixer,        "Mixer")

    # Topology — unchanged
    cycle.link_components("Expander",     "m-ex",   "RecupHT",    "m-su_H")
    cycle.link_components("RecupHT",      "m-ex_H", "RecupLT",    "m-su_H")
    cycle.link_components("RecupLT",      "m-ex_H", "Spliter",    "m-su")
    cycle.link_components("Spliter",      "m-ex_1", "Gas_Cooler", "m-su_H")
    cycle.link_components("Gas_Cooler",   "m-ex_H", "Compressor", "m-su")
    cycle.link_components("Compressor",   "m-ex",   "Mixer",      "m-su_1")
    cycle.link_components("Spliter",      "m-ex_2", "Recompressor","m-su")
    cycle.link_components("Recompressor", "m-ex",   "Mixer",      "m-su_2")
    cycle.link_components("Mixer",        "m-ex",   "RecupLT",    "m-su_C")
    cycle.link_components("RecupLT",      "m-ex_C", "RecupHT",    "m-su_C")
    cycle.link_components("RecupHT",      "m-ex_C", "Gas_Heater", "m-su_C")
    cycle.link_components("Gas_Heater",   "m-ex_C", "Expander",   "m-su")

    cycle.add_source("Hot_source",  HSource, cycle.components["Gas_Heater"], "m-su_H")
    cycle.add_source("Cold_source", CSource, cycle.components["Gas_Cooler"], "m-su_C")

    # -------------------------------------------------------------------------
    # Guesses — account for pressure drops
    # -------------------------------------------------------------------------
    T_recupHT_ex_H = T_c_su + 165
    T_recupHT_ex_C = T_exp_su_guess - 160

    cycle.set_cycle_guess(target='Expander:su',
                          m_dot=m_dot_CO2, T=T_exp_su_guess, p=P_high)
    cycle.set_cycle_guess(target='Expander:ex',
                          T=T_exp_ex_guess, p=P_low)

    cycle.set_cycle_guess(target='Compressor:su',
                          m_dot=m_dot_CO2 * split_ratio, T=T_c_su,
                          p=P_low - (DP_BP_total if with_dp else 0))
    cycle.set_cycle_guess(target='Compressor:ex',
                          T=T_c_su + 60, p=P_comp_ex)

    cycle.set_cycle_guess(target='Recompressor:su',
                          m_dot=m_dot_CO2 * (1-split_ratio), T=T_c_su + 60,
                          p=P_low - (DP['recupHT_hot'] + DP['recupLT_hot'] if with_dp else 0))
    cycle.set_cycle_guess(target='Recompressor:ex',
                          T=T_c_su + 110, p=P_comp_ex)

    cycle.set_cycle_guess(target='RecupHT:su_H',
                          m_dot=m_dot_CO2, T=T_exp_ex_guess, p=P_low)
    cycle.set_cycle_guess(target='RecupHT:ex_H',
                          m_dot=m_dot_CO2, T=T_recupHT_ex_H,
                          p=P_low - (DP['recupHT_hot'] if with_dp else 0))
    cycle.set_cycle_guess(target='RecupHT:su_C',
                          m_dot=m_dot_CO2, T=T_c_su + 100,
                          p=P_comp_ex - (DP['recupLT_cold'] if with_dp else 0))
    cycle.set_cycle_guess(target='RecupHT:ex_C',
                          T=T_recupHT_ex_C,
                          p=P_comp_ex - (DP['recupLT_cold'] + DP['recupHT_cold'] if with_dp else 0))

    cycle.set_cycle_guess(target='RecupLT:su_H',
                          m_dot=m_dot_CO2, T=T_recupHT_ex_H,
                          p=P_low - (DP['recupHT_hot'] if with_dp else 0))
    cycle.set_cycle_guess(target='RecupLT:ex_H',
                          T=T_c_su + 50,
                          p=P_low - (DP['recupHT_hot'] + DP['recupLT_hot'] if with_dp else 0))
    cycle.set_cycle_guess(target='RecupLT:su_C',
                          m_dot=m_dot_CO2, T=T_c_su + 50, p=P_comp_ex)
    cycle.set_cycle_guess(target='RecupLT:ex_C',
                          T=T_c_su + 150,
                          p=P_comp_ex - (DP['recupLT_cold'] if with_dp else 0))

    for tgt in ['Expander:ex', 'Compressor:ex', 'Gas_Heater:ex_C', 'Gas_Cooler:ex_H',
                'RecupHT:su_H', 'RecupHT:ex_H', 'RecupLT:ex_H', 'Spliter:ex_2']:
        cycle.set_residual_variable(target=tgt, variable='h', tolerance=1e-6)
        cycle.set_residual_variable(target=tgt, variable='p', tolerance=1e-6)

    return cycle, compressor, recompressor, expander, gas_heater, gas_cooler, recup_lt, recup_ht


# =============================================================================
# PERFORMANCE
# =============================================================================
def compute_cycle_performance(compressor, expander, gas_heater, gas_cooler,
                               recompressor=None, recup_lt=None, recup_ht=None):
    m_dot_cp  = compressor.su.m_dot
    m_dot_tot = expander.su.m_dot

    W_comp   = m_dot_cp  * (compressor.ex.h - compressor.su.h) / 1000
    W_exp    = m_dot_tot * (expander.su.h   - expander.ex.h)   / 1000
    W_recomp = (recompressor.su.m_dot *
                (recompressor.ex.h - recompressor.su.h) / 1000
                if recompressor is not None else 0.0)

    W_net    = W_exp - W_comp - W_recomp
    Q_heater = m_dot_tot * (gas_heater.ex_C.h - gas_heater.su_C.h) / 1000

    return {
        'W_comp':    W_comp,
        'W_recomp':  W_recomp,
        'W_exp':     W_exp,
        'W_net':     W_net,
        'Q_heater':  Q_heater,
        'Q_cooler':  Q_heater - W_net,
        'Q_recup_LT': m_dot_tot * (recup_lt.ex_C.h - recup_lt.su_C.h) / 1000,
        'Q_recup_HT': m_dot_tot * (recup_ht.ex_C.h - recup_ht.su_C.h) / 1000,
        'eta': W_net / Q_heater,
    }


def print_results(perf, label=""):
    print(f"\n=== Cycle performance {label} ===")
    print(f"  W_comp    : {perf['W_comp']:8.2f} kW")
    print(f"  W_recomp  : {perf['W_recomp']:8.2f} kW")
    print(f"  W_exp     : {perf['W_exp']:8.2f} kW")
    print(f"  W_net     : {perf['W_net']:8.2f} kW")
    print(f"  Q_heater  : {perf['Q_heater']:8.2f} kW")
    print(f"  Q_cooler  : {perf['Q_cooler']:8.2f} kW")
    print(f"  Q_recupLT : {perf['Q_recup_LT']:8.2f} kW")
    print(f"  Q_recupHT : {perf['Q_recup_HT']:8.2f} kW")
    print(f"  {'─'*30}")
    print(f"  eta       : {perf['eta']*100:8.2f} %")
    print(f"  {'═'*30}")


def print_pressures(compressor, recompressor, expander,
                    recup_lt, recup_ht, gas_heater, gas_cooler, label=""):
    print(f"\n=== Pressure levels {label} ===")
    print(f"  Expander su        : {expander.su.p/1e5:7.2f} bar")
    print(f"  Expander ex        : {expander.ex.p/1e5:7.2f} bar")
    print(f"  RecupHT  su_H      : {recup_ht.su_H.p/1e5:7.2f} bar")
    print(f"  RecupHT  ex_H      : {recup_ht.ex_H.p/1e5:7.2f} bar")
    print(f"  RecupLT  su_H      : {recup_lt.su_H.p/1e5:7.2f} bar")
    print(f"  RecupLT  ex_H      : {recup_lt.ex_H.p/1e5:7.2f} bar")
    print(f"  GasCooler su_H     : {gas_cooler.su_H.p/1e5:7.2f} bar")
    print(f"  Compressor su      : {compressor.su.p/1e5:7.2f} bar")
    print(f"  Compressor ex      : {compressor.ex.p/1e5:7.2f} bar")
    print(f"  RecupLT  su_C      : {recup_lt.su_C.p/1e5:7.2f} bar")
    print(f"  RecupLT  ex_C      : {recup_lt.ex_C.p/1e5:7.2f} bar")
    print(f"  RecupHT  su_C      : {recup_ht.su_C.p/1e5:7.2f} bar")
    print(f"  RecupHT  ex_C      : {recup_ht.ex_C.p/1e5:7.2f} bar")
    print(f"  GasHeater su_C     : {gas_heater.su_C.p/1e5:7.2f} bar")
    print(f"  GasHeater ex_C     : {gas_heater.ex_C.p/1e5:7.2f} bar")


# =============================================================================
# MAIN — compare with and without pressure drops
# =============================================================================
if __name__ == "__main__":

    P_low    = 89.6e5
    P_high   = 215.6e5
    T_cold   = 37.1 + 273.15
    T_hot    = 600 + 273.15
    m_dot    = 1.0
    split_ratio = 47.68 / 68.11

    T_salt_limit = 600.0 + 273.15
    if T_hot <= T_salt_limit:
        def make_sources():
            H = SolarSaltConnector()
            H.set_properties(T=T_hot, p=2e5, m_dot=100)
            C = MassConnector()
            C.set_properties(fluid='air', T=T_cold, p=101325, m_dot=500)
            return H, C
        PINCH_H = 1.2
        eta_h   = 0.999
        print(f"Hot source : Solar Salt at {T_hot-273.15:.1f} °C")
    else:
        def make_sources():
            H = MassConnector()
            H.set_properties(fluid='air', T=T_hot, p=1.2e5, m_dot=500)
            C = MassConnector()
            C.set_properties(fluid='air', T=T_cold, p=101325, m_dot=500)
            return H, C
        PINCH_H = 5.0
        eta_h   = 0.999
        print(f"Hot source : Air at {T_hot-273.15:.1f} °C")

    T_su_guess = T_hot - PINCH_H
    T_ex_guess = T_su_guess - 100.0

    common = dict(
        eta_cp=0.8824, eta_rc=0.7544, eta_tb=0.9102,
        eta_heater=eta_h, eta_cooler=0.999,
        eta_recup_lt=0.8987, eta_recup_ht=0.9385,   # ← mis à jour
        split_ratio=split_ratio,
        P_low=P_low, P_high=P_high,
        T_c_su=T_cold, m_dot_CO2=m_dot,
        pinch_heater=PINCH_H, pinch_recup=2.0, n_disc=20,
        T_exp_su_guess=T_su_guess, T_exp_ex_guess=T_ex_guess,
    )

    # ------------------------------------------------------------------
    # Run 1 : no pressure drops (reference)
    # ------------------------------------------------------------------
    H, C = make_sources()
    cycle1, cp1, rc1, exp1, gh1, gc1, rlt1, rht1 = sCO2_recompression_dp(
        HSource=H, CSource=C, with_dp=False, **common)
    cycle1.solve()
    perf1 = compute_cycle_performance(cp1, exp1, gh1, gc1, rc1, rlt1, rht1)
    print_results(perf1, "(no ΔP)")
    print_pressures(cp1, rc1, exp1, rlt1, rht1, gh1, gc1, "(no ΔP)")

    # ------------------------------------------------------------------
    # Run 2 : with Hanwha pressure drops
    # ------------------------------------------------------------------
    H, C = make_sources()
    cycle2, cp2, rc2, exp2, gh2, gc2, rlt2, rht2 = sCO2_recompression_dp(
        HSource=H, CSource=C, with_dp=True, **common)
    cycle2.solve()
    perf2 = compute_cycle_performance(cp2, exp2, gh2, gc2, rc2, rlt2, rht2)
    print_results(perf2, "(with Hanwha ΔP)")
    print_pressures(cp2, rc2, exp2, rlt2, rht2, gh2, gc2, "(with Hanwha ΔP)")
    
    # Colle ça après cycle2.solve()

    print("\n" + "="*60)
    print("STATE POINTS — comparison with Hanwha PFD")
    print("="*60)
    print(f"{'Point':<20} {'T_model':>10} {'T_PFD':>8} {'P_model':>10} {'P_PFD':>8}")
    print("-"*60)
    
    # Turbine
    print(f"{'Exp su (P1)':<20} {exp2.su.T-273.15:>10.1f} {'598.8':>8} {exp2.su.p/1e5:>10.1f} {'215.6':>8}")
    print(f"{'Exp ex (P2)':<20} {exp2.ex.T-273.15:>10.1f} {'488.4':>8} {exp2.ex.p/1e5:>10.1f} {'89.6':>8}")
    
    # RecupHT hot side
    print(f"{'RecupHT su_H (P2)':<20} {rht2.su_H.T-273.15:>10.1f} {'488.4':>8} {rht2.su_H.p/1e5:>10.1f} {'89.6':>8}")
    print(f"{'RecupHT ex_H (P4)':<20} {rht2.ex_H.T-273.15:>10.1f} {'202.4':>8} {rht2.ex_H.p/1e5:>10.1f} {'88.3':>8}")
    
    # RecupLT hot side
    print(f"{'RecupLT su_H (P4)':<20} {rlt2.su_H.T-273.15:>10.1f} {'202.4':>8} {rlt2.su_H.p/1e5:>10.1f} {'88.3':>8}")
    print(f"{'RecupLT ex_H (P5)':<20} {rlt2.ex_H.T-273.15:>10.1f} {'86.7':>8} {rlt2.ex_H.p/1e5:>10.1f} {'87.5':>8}")
    
    # GasCooler + Compressor
    print(f"{'Cooler ex / Comp su (P7)':<20} {cp2.su.T-273.15:>10.1f} {'37.1':>8} {cp2.su.p/1e5:>10.1f} {'85.5':>8}")
    print(f"{'Comp ex (P9)':<20} {cp2.ex.T-273.15:>10.1f} {'76.6':>8} {cp2.ex.p/1e5:>10.1f} {'224.6':>8}")
    
    # Recompressor
    print(f"{'Recomp su (P12)':<20} {rc2.su.T-273.15:>10.1f} {'86.2':>8} {rc2.su.p/1e5:>10.1f} {'86.7':>8}")
    print(f"{'Recomp ex (P13)':<20} {rc2.ex.T-273.15:>10.1f} {'184.2':>8} {rc2.ex.p/1e5:>10.1f} {'223.0':>8}")
    
    # RecupLT cold side
    print(f"{'RecupLT su_C (P10≈P9)':<20} {rlt2.su_C.T-273.15:>10.1f} {'76.5':>8} {rlt2.su_C.p/1e5:>10.1f} {'224.2':>8}")
    print(f"{'RecupLT ex_C (P11)':<20} {rlt2.ex_C.T-273.15:>10.1f} {'183.7':>8} {rlt2.ex_C.p/1e5:>10.1f} {'222.7':>8}")
    
    # RecupHT cold side
    print(f"{'RecupHT su_C (P11)':<20} {rht2.su_C.T-273.15:>10.1f} {'183.7':>8} {rht2.su_C.p/1e5:>10.1f} {'222.7':>8}")
    print(f"{'RecupHT ex_C (P14)':<20} {rht2.ex_C.T-273.15:>10.1f} {'440.9':>8} {rht2.ex_C.p/1e5:>10.1f} {'220.7':>8}")
    
    # GasHeater cold side
    print(f"{'GasHeater su_C (P14)':<20} {gh2.su_C.T-273.15:>10.1f} {'440.9':>8} {gh2.su_C.p/1e5:>10.1f} {'220.7':>8}")
    print(f"{'GasHeater ex_C (P1)':<20} {gh2.ex_C.T-273.15:>10.1f} {'598.8':>8} {gh2.ex_C.p/1e5:>10.1f} {'215.6':>8}")
    print("="*60)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "="*50)
    print("COMPARISON SUMMARY")
    print("="*50)
    print(f"\n  Hanwha PFD target : 37.2 % (net, with mech. losses)")
    print(f"  η with ΔP    : {perf2['eta']*100:.2f} %")
    print(f"  Δη           : {37.2-perf2['eta']*100:+.2f} pts")

    print("="*50)
    
    """
    η_recupHT = 0.9385
    η_recupLT = 0.8987
    η_turbine  = 0.9102
    η_maincomp = 0.8824
    η_recomp   = 0.7544
    """