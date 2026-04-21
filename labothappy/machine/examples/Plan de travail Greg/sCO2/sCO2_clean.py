# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:50:41 2026

@author: gregoire.hendrix
"""

# -*- coding: utf-8 -*-
"""
sCO2 Recompression Cycle - simplified version (fixed split = Hanwha)
Calibrated on Hanwha HSC-90-RC-600 (5 MWe, 600°C, 37.2% net electrical)

Designed to be imported and called from a combined-cycle script.
"""

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

CP = 'CO2'

# ─────────────────────────────────────────────
# CoolProp helpers
# ─────────────────────────────────────────────
def hPT(T_C, P_bar):  return PropsSI('H', 'T', T_C + 273.15, 'P', P_bar * 1e5, CP)
def Tph(P_bar, h_J):  return PropsSI('T', 'P', P_bar * 1e5,  'H', h_J,          CP) - 273.15
def hPS(P_bar, s_J):  return PropsSI('H', 'P', P_bar * 1e5,  'S', s_J,          CP)
def sPT(T_C, P_bar):  return PropsSI('S', 'T', T_C + 273.15, 'P', P_bar * 1e5,  CP)

# ─────────────────────────────────────────────
# Cycle constants (Hanwha PFD)
# ─────────────────────────────────────────────
# Pressures [bar]
P_TB_SU  = 215.6   # turbine inlet
P_TB_EX  =  89.6   # turbine exit
P_MC_SU  =  85.7   # main compressor suction
P_MC_EX  = 224.6   # main compressor discharge
P_RC_SU  =  86.7   # recompressor suction
P_RC_EX  = 223.0   # recompressor discharge
P_RHT_HI = 220.0   # RecupHT cold side
P_RHT_LO =  89.0   # RecupHT hot side
P_RLT_HI = 222.0   # RecupLT cold side
P_RLT_LO =  88.0   # RecupLT hot side

# Component efficiencies (back-calculated from Hanwha PFD)
ETA_TB  = 0.9102
ETA_MC  = 0.8824
ETA_RC  = 0.7544
ETA_RHT = 0.9385   # RecupHT effectiveness
ETA_RLT = 0.8987   # RecupLT effectiveness

# Losses and auxiliaries
ETA_GEN  = 0.97
ETA_MECH = 0.94
W_PAR_kW = 115 + 55.1         # air cooler + BOP parasitics [kW]

# Operating point
T_COMP_SU = 37.1              # °C  MC suction
PINCH_HTR =  1.2              # K   main heater pinch
SPLIT     = 47.68 / 68.11     # = 0.6999  (Hanwha reference)


# ─────────────────────────────────────────────
# Main cycle function
# ─────────────────────────────────────────────
def sco2_cycle(T_hot=600.0, m_dot=1.0, split=SPLIT, verbose=False):
    """
    Solve Hanwha-type sCO2 recompression cycle.

    Parameters
    ----------
    T_hot   : HTF (salt/oil/gas) inlet temperature [°C]
    m_dot   : CO2 mass flow rate [kg/s]  — scales absolute powers
    split   : MC branch fraction (default = Hanwha 0.70)
    verbose : print summary line

    Returns
    -------
    dict with:
        eta          : thermal efficiency [-]
        eta_net_elec : net electrical efficiency (incl. gen/mech/parasitics) [-]
        W_net        : net thermodynamic power [kW]  (= m_dot * W_net_specific)
        W_net_elec   : net electrical power [kW]
        Q_heater     : heat input to main heater [kW]
        W_tb, W_mc, W_rc : individual component powers [kW]
        T1..T15      : key cycle temperatures [°C]
    """
    alpha = split
    beta  = 1.0 - split
    T1    = T_hot - PINCH_HTR

    # Turbine
    h1  = hPT(T1, P_TB_SU)
    h2s = hPS(P_TB_EX, sPT(T1, P_TB_SU))
    h2  = h1 - ETA_TB * (h1 - h2s)
    T2  = Tph(P_TB_EX, h2)

    # Main compressor
    h7  = hPT(T_COMP_SU, P_MC_SU)
    h9s = hPS(P_MC_EX, sPT(T_COMP_SU, P_MC_SU))
    h9  = h7 + (h9s - h7) / ETA_MC
    T9  = Tph(P_MC_EX, h9)

    # Coupled residuals: h4 (RecupHT hot exit), h_mix (mixer exit), h_rc (RC exit)
    def residuals(x):
        h4, h_mix, h_rc = x
        T4   = Tph(P_RHT_LO, h4)
        T_mx = Tph(P_RLT_HI, h_mix)

        # RecupLT (asymmetric: m_tot hot, alpha cold)
        Qmax_lt = min(1.0   * (h4 - hPT(T_mx, P_RLT_LO)),
                      alpha * (hPT(T4, P_RLT_HI) - h_mix))
        Q_lt    = ETA_RLT * Qmax_lt
        h5      = h4    - Q_lt
        h11     = h_mix + Q_lt / alpha
        T5      = Tph(P_RLT_LO, h5)
        T11     = Tph(P_RLT_HI, h11)

        # RecupHT (symmetric: m_tot both sides)
        Qmax_ht = min(h2 - hPT(T11, P_RHT_LO),
                      hPT(T2, P_RHT_HI) - h11)
        Q_ht    = ETA_RHT * Qmax_ht
        h4_c    = h2 - Q_ht

        # Recompressor
        h_rc_c  = h5 + (hPS(P_RC_EX, sPT(T5, P_RC_SU)) - h5) / ETA_RC

        # Mixer energy balance
        h_mix_c = alpha * h9 + beta * h_rc

        return [h4 - h4_c, h_mix - h_mix_c, h_rc - h_rc_c]

    x0 = [hPT(202.4, P_RHT_LO), hPT(103.1, P_RLT_HI), hPT(184.2, P_RC_EX)]
    sol, _, ier, _ = fsolve(residuals, x0, full_output=True)
    if ier != 1:
        raise RuntimeError(f"sCO2 cycle solver failed at T_hot={T_hot}°C")

    h4, h_mix, h_rc_ex = sol

    # Reconstruct all states
    T4   = Tph(P_RHT_LO, h4)
    T_mx = Tph(P_RLT_HI, h_mix)
    T_rc = Tph(P_RC_EX,  h_rc_ex)

    Qmax_lt = min(h4 - hPT(T_mx, P_RLT_LO),
                  alpha * (hPT(T4, P_RLT_HI) - h_mix))
    Q_lt    = ETA_RLT * Qmax_lt
    h5      = h4    - Q_lt
    h11     = h_mix + Q_lt / alpha
    T5      = Tph(P_RLT_LO, h5)
    T11     = Tph(P_RLT_HI, h11)

    Qmax_ht = min(h2 - hPT(T11, P_RHT_LO),
                  hPT(T2, P_RHT_HI) - h11)
    Q_ht    = ETA_RHT * Qmax_ht
    h15     = h11 + Q_ht
    T15     = Tph(P_RHT_HI, h15)

    # Specific quantities [J/kg]
    Q_heater_spec = hPT(T1, P_TB_SU) - h15
    W_tb_spec     = h1 - h2
    W_mc_spec     = alpha * (h9 - h7)
    W_rc_spec     = beta  * (h_rc_ex - h5)
    W_net_spec    = W_tb_spec - W_mc_spec - W_rc_spec

    eta = W_net_spec / Q_heater_spec

    # Scale by m_dot to get kW
    W_net    = W_net_spec    * m_dot / 1000.0
    Q_heater = Q_heater_spec * m_dot / 1000.0
    W_tb     = W_tb_spec     * m_dot / 1000.0
    W_mc     = W_mc_spec     * m_dot / 1000.0
    W_rc     = W_rc_spec     * m_dot / 1000.0

    # Net electrical (incl. generator + mechanical losses + parasitics)
    # Parasitics scale proportionally to Q_heater (same assumption as PFD)
    W_par        = W_PAR_kW * (Q_heater / 13450.0)   # 13450 kW = Hanwha Q_heater ref
    W_net_elec   = W_net * ETA_GEN * ETA_MECH - W_par
    eta_net_elec = W_net_elec / Q_heater

    if verbose:
        print(f"sCO2 @ T_hot={T_hot:.1f}°C, m_dot={m_dot:.3f} kg/s [Hanwha ~68 kg/s], split={split:.3f}")
        print(f"  η_thermo   = {eta*100:.2f}%")
        print(f"  η_net_elec = {eta_net_elec*100:.2f}% [Hanwha 37.2%]")
        print(f"  W_net      = {W_net:.1f} kW")
        print(f"  W_net_elec = {W_net_elec:.1f} kW")
        print(f"  Q_heater   = {Q_heater:.1f} kW")

    return {
        'eta':          eta,
        'eta_net_elec': eta_net_elec,
        'W_net':        W_net,
        'W_net_elec':   W_net_elec,
        'Q_heater':     Q_heater,
        'W_tb':         W_tb,
        'W_mc':         W_mc,
        'W_rc':         W_rc,
        'T1':  T1,  'T2':  T2,  'T4':  T4,  'T5':  T5,
        'T9':  T9,  'T11': T11, 'T15': T15,
        'T_mix': T_mx, 'T_rc': T_rc,
    }


# ─────────────────────────────────────────────
# Standalone test (validation vs Hanwha PFD)
# ─────────────────────────────────────────────
if __name__ == "__main__":
    # Validation: find m_dot that gives 5 MWe net electrical (Hanwha design point)
    r = sco2_cycle(T_hot=600, m_dot=71.896, verbose=True)
    print(f"\nExpected (Hanwha): η_net_elec = 37.2%, W_net_elec = 5000 kW")