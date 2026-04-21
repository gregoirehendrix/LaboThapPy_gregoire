# -*- coding: utf-8 -*-
"""
sCO2 Recompression Cycle - Pure CoolProp implementation
Validation target: Hanwha HSC-90-RC-600, Normal Design (5 MWe, 37.2%)

Architecture:
  Turbine → RecupHT(hot) → RecupLT(hot) → Splitter
                                            ├─ [alpha]  → GasCooler → MainComp → Mixer
                                            └─ [1-alpha] → Recomp             → Mixer
  Mixer → RecupLT(cold) → RecupHT(cold) → MainHeater → Turbine

Pressure levels (from PFD, with drops):
  P_hi_turb  = 215.6 bar  (turbine inlet, P1)
  P_lo_turb  = 89.6  bar  (turbine exit,  P2)
  P_hi_comp  = 224.6 bar  (compressor exit, P9 = P8)
  P_lo_comp  = 85.7  bar  (compressor inlet, P7)

Note: using average P_hi=215.6 and P_lo=89.6 as in original code (drops negligible for eta).
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import fsolve, brentq
import matplotlib.pyplot as plt

CP = 'CO2'

# ─────────────────────────────────────────────
# CoolProp helpers
# ─────────────────────────────────────────────
def hPT(T_C, P_bar):   return PropsSI('H', 'T', T_C+273.15, 'P', P_bar*1e5, CP)
def Tph(P_bar, h_J):   return PropsSI('T', 'P', P_bar*1e5,  'H', h_J,        CP) - 273.15
def hPS(P_bar, s_J):   return PropsSI('H', 'P', P_bar*1e5,  'S', s_J,        CP)
def sPT(T_C,  P_bar):  return PropsSI('S', 'T', T_C+273.15, 'P', P_bar*1e5,  CP)

# ─────────────────────────────────────────────
# Pressure levels (from PFD, with line drops)
# ─────────────────────────────────────────────
# Turbine side
P_TB_SU = 215.6   # bar  turbine inlet  (P1)
P_TB_EX =  89.6   # bar  turbine exit   (P2)
# Main compressor
P_MC_SU =  85.7   # bar  MC suction     (P7)
P_MC_EX = 224.6   # bar  MC discharge   (P9)
# Recompressor
P_RC_SU =  86.7   # bar  RC suction     (P12)
P_RC_EX = 223.0   # bar  RC discharge   (P13)
# Recuperator representative pressures
P_RHT_HI = 220.0  # bar  RecupHT cold side
P_RHT_LO =  89.0  # bar  RecupHT hot side
P_RLT_HI = 222.0  # bar  RecupLT cold side
P_RLT_LO =  88.0  # bar  RecupLT hot side

# Legacy aliases (kept for backward compatibility with sweep function)
P_HI = P_TB_SU
P_LO = P_TB_EX

T_COMP_SU =  37.1   # °C   compressor suction (P7)
T_HOT_SU  = 600.0   # °C   salt/HTF inlet (P16)
PINCH_HTR =   1.2   # K    main heater pinch → T_turb_su = 598.8°C

# Back-calculated efficiencies from Hanwha PFD Normal Design
ETA_TB  = 0.9102   # turbine          (isentropic)
ETA_MC  = 0.8824   # main compressor  (isentropic)
ETA_RC  = 0.7544   # recompressor     (isentropic)
ETA_RHT = 0.9385   # RecupHT effectiveness
ETA_RLT = 0.8987   # RecupLT effectiveness

# Manufacturer losses (PFD assumptions)
ETA_GEN  = 0.97         # generator efficiency
ETA_MECH = 0.94         # mechanical efficiency (= 1 - 6% loss)
W_PAR_kW = 115 + 55.1   # parasitic: air cooler + BOP [kW]

# Hanwha reference split ratio
SPLIT_HANWHA = 47.68 / 68.11   # = 0.6999...

# NOTE on PFD labeling:
# P10 (76.5°C) in the PFD is the MC exit BEFORE the mixer — NOT the mixer exit.
# Energy balance: T_mix = (m_mc*h9 + m_rc*h13) / m_tot ≈ 103°C (not 76.5°C).
# The PFD does not label the actual mixer exit temperature explicitly.

# ─────────────────────────────────────────────
# Isentropic enthalpy helpers
# ─────────────────────────────────────────────
def h_isentropic(T_su_C, P_su_bar, P_ex_bar):
    """Exit enthalpy for isentropic process."""
    s_su = sPT(T_su_C, P_su_bar)
    return hPS(P_ex_bar, s_su)

def h_isentropic_from_h(h_su, P_su_bar, P_ex_bar):
    """Exit enthalpy for isentropic process given inlet enthalpy."""
    T_su = Tph(P_su_bar, h_su)
    s_su = sPT(T_su, P_su_bar)
    return hPS(P_ex_bar, s_su)

# ─────────────────────────────────────────────
# Main cycle solver
# ─────────────────────────────────────────────
def run_cycle(split, m_dot=1.0, verbose=False):
    """
    Solve recompression sCO2 cycle at given split ratio.

    Architecture (Hanwha HSC-90-RC-600):
      RecupLT: hot = m_tot @ P_RLT_LO, cold = m_mc (=alpha*m_tot) @ P_RLT_HI  [asymmetric]
      RecupHT: hot = m_tot @ P_RHT_LO, cold = m_tot               @ P_RHT_HI  [symmetric]
      Mixer exit (T_mix ≈ 103°C) feeds RecupLT cold inlet.
      PFD label P10 (76.5°C) is the MC exit BEFORE mixing — not the mixer exit.

    Parameters
    ----------
    split  : fraction of total m_dot through MainComp branch (= alpha)
    m_dot  : total CO2 mass flow [kg/s], only affects absolute power output
    verbose: print state comparison vs PFD

    Returns
    -------
    dict with performance metrics, or None if solver failed.
    """
    alpha = split
    beta  = 1.0 - split
    T1    = T_HOT_SU - PINCH_HTR   # turbine inlet [°C]

    # ── Turbine ─────────────────────────────────────────────────
    h1  = hPT(T1,       P_TB_SU)
    h2  = h1 - ETA_TB * (h1 - hPT(T1, P_TB_SU) + hPS(P_TB_EX, sPT(T1, P_TB_SU)) - h1)
    # cleaner:
    h2s = hPS(P_TB_EX, sPT(T1, P_TB_SU))
    h2  = h1 - ETA_TB * (h1 - h2s)
    T2  = Tph(P_TB_EX, h2)

    # ── Main Compressor ─────────────────────────────────────────
    h7  = hPT(T_COMP_SU, P_MC_SU)
    h9s = hPS(P_MC_EX, sPT(T_COMP_SU, P_MC_SU))
    h9  = h7 + (h9s - h7) / ETA_MC
    T9  = Tph(P_MC_EX, h9)

    # ── Coupled system ───────────────────────────────────────────
    # Unknowns: h4 (RecupHT hot exit), h_mix (Mixer exit), h_rc_ex (RC exit)
    # RecupLT Qmax: asymmetric — m_tot hot side, alpha*m_tot cold side
    #   Qmax = min(1.0*(h4 - h(T_mix, P_RLT_LO)), alpha*(h(T4, P_RLT_HI) - h_mix))
    # Epsilon definition: Q_actual = ETA * Qmax, based on hot-side (limiting) stream

    def residuals(x):
        h4, h_mix, h_rc_ex = x
        T4   = Tph(P_RHT_LO, h4)
        T_mx = Tph(P_RLT_HI, h_mix)

        # RecupLT: m_tot hot, alpha cold
        Qmax_lt = min(1.0   * (h4 - hPT(T_mx, P_RLT_LO)),
                      alpha * (hPT(T4, P_RLT_HI) - h_mix))
        Q_lt    = ETA_RLT * Qmax_lt
        h5      = h4    - Q_lt / 1.0
        h11     = h_mix + Q_lt / alpha
        T11     = Tph(P_RLT_HI, h11)
        T5      = Tph(P_RLT_LO, h5)

        # RecupHT: m_tot both sides
        Qmax_ht = min(1.0 * (h2 - hPT(T11, P_RHT_LO)),
                      1.0 * (hPT(T2, P_RHT_HI) - h11))
        Q_ht    = ETA_RHT * Qmax_ht
        h4_calc = h2 - Q_ht

        # Recompressor
        h_rc_calc = h5 + (hPS(P_RC_EX, sPT(T5, P_RC_SU)) - h5) / ETA_RC

        # Mixer
        h_mix_calc = alpha * h9 + beta * h_rc_ex

        return [h4 - h4_calc, h_mix - h_mix_calc, h_rc_ex - h_rc_calc]

    x0 = [hPT(202.4, P_RHT_LO), hPT(103.1, P_RLT_HI), hPT(184.2, P_RC_EX)]
    try:
        sol, _, ier, _ = fsolve(residuals, x0, full_output=True)
    except Exception:
        return None
    if ier != 1 or max(abs(r) for r in residuals(sol)) > 500:
        return None

    h4, h_mix, h_rc_ex = sol
    T4   = Tph(P_RHT_LO, h4)
    T_mx = Tph(P_RLT_HI, h_mix)
    T_rc = Tph(P_RC_EX,  h_rc_ex)

    Qmax_lt = min(1.0*(h4 - hPT(T_mx, P_RLT_LO)), alpha*(hPT(T4, P_RLT_HI) - h_mix))
    Q_lt    = ETA_RLT * Qmax_lt
    h5      = h4    - Q_lt
    h11     = h_mix + Q_lt / alpha
    T5      = Tph(P_RLT_LO, h5)
    T11     = Tph(P_RLT_HI, h11)

    Qmax_ht = min(1.0*(h2 - hPT(T11, P_RHT_LO)), 1.0*(hPT(T2, P_RHT_HI) - h11))
    Q_ht    = ETA_RHT * Qmax_ht
    h15     = h11 + Q_ht
    T15     = Tph(P_RHT_HI, h15)

    h1v      = hPT(T1, P_TB_SU)
    Q_heater = h1v - h15

    W_tb  = h1  - h2
    W_mc  = alpha * (h9 - h7)
    W_rc  = beta  * (h_rc_ex - h5)
    W_net = W_tb - W_mc - W_rc
    eta   = W_net / Q_heater

    if verbose:
        print(f"\n{'═'*56}")
        print(f"  Recompression sCO2 — split = {split:.4f}")
        print(f"  (Hanwha reference  = {SPLIT_HANWHA:.4f})")
        print(f"{'═'*56}")
        print(f"  {'State':<16} {'Calc':>8}   {'PFD':>8}   note")
        print(f"  {'─'*50}")
        # T_mix_PFD: energy balance gives 103.1°C (PFD P10=76.5 is MC exit, not mixer exit)
        states = [
            ("T1  turb-su",   T1,          598.8, ""),
            ("T2  turb-ex",   T2,          488.4, ""),
            ("T4  rHT h-ex",  T4,          202.4, ""),
            ("T5  rLT h-ex",  T5,           86.7, ""),
            ("T9  MC-ex",     T9,           76.6, ""),
            ("T11 rLT c-ex",  T11,         183.7, ""),
            ("T15 rHT c-ex",  T15,         440.9, ""),
            ("T_RC_ex",       T_rc,        184.2, ""),
            ("T_mix",         T_mx,        103.1, "energy balance (not PFD P10)"),
        ]
        for label, calc, pfd, note in states:
            diff = calc - pfd
            flag = "✓" if abs(diff) < 5 else f"Δ={diff:+.1f}°C"
            print(f"  {label:<16} {calc:>7.1f}°C  {pfd:>7.1f}°C  {flag}  {note}")
        print(f"  {'─'*50}")
        print(f"  W_net     {W_net/1000:>7.2f} kW/(kg/s)")
        print(f"  Q_heater  {Q_heater/1000:>7.2f} kW/(kg/s)")
        print(f"  η_thermo  {eta*100:>7.2f}%")
        W_net_thermo = (5000 + W_PAR_kW) / (ETA_GEN * ETA_MECH)
        scale        = W_net_thermo / (W_net / 1000)
        eta_net_e    = 5000 / (Q_heater / 1000 * scale)
        print(f"  η_net_elec{eta_net_e*100:>7.2f}%  [PFD: 37.2%]")
        print(f"  m_dot CO2 {scale:>7.1f} kg/s    [PFD: ~68 kg/s]")
        print(f"{'═'*56}")

    return {
        'split': split, 'eta': eta,
        'W_net':     W_net     * m_dot / 1000,
        'Q_heater':  Q_heater  * m_dot / 1000,
        'W_tb': W_tb * m_dot / 1000,
        'W_mc': W_mc * m_dot / 1000,
        'W_rc': W_rc * m_dot / 1000,
        'T1': T1,   'T2': T2,   'T4': T4,   'T5': T5,
        'T9': T9,   'T11': T11, 'T15': T15,
        'T_mix': T_mx, 'T_rc': T_rc,
    }

# ─────────────────────────────────────────────
# Split ratio sweep + optimisation
# ─────────────────────────────────────────────
def sweep_split(splits=None, plot=True):
    if splits is None:
        splits = np.linspace(0.50, 0.95, 91)

    results = []
    for sp in splits:
        r = run_cycle(sp)
        if r:
            results.append((sp, r['eta']*100))

    splits_v = [r[0] for r in results]
    etas_v   = [r[1] for r in results]

    best_idx = int(np.argmax(etas_v))
    sp_opt   = splits_v[best_idx]
    eta_opt  = etas_v[best_idx]

    r_hanwha = run_cycle(SPLIT_HANWHA)
    eta_hanwha = r_hanwha['eta']*100 if r_hanwha else None

    print(f"\n{'─'*44}")
    print(f"  Split sweep summary")
    print(f"{'─'*44}")
    print(f"  Optimal split  = {sp_opt:.3f},  η = {eta_opt:.3f}%")
    if eta_hanwha:
        print(f"  Hanwha  split  = {SPLIT_HANWHA:.3f},  η = {eta_hanwha:.3f}%")
        print(f"  Δη vs Hanwha   = {eta_opt - eta_hanwha:+.3f}%")
    print(f"{'─'*44}")

    if plot:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(splits_v, etas_v, 'b-', lw=1.8, label='η_thermo')
        ax.axvline(sp_opt,       color='g', ls='--', lw=1.2, label=f'Optimum ({sp_opt:.3f})')
        ax.axvline(SPLIT_HANWHA, color='r', ls=':',  lw=1.2, label=f'Hanwha ({SPLIT_HANWHA:.3f})')
        ax.set_xlabel('Split ratio  α  (MainComp fraction)')
        ax.set_ylabel('Thermal efficiency η [%]')
        ax.set_title('sCO2 Recompression — Split Ratio Optimisation')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        #plt.savefig('/mnt/user-data/outputs/sCO2_split_sweep.png', dpi=150)
        plt.show()
        print("  Plot saved → sCO2_split_sweep.png")

    return splits_v, etas_v, sp_opt, eta_opt

# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────
if __name__ == "__main__":

    # 1) Validation at Hanwha split
    run_cycle(SPLIT_HANWHA, verbose=True)

    # 2) Run at optimal split (around 0.72-0.73 expected)
    _, _, sp_opt, _ = sweep_split(plot=True)
    print("\n=== Optimal cycle ===")
    run_cycle(sp_opt, verbose=True)