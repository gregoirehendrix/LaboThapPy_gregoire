# -*- coding: utf-8 -*-
"""
Drainage time calculator for molten salt piping (Solar Salt).

Model: Darcy-Weisbach + Colebrook-White for partially filled open-channel flow.
Fluid properties (rho, mu) fully taken into account via Reynolds number.

Reference: DN70, slope=0.15%, L=300m -> t~30min (field data).

Physics:
    Force balance (open-channel, uniform flow):
        rho*g*S = f_D * rho*v^2 / (8*R_h)     [Pa/m]

        Note: uses 8*R_h instead of 2*D because for open-channel flow the
        hydraulic diameter is D_h = 4*R_h, and Darcy: dP/L = f_D*rho*v^2/(2*D_h)
        => f_D*rho*v^2/(8*R_h)

    Velocity:
        v = sqrt(8*g*R_h*S / f_D)              [m/s]

    Colebrook-White (turbulent, Re > 2300):
        1/sqrt(f_D) = -2*log10(eps/(3.7*D_h) + 2.51/(Re_h*sqrt(f_D)))

        where D_h = 4*R_h  (hydraulic diameter)
              Re_h = rho*v*D_h/mu  (Reynolds based on hydraulic diameter)

    Laminar (Re_h < 2300):
        f_D = 64 / Re_h

    Coupling: v depends on f_D, Re_h depends on v -> fixed-point iteration.

    Flow rate:  Q = v * A_wet                  [m^3/s]
    Pipe vol.:  V = pi*D^2/4 * L               [m^3]
    Time:       t = V / Q                      [s]

Calibrated parameters (adjust in __main__):
    FILL  = 0.60    fill ratio h/D [-]
    EPS   = 46e-6   pipe wall roughness [m]  (commercial steel/stainless)
"""

import math
from labothappy.connector.solar_salt_connector import SolarSaltConnector

G = 9.81   # gravitational acceleration [m/s^2]


# =============================================================================
# Geometry of a partially filled circular pipe
# =============================================================================

def circular_pipe_geometry(D, fill_ratio):
    """
    Cross-section geometry for a circular pipe at fill ratio h/D.

    Parameters
    ----------
    D          : internal diameter [m]
    fill_ratio : h/D  (0 < fill_ratio < 1)

    Returns
    -------
    A   : wetted flow area [m^2]
    P_w : wetted perimeter [m]
    R_h : hydraulic radius A/P_w [m]
    T   : free surface width [m]
    D_h : hydraulic diameter 4*R_h [m]
    """
    h     = fill_ratio * D
    r     = D / 2.0
    theta = math.acos((r - h) / r)   # half-angle subtended by free surface [rad]
    A     = r**2 * (theta - math.sin(theta) * math.cos(theta))
    P_w   = 2.0 * r * theta
    R_h   = A / P_w
    T     = 2.0 * r * math.sin(theta)
    D_h   = 4.0 * R_h
    return A, P_w, R_h, T, D_h


# =============================================================================
# Colebrook-White solver
# =============================================================================

def colebrook_white(Re_h, D_h, eps, max_iter=200, tol=1e-10):
    """
    Darcy friction factor via Colebrook-White (implicit, successive substitution).

    Uses Swamee-Jain as initial guess.

    Parameters
    ----------
    Re_h : Reynolds number based on hydraulic diameter [-]
    D_h  : hydraulic diameter [m]
    eps  : pipe wall roughness [m]

    Returns
    -------
    f_D : Darcy friction factor [-]
    """
    rel_rough = eps / D_h
    # Swamee-Jain initial guess
    f = 0.25 / (math.log10(rel_rough / 3.7 + 5.74 / Re_h**0.9))**2
    for _ in range(max_iter):
        rhs   = -2.0 * math.log10(rel_rough / 3.7 + 2.51 / (Re_h * math.sqrt(f)))
        f_new = 1.0 / rhs**2
        if abs(f_new - f) < tol:
            return f_new
        f = f_new
    return f


# =============================================================================
# Core function
# =============================================================================

def drainage_time(D_int, slope, L, T_C, fill_ratio, eps=46e-6, verbose=True):
    """
    Gravity-driven drainage time of a Solar Salt filled pipe.

    Darcy-Weisbach + Colebrook-White for partially filled open-channel flow.
    Fluid properties rho and mu are fully accounted for via Re_h.

    Parameters
    ----------
    D_int      : internal pipe diameter [m]
    slope      : pipe slope [m/m]
    L          : pipe length [m]
    T_C        : salt temperature [C]
    fill_ratio : flow depth ratio h/D [-]
    eps        : pipe wall roughness [m]  (default 46 um, commercial steel)
    verbose    : print detailed results

    Returns
    -------
    dict: rho, mu, A, R_h, D_h, v, f_D, Re_h, Fr, Q_m3s, V_m3, t_s, t_min
    """
    T_K = T_C + 273.15

    # --- Fluid properties ---
    rho = SolarSaltConnector._rho(T_K)   # [kg/m^3]
    mu  = SolarSaltConnector._mu(T_K)    # [Pa.s]

    # --- Section geometry ---
    A, P_w, R_h, T_surf, D_h = circular_pipe_geometry(D_int, fill_ratio)
    h = fill_ratio * D_int

    # --- Coupled iteration: v <-> f_D via Re_h ---
    # Initial guess: f_D = 0.02 (typical turbulent)
    f_D = 0.02
    for _ in range(300):
        v     = math.sqrt(8.0 * G * R_h * slope / f_D)
        Re_h  = rho * v * D_h / mu
        if Re_h < 2300:
            f_D_new = 64.0 / Re_h
        else:
            f_D_new = colebrook_white(Re_h, D_h, eps)
        if abs(f_D_new - f_D) < 1e-10:
            f_D = f_D_new
            break
        f_D = f_D_new

    # Final consistent values
    v    = math.sqrt(8.0 * G * R_h * slope / f_D)
    Re_h = rho * v * D_h / mu
    Fr   = v / math.sqrt(G * h)

    if Re_h < 2300:
        regime = "laminar"
    elif Re_h < 4000:
        regime = "transition"
    else:
        regime = "turbulent"

    # --- Flow rate and drainage time ---
    Q     = v * A
    V     = (math.pi * D_int**2 / 4.0) * L
    t_s   = V / Q
    t_min = t_s / 60.0

    if verbose:
        print("=" * 61)
        print(f"  DRAINAGE TIME -- Solar Salt @ {T_C:.0f}C")
        print("=" * 61)
        print(f"  D_int      = {D_int*1000:.1f} mm")
        print(f"  Slope      = {slope*100:.3f} %")
        print(f"  Length     = {L:.0f} m")
        print(f"  Fill ratio = {fill_ratio:.2f}  (h = {h*1000:.1f} mm)")
        print(f"  Roughness  = {eps*1e6:.1f} um")
        print("-" * 61)
        print(f"  rho        = {rho:.1f} kg/m^3")
        print(f"  mu         = {mu*1e3:.3f} mPa.s")
        print("-" * 61)
        print(f"  R_h        = {R_h*1000:.2f} mm")
        print(f"  D_h        = {D_h*1000:.2f} mm  (hydraulic diameter)")
        print(f"  A_wet      = {A*1e4:.3f} cm^2")
        print(f"  f_D        = {f_D:.5f}  ({regime})")
        print(f"  Re_h       = {Re_h:.0f}")
        print(f"  Froude     = {Fr:.3f}  ({'subcritical' if Fr < 1 else 'supercritical'})")
        print(f"  Velocity   = {v:.4f} m/s")
        print(f"  Flow rate  = {Q*1e3:.3f} L/s")
        print(f"  Pipe vol.  = {V*1000:.2f} L")
        print(f"  Drainage   = {t_s:.0f} s  -->  {t_min:.1f} min")
        print("=" * 61)

    return {
        "rho": rho, "mu": mu, "A": A, "R_h": R_h, "D_h": D_h,
        "v": v, "f_D": f_D, "Re_h": Re_h, "Fr": Fr,
        "Q_m3s": Q, "V_m3": V,
        "t_s": t_s, "t_min": t_min,
    }


# =============================================================================
# Validation against field reference
# =============================================================================

def validate_reference(fill_ratio, eps):
    """DN70 | slope=0.15% | L=300m | t_ref=30min | T=430C (avg during drain)."""
    print("\n>>> FIELD REFERENCE VALIDATION")
    print(f"    DN70 | slope=0.15% | L=300m | t_ref=30 min")
    print(f"    fill_ratio={fill_ratio:.2f} | eps={eps*1e6:.0f} um")
    res   = drainage_time(D_int=0.070, slope=0.0015, L=300.0, T_C=430.0,
                          fill_ratio=fill_ratio, eps=eps)
    t_ref = 30.0
    err   = (res["t_min"] - t_ref) / t_ref * 100
    print(f"    Model: {res['t_min']:.1f} min  |  Field: {t_ref:.0f} min  |  Error: {err:+.1f}%\n")
    return res


# =============================================================================
# Parametric study
# =============================================================================

def parametric_study(fill_ratio, eps):
    """Sweep over diameters and temperatures for L=200m, slope=0.15%."""
    L     = 200.0
    slope = 0.0015

    diameters_mm = [50, 70, 100, 150, 200]
    temps_C      = [300, 400, 565]

    print(f"\n>>> PARAMETRIC STUDY -- L=200m, slope=0.15%")
    print(f"    fill_ratio={fill_ratio:.2f} | eps={eps*1e6:.0f} um")
    print(f"{'D_int':^12}", end="")
    for T in temps_C:
        print(f"{'T='+str(T)+'C':^22}", end="")
    print()
    print("-" * 78)

    for D_mm in diameters_mm:
        D_int = D_mm / 1000.0
        print(f"  DN{D_mm:<8}", end="")
        for T in temps_C:
            res = drainage_time(D_int=D_int, slope=slope, L=L, T_C=T,
                                fill_ratio=fill_ratio, eps=eps, verbose=False)
            print(f"  {res['t_min']:>7.1f} min ({res['Re_h']:.0f})  ", end="")
        print()
    print("  (Re_h shown in parentheses)")
    print()


# =============================================================================
# Main  --  only place to change parameters
# =============================================================================

if __name__ == "__main__":

    FILL = 0.55    # fill ratio h/D  -- calibrated on field data
    EPS  = 46e-6   # pipe roughness [m] -- commercial steel/stainless

    validate_reference(FILL, EPS)

    print("\n>>> CASE 1: L=200m | slope=0.15% | T=565C")
    drainage_time(D_int=0.070, slope=0.0015, L=200.0, T_C=565.0,
                  fill_ratio=FILL, eps=EPS)

    print("\n>>> CASE 2: L=200m | slope=0.15% | T=300C")
    drainage_time(D_int=0.070, slope=0.0015, L=200.0, T_C=300.0,
                  fill_ratio=FILL, eps=EPS)

    parametric_study(FILL, EPS)