# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:07:49 2026

@author: gregoire.hendrix
"""

from CoolProp.CoolProp import PropsSI
import numpy as np

fluid = 'CO2'

def h(T, P):     return PropsSI('H','T',T,'P',P,fluid)
def s(T, P):     return PropsSI('S','T',T,'P',P,fluid)
def T_hp(H,P):   return PropsSI('T','H',H,'P',P,fluid)
def h_sp(S,P):   return PropsSI('H','P',P,'S',S,fluid)

# =============================================================================
# INPUTS
# =============================================================================

# --- Hot source ---
T_hot_C = 600.0          # °C  ← change ici (600 = sels, >600 = air)
T_hot   = T_hot_C + 273.15

if T_hot_C <= 600:
    pinch_heater = 1.2
    source_label = f"Solar Salt at {T_hot_C:.1f}°C"
else:
    pinch_heater = 5.0
    source_label = f"Air at {T_hot_C:.1f}°C"

# --- Cycle ---
T_cold  = 37.1  + 273.15   # K
P_low   = 89.6e5            # Pa  (turbine exit, P2)
P_high  = 215.6e5           # Pa  (turbine inlet, P1)

# --- Efficiencies (from Hanwha PFD back-calculation) ---
eta_tb  = 0.9102
eta_cp  = 0.8824
eta_rc  = 0.7544

# --- Split ratio ---
split   = 47.68 / 68.11    # fraction → Main Compressor
m_tot   = 1.0              # kg/s reference (scale freely)
m_main  = m_tot * split
m_rc_   = m_tot * (1 - split)

# --- Pressure drops (from Hanwha PFD, Pa) ---
dP_recupHT_hot  = 1.3e5
dP_recupLT_hot  = 0.8e5
dP_cooler_hot   = 1.0e5
dP_recupLT_cold = 1.9e5
dP_recupHT_cold = 2.0e5
dP_heater_cold  = 4.7e5

# --- Recuperator effectivenesses (from Hanwha PFD back-calculation) ---
eta_recupHT = 0.9385
eta_recupLT = 0.8987

# =============================================================================
# STATE POINTS
# =============================================================================

# --- P1 : Turbine inlet ---
T1 = T_hot - pinch_heater
P1 = P_high
h1 = h(T1, P1);  s1 = s(T1, P1)

# --- P2 : Turbine exit ---
P2   = P_low
h2is = h_sp(s1, P2)
h2   = h1 - eta_tb * (h1 - h2is)
T2   = T_hp(h2, P2)

# --- P4 : RecupHT hot exit ---
# Hot side : m_tot, P2→P4 (dP_recupHT_hot)
# Cold side: m_tot, P11→P14 (dP_recupHT_cold)
# Q_max determined by the limiting stream
P4   = P2 - dP_recupHT_hot
P14  = P_high + dP_heater_cold + dP_recupHT_cold + dP_recupLT_cold  # comp must supply this

# Cold side inlet guess: after RecupLT (we'll iterate, but approximate first)
# Use effectiveness definition on hot side as limiting stream
h2_cold_limit = h(T2, P14)           # if cold reaches T_hot_side_in
h_P4_cold_in_guess = h(T_cold + 50 + 273.15, P14)  # rough guess for RecupHT cold inlet

# Q_max hot  = m_tot * (h2 - h(T_cold_in_approx, P2))  — limited by cold reaching T2
# We solve iteratively: RecupHT and RecupLT are coupled through the mixer
# Strategy: assume RecupLT cold outlet (=RecupHT cold inlet) temperature,
#           solve RecupHT, then RecupLT, then check mixer energy balance.

# --- P7 : Compressor inlet (after gas cooler) ---
T7  = T_cold
P7  = P2 - dP_recupHT_hot - dP_recupLT_hot - dP_cooler_hot
h7  = h(T7, P7);  s7 = s(T7, P7)

# --- P9 : Main compressor exit ---
P9   = P_high + dP_heater_cold + dP_recupHT_cold + dP_recupLT_cold
h9is = h_sp(s7, P9)
h9   = h7 + (h9is - h7) / eta_cp
T9   = T_hp(h9, P9)

# --- P12 : Recompressor inlet = RecupLT hot exit ---
P12  = P2 - dP_recupHT_hot - dP_recupLT_hot
# T12 determined by RecupLT — solved below

# --- P13 : Recompressor exit ---
P13  = P9   # same HP rail as main comp exit (mixer)

# =============================================================================
# ITERATIVE SOLVE : RecupHT + RecupLT coupled through mixer
# =============================================================================
# Unknowns: T4 (RecupHT hot exit) = T_recupLT_hot_in
#           T12 (RecupLT hot exit) = T_recomp_su
#           T_mixer_ex             = RecupLT cold inlet
#
# Equations:
#   (1) RecupHT : Q_HT = eta_HT * Q_max_HT
#   (2) RecupLT : Q_LT = eta_LT * Q_max_LT
#   (3) Mixer   : m_main*h9 + m_rc*h13 = m_tot*h_mixer
#   (4) RecupLT cold inlet = mixer exit

# Iterative loop on T_mixer (= RecupLT cold su)
T_mixer = T9 + 10  # initial guess [K] — mixer is warmer than comp exit

for iteration in range(50):

    P_mixer = P9   # mixer exit pressure ≈ P9 (small difference neglected)
    h_mixer = m_main/m_tot * h9  # recomp contribution added after T12 known

    # --- RecupHT ---
    P_recupHT_cold_in  = P9 - dP_recupLT_cold   # after RecupLT cold
    P_recupHT_cold_out = P9 - dP_recupLT_cold - dP_recupHT_cold

    h_recupHT_cold_in  = h(T_mixer, P_mixer)     # = mixer exit ≈ RecupLT cold ex

    # Q_max_HT : min of (hot stream cooling to T_mixer, cold stream heating to T2)
    h_HT_hot_cold_limit  = h(T_mixer, P2)        # hot cools to cold inlet T
    h_HT_cold_hot_limit  = h(T2,     P_recupHT_cold_in)  # cold heats to T2
    Q_max_HT = min(m_tot * (h2 - h_HT_hot_cold_limit),
                   m_tot * (h_HT_cold_hot_limit - h_recupHT_cold_in))

    Q_HT     = eta_recupHT * Q_max_HT
    h4       = h2  - Q_HT / m_tot
    T4       = T_hp(h4, P4)
    h_recupHT_cold_out = h_recupHT_cold_in + Q_HT / m_tot
    T14      = T_hp(h_recupHT_cold_out, P_recupHT_cold_out)

    # --- RecupLT ---
    P5   = P4 - dP_recupLT_hot
    P11  = P9 - dP_recupLT_cold

    h_recupLT_cold_in = h(T_mixer, P_mixer)  # mixer exit = RecupLT cold su

    # Q_max_LT
    h_LT_hot_cold_limit  = h(T_mixer, P4)
    h_LT_cold_hot_limit  = h(T4,     P11)
    Q_max_LT = min(m_tot * (h4  - h_LT_hot_cold_limit),
                   m_tot * (h_LT_cold_hot_limit - h_recupLT_cold_in))

    Q_LT     = eta_recupLT * Q_max_LT
    h12      = h4  - Q_LT / m_tot
    T12      = T_hp(h12, P12)
    h_recupLT_cold_out = h_recupLT_cold_in + Q_LT / m_tot
    T11      = T_hp(h_recupLT_cold_out, P11)

    # --- Recompressor ---
    s12_  = s(T12, P12)
    h13is = h_sp(s12_, P13)
    h13   = h12 + (h13is - h12) / eta_rc
    T13   = T_hp(h13, P13)

    # --- Mixer energy balance ---
    h_mixer_new = (m_main * h9 + m_rc_ * h13) / m_tot
    T_mixer_new = T_hp(h_mixer_new, P_mixer)

    if abs(T_mixer_new - T_mixer) < 0.01:
        T_mixer = T_mixer_new
        break
    T_mixer = 0.5 * T_mixer + 0.5 * T_mixer_new  # damped update

# =============================================================================
# GAS HEATER
# =============================================================================
P_heater_cold_in  = P9 - dP_recupLT_cold - dP_recupHT_cold
P_heater_cold_out = P_high
h_heater_cold_in  = h_recupHT_cold_out
Q_heater          = m_tot * (h1 - h_heater_cold_in)

# =============================================================================
# PERFORMANCE
# =============================================================================
W_turb   = m_tot  * (h1  - h2)
W_comp   = m_main * (h9  - h7)
W_recomp = m_rc_  * (h13 - h12)
W_net    = W_turb - W_comp - W_recomp
eta      = W_net / Q_heater

# =============================================================================
# PRINT
# =============================================================================
print(f"Hot source : {source_label}")
print(f"Iterations : {iteration+1}")
print()
print("="*62)
print(f"{'Point':<25} {'T_model':>9} {'T_PFD':>7} {'P_model':>9} {'P_PFD':>7}")
print("-"*62)

PFD = {
    'P1':  (598.8, 215.6), 'P2':  (488.4,  89.6),
    'P4':  (202.4,  88.3), 'P5':  ( 86.7,  87.5),
    'P7':  ( 37.1,  85.5), 'P9':  ( 76.6, 224.6),
    'P12': ( 86.2,  86.7), 'P13': (184.2, 223.0),
    'P11': (183.7, 222.7), 'P14': (440.9, 220.7),
}

def row(label, T_K, P_Pa, pfd_key):
    Tp, Pp = PFD[pfd_key]
    print(f"{label:<25} {T_K-273.15:>9.1f} {Tp:>7.1f} {P_Pa/1e5:>9.1f} {Pp:>7.1f}")

row("Exp su       (P1)",  T1,  P1,  'P1')
row("Exp ex       (P2)",  T2,  P2,  'P2')
row("RecupHT ex_H (P4)",  T4,  P4,  'P4')
row("RecupLT ex_H (P5)",  T12, P12, 'P5')   # P5 = recomp inlet in PFD
row("Comp su      (P7)",  T7,  P7,  'P7')
row("Comp ex      (P9)",  T9,  P9,  'P9')
row("Recomp su    (P12)", T12, P12, 'P12')
row("Recomp ex    (P13)", T13, P13, 'P13')
row("RecupLT ex_C (P11)", T_hp(h_recupLT_cold_out, P11), P11, 'P11')
row("RecupHT ex_C (P14)", T14, P_recupHT_cold_out, 'P14')
print("="*62)
print()
print(f"  W_turb   = {W_turb/1000:8.2f} kW")
print(f"  W_comp   = {W_comp/1000:8.2f} kW")
print(f"  W_recomp = {W_recomp/1000:8.2f} kW")
print(f"  W_net    = {W_net/1000:8.2f} kW")
print(f"  Q_heater = {Q_heater/1000:8.2f} kW")
print(f"  {'─'*30}")
print(f"  η cycle  = {eta*100:8.2f} %")
print(f"  {'═'*30}")
print(f"\n  Hanwha target (net elec) : 37.2 %")
print(f"  Δ (brut cycle vs net)    : Mech 6% + Gen 97% → brut ≈ 40.8%")