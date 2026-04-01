# -*- coding: utf-8 -*-
"""
Standalone validation script — Molten salt piping losses
Energy-equivalent thermal loss model (per diameter class)

Author: Grégoire Hendrix
Date  : Apr 2026
"""

import numpy as np


# --- Temperatures ---
T_hot = 565 + 273.15      # K  (solar field outlet)
T_amb = 25 + 273.15       # K  (ambient temperature)

# --- Molten salt properties (temperature dependent Cp) ---
T_C = T_hot - 273.15
Cp_salt = 1443.0 + 0.172 * T_C    # J/kg/K

# --- Network geometry (from layout model) ---
# RMS length per diameter [m]
rms_per_diam = {
    '3"'     : 495.0,
    '2 1/2"' : 200.0,
    '2"'     : 200.0,
    '1 1/4"' : 200.0
}

# Total length per diameter [m]
lengths_per_diam = {
    '3"'     : 16000.0,
    '2 1/2"' : 4000.0,
    '2"'     : 4000.0,
    '1 1/4"' : 4000.0
}

# --- Equivalent thermal flow per diameter ---
# IMPORTANT:
# This does NOT represent the flow in a single pipe.
# It represents the TOTAL equivalent number of CSP units
# transported by all pipes of this diameter together.
equivalent_units_per_diam = {
    '3"'     : 100,   # e.g. ~20 pipes × 5 units each
    '2 1/2"' : 25,
    '2"'     : 10,
    '1 1/4"' : 5
}

# --- Flow per CSP unit ---
m_dot_unit = 6.39   # kg/s per CSP unit

# --- Pipe external diameters [m] ---
DIAMETER_M = {
    '3"'     : 0.089,
    '2 1/2"' : 0.075,
    '2"'     : 0.060,
    '1 1/4"' : 0.042
}

# --- Overall heat transfer coefficients [W/m²/K] ---
U_PIPE = {
    '3"'     : 0.9,
    '2 1/2"' : 1.0,
    '2"'     : 1.1,
    '1 1/4"' : 1.2
}


def compute_total_molten_salt_piping_losses(
    T_hot, T_amb,
    rms_per_diam, lengths_per_diam,
    equivalent_units_per_diam, m_dot_unit,
    DIAMETER_M, U_PIPE,
    Cp_salt
):
    """
    Computes total molten salt piping losses using an
    exponential (NTU-based) energy-equivalent model.

    Returns:
        Q_loss_tot [W]
    """

    Q_loss_tot = 0.0

    print("\n=== Molten salt piping losses per diameter ===")

    for d in rms_per_diam:

        # Geometry
        L_rms = rms_per_diam[d]
        L_tot = lengths_per_diam[d]
        n_k   = L_tot / L_rms                   # number of equivalent segments

        # Equivalent mass flow
        m_dot_k = equivalent_units_per_diam[d] * m_dot_unit

        # Thermal parameters
        D_k = DIAMETER_M[d]
        U_k = U_PIPE[d]

        A_rms = np.pi * D_k * L_rms
        NTU   = U_k * A_rms / (m_dot_k * Cp_salt)

        # Thermal losses
        Q_k = (
            n_k
            * m_dot_k
            * Cp_salt
            * (T_hot - T_amb)
            * (1.0 - np.exp(-NTU))
        )

        Q_loss_tot += Q_k

        print(f"\nDiameter {d}")
        print(f"  n_segments (equiv) : {n_k:.2f}")
        print(f"  m_dot_equiv        : {m_dot_k:.2f} kg/s")
        print(f"  L_rms              : {L_rms:.1f} m")
        print(f"  NTU                : {NTU:.4f}")
        print(f"  Q_loss             : {Q_k/1e6:.2f} MW")

    print("\n============================================")
    return Q_loss_tot



Q_loss_total = compute_total_molten_salt_piping_losses(
    T_hot=T_hot,
    T_amb=T_amb,
    rms_per_diam=rms_per_diam,
    lengths_per_diam=lengths_per_diam,
    equivalent_units_per_diam=equivalent_units_per_diam,
    m_dot_unit=m_dot_unit,
    DIAMETER_M=DIAMETER_M,
    U_PIPE=U_PIPE,
    Cp_salt=Cp_salt
)

print(f"\n>>> TOTAL piping losses: {Q_loss_total/1e6:.2f} MW")

N_units = 100
m_dot_total = N_units * m_dot_unit

T_salt_PB = T_hot - Q_loss_total / (m_dot_total * Cp_salt)

print(f">>> Salt temperature at Steam Gen: {T_salt_PB - 273.15:.1f} °C")