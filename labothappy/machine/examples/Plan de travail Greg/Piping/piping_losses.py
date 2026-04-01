# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 15:29:21 2026

@author: gregoire.hendrix
"""

# =========================================================
# Test standalone — Molten salt piping losses
# =========================================================

import numpy as np

# -------------------------------
# INPUT DATA (TEST CASE)
# -------------------------------

# Temperatures
T_hot = 565 + 273.15      # K  (solar field outlet)
T_amb = 25 + 273.15       # K  (ambient)

# Molten salt properties
T_C = T_hot - 273.15
Cp_salt = 1443.0 + 0.172 * T_C

# Network data (example from your case)
rms_per_diam = {
    '3"'     : 495.0,
    '2 1/2"' : 200.0,
    '2"'     : 200.0,
    '1 1/4"' : 200.0
}

lengths_per_diam = {
    '3"'     : 16000.0,
    '2 1/2"' : 4000.0,
    '2"'     : 4000.0,
    '1 1/4"' : 4000.0
}

# Representative charge (number of towers carried by pipes of this diameter)
charge_per_diam = {
    '3"'     : 100,
    '2 1/2"' : 25,
    '2"'     : 10,
    '1 1/4"' : 5
}

# Flow per tower
m_dot_unit = 0.8           # kg/s per CSP unit

# Pipe diameters [m]
DIAMETER_M = {
    '3"'     : 0.089,
    '2 1/2"' : 0.075,
    '2"'     : 0.060,
    '1 1/4"' : 0.042
}

# Overall heat transfer coefficients [W/m²/K]
U_PIPE = {
    '3"'     : 0.9,
    '2 1/2"' : 1.0,
    '2"'     : 1.1,
    '1 1/4"' : 1.2
}

# -------------------------------
# FUNCTION
# -------------------------------

def compute_total_molten_salt_piping_losses(
    T_hot, T_amb,
    rms_per_diam, lengths_per_diam,
    charge_per_diam, m_dot_unit,
    DIAMETER_M, U_PIPE,
    Cp_salt
):

    Q_loss_tot = 0.0

    print("\n=== Molten salt piping losses per diameter ===")

    for d in rms_per_diam:

        L_rms = rms_per_diam[d]
        L_tot = lengths_per_diam[d]
        n_k   = L_tot / L_rms

        m_dot_k = charge_per_diam[d] * m_dot_unit

        D_k   = DIAMETER_M[d]
        U_k   = U_PIPE[d]

        A_rms = np.pi * D_k * L_rms
        NTU   = U_k * A_rms / (m_dot_k * Cp_salt)

        Q_k = n_k * m_dot_k * Cp_salt * (T_hot - T_amb) * (1 - np.exp(-NTU))

        Q_loss_tot += Q_k

        print(f"\nDiameter {d}")
        print(f"  n_segments   = {n_k:.2f}")
        print(f"  m_dot_k      = {m_dot_k:.2f} kg/s")
        print(f"  L_rms        = {L_rms:.1f} m")
        print(f"  NTU          = {NTU:.3f}")
        print(f"  Q_loss_k     = {Q_k/1e6:.2f} MW")

    print("\n============================================")
    return Q_loss_tot


# -------------------------------
# RUN TEST
# -------------------------------

Q_loss_total = compute_total_molten_salt_piping_losses(
    T_hot=T_hot,
    T_amb=T_amb,
    rms_per_diam=rms_per_diam,
    lengths_per_diam=lengths_per_diam,
    charge_per_diam=charge_per_diam,
    m_dot_unit=m_dot_unit,
    DIAMETER_M=DIAMETER_M,
    U_PIPE=U_PIPE,
    Cp_salt=Cp_salt
)

print(f"\n>>> TOTAL piping losses: {Q_loss_total/1e6:.2f} MW")

# -------------------------------
# Equivalent salt temperature at PB
# -------------------------------

N_units = 100
m_dot_total = N_units * m_dot_unit

T_salt_PB = T_hot - Q_loss_total / (m_dot_total * Cp_salt)

print(f">>> Salt temperature at PB inlet: {T_salt_PB - 273.15:.1f} °C")
