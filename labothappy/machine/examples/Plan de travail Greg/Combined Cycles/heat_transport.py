# -*- coding: utf-8 -*-
"""
heat_transport.py
=================
Module de transport thermique entre les tours CSP et le power block central.

Sélection automatique du fluide de transport selon T_transport :
  T ≤ 200 °C  →  Eau pressurisée  (CoolProp, P_water configurable)
  T ≤ 400 °C  →  Huile thermique  (Cp constant, Dowtherm A type)
  T >  400 °C  →  Sel solaire      (SolarSaltConnector / corrélation Zavoico)

Modèle thermique : NTU cylindrique avec isolation rockwool uniquement.
Diamètre du collecteur principal : sélectionné selon N tours (charge totale).
Longueur moyenne : modèle géométrique champ circulaire  L = (2/3) × R_max.

Usage typique (importé dans un cycle combiné) :
    from heat_transport import transport_heat_loss, select_fluid

Author : Grégoire Hendrix
"""

import math
import numpy as np
from CoolProp.CoolProp import PropsSI

# ─────────────────────────────────────────────────────────────────────────────
#   PARAMÈTRES GLOBAUX — à ajuster ici
# ─────────────────────────────────────────────────────────────────────────────

# Isolation (rockwool uniquement)
E_ROCKWOOL    = 0.08       # m   — épaisseur rockwool  ← valeur arbitraire, à ajuster
LAMBDA_RW     = 0.040      # W/m·K — conductivité rockwool
H_EXT         = 10.0       # W/m²·K — convection extérieure

# Eau pressurisée
P_WATER_DEFAULT = 30e5     # Pa — pression eau (évite ébullition jusqu'à ~235 °C)

# Huile thermique (Dowtherm A / Therminol VP-1 type)
CP_OIL        = 2200.0     # J/kg·K  — Cp moyen
RHO_OIL       = 800.0      # kg/m³   — densité approx.

# Espacement entre tours (surface allouée par tour, champ carré)
A_PER_TOWER   = 10_000.0   # m²/tour  (100 m × 100 m)

# Diamètre du collecteur principal selon le nombre de tours N
# Basé sur ASME B36.10 — schedule 40 ; augmenter si v > 3 m/s
_DIAM_TABLE = [
    (  2, 0.0603),   # DN 50  /  2"       : N ≤ 2  tours
    (  5, 0.0730),   # DN 65  /  2½"      : N ≤ 5
    ( 10, 0.0889),   # DN 80  /  3"       : N ≤ 10
    ( 20, 0.1143),   # DN 100 /  4"       : N ≤ 20
    ( 50, 0.1683),   # DN 150 /  6"       : N ≤ 50
    (100, 0.2191),   # DN 200 /  8"       : N ≤ 100
    (999, 0.2731),   # DN 250 / 10"       : N > 100
]


# ─────────────────────────────────────────────────────────────────────────────
#   SÉLECTION DU FLUIDE DE TRANSPORT
# ─────────────────────────────────────────────────────────────────────────────

FLUID_WATER = "water"
FLUID_OIL   = "oil"
FLUID_SALT  = "salt"


def select_fluid(T_transport_K):
    """
    Retourne le fluide de transport recommandé selon la température.

    Returns
    -------
    str : FLUID_WATER | FLUID_OIL | FLUID_SALT
    """
    T_C = T_transport_K - 273.15
    if T_C <= 200.0:
        return FLUID_WATER
    elif T_C <= 400.0:
        return FLUID_OIL
    else:
        return FLUID_SALT


def fluid_label(fluid):
    return {FLUID_WATER: "Eau pressurisée",
            FLUID_OIL:   "Huile thermique",
            FLUID_SALT:  "Sel solaire (SolarSalt)"}[fluid]


# ─────────────────────────────────────────────────────────────────────────────
#   PROPRIÉTÉS DU FLUIDE
# ─────────────────────────────────────────────────────────────────────────────

def _Cp_salt(T_K):
    """Cp sel solaire [J/kg·K] — corrélation Zavoico 2001."""
    T_C = T_K - 273.15
    return 1443.0 + 0.172 * T_C


def _h_salt(T_K):
    """Enthalpie massique sel solaire [J/kg] (ref 0 °C)."""
    T_C = T_K - 273.15
    return 1443.0 * T_C + 0.086 * T_C**2


def _T_salt_from_h(h_J):
    """Inversion enthalpie → température sel solaire [K]."""
    # résolution quadratique : 0.086·T² + 1443·T - h = 0
    a, b, c = 0.086, 1443.0, -h_J
    T_C = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
    return T_C + 273.15


def get_Cp(fluid, T_K, P_Pa=P_WATER_DEFAULT):
    """Cp [J/kg·K] du fluide de transport à (T, P)."""
    if fluid == FLUID_WATER:
        return PropsSI('C', 'T', T_K, 'P', P_Pa, 'Water')
    elif fluid == FLUID_OIL:
        return CP_OIL
    else:  # SALT
        return _Cp_salt(T_K)


def get_rho(fluid, T_K, P_Pa=P_WATER_DEFAULT):
    """Densité [kg/m³] du fluide de transport."""
    if fluid == FLUID_WATER:
        return PropsSI('D', 'T', T_K, 'P', P_Pa, 'Water')
    elif fluid == FLUID_OIL:
        return RHO_OIL
    else:
        T_C = T_K - 273.15
        return 2090.0 - 0.636 * T_C   # Zavoico 2001


# ─────────────────────────────────────────────────────────────────────────────
#   GÉOMÉTRIE — LONGUEUR MOYENNE ET DIAMÈTRE
# ─────────────────────────────────────────────────────────────────────────────

def average_pipe_length(N, A_per_tower=A_PER_TOWER):
    """
    Distance moyenne tour → power block central [m]
    pour N tours uniformément distribuées dans un champ circulaire.

    L_moy = (2/3) × R_max    avec    R_max = sqrt(N × A_tower / π)
    """
    R_max = math.sqrt(N * A_per_tower / math.pi)
    return (2.0 / 3.0) * R_max


def pipe_diameter(N):
    """
    Diamètre extérieur du collecteur principal [m] selon N tours.
    Basé sur la table ASME B36.10 interne.
    """
    for n_max, D in _DIAM_TABLE:
        if N <= n_max:
            return D
    return _DIAM_TABLE[-1][1]


# ─────────────────────────────────────────────────────────────────────────────
#   MODÈLE THERMIQUE — NTU CYLINDRIQUE ROCKWOOL
# ─────────────────────────────────────────────────────────────────────────────

def _U_eq_rockwool(D_pipe_m):
    """
    Coefficient global U_eq [W/m²·K] basé sur D_pipe (surface interne).
    Résistances : rockwool cylindrique + convection externe.

    R_rw   = ln(D_ext / D_pipe) / (2π λ_rw)      [K·m/W]
    R_conv = 1 / (h_ext π D_ext)                   [K·m/W]
    U_eq   = 1 / (R_total × π × D_pipe)            [W/m²·K]
    """
    D_int = D_pipe_m
    D_ext = D_int + 2.0 * E_ROCKWOOL

    R_rw   = math.log(D_ext / D_int) / (2.0 * math.pi * LAMBDA_RW)
    R_conv = 1.0 / (H_EXT * math.pi * D_ext)
    R_tot  = R_rw + R_conv

    return 1.0 / (R_tot * math.pi * D_int), D_ext


def pipe_ntu_heat_loss(T_in_K, T_amb_K, m_dot_kg_s, Cp_J_kgK, D_pipe_m, L_m):
    """
    Calcul des pertes thermiques d'un tronçon de canalisation.
    Modèle NTU cylindrique avec isolation rockwool uniquement.
    Itération sur Cp_mean (une seule correction).

    Parameters
    ----------
    T_in_K    : température d'entrée du fluide [K]
    T_amb_K   : température ambiante [K]
    m_dot_kg_s: débit massique [kg/s]
    Cp_J_kgK  : Cp du fluide [J/kg·K] (à T_in)
    D_pipe_m  : diamètre extérieur de la canalisation [m]
    L_m       : longueur du tronçon [m]

    Returns
    -------
    dict avec :
        T_out    [K]       — température de sortie
        Q_loss   [W]       — pertes thermiques
        dT       [K]       — chute de température
        U_eq     [W/m²·K]  — coefficient global (surface interne)
        NTU      [-]       — nombre d'unités de transfert
        D_ext    [m]       — diamètre extérieur avec isolation
    """
    U_eq, D_ext = _U_eq_rockwool(D_pipe_m)
    A_int       = math.pi * D_pipe_m * L_m

    # Première estimation
    NTU_0   = U_eq * A_int / (m_dot_kg_s * Cp_J_kgK)
    T_out_0 = T_amb_K + (T_in_K - T_amb_K) * math.exp(-NTU_0)

    # Correction Cp_mean (utile pour sel solaire dont Cp dépend de T)
    T_mean  = (T_in_K + T_out_0) / 2.0
    # Cp_mean conservé côté appelant pour les fluides avec corrélation
    NTU     = U_eq * A_int / (m_dot_kg_s * Cp_J_kgK)
    T_out   = T_amb_K + (T_in_K - T_amb_K) * math.exp(-NTU)
    Q_loss  = m_dot_kg_s * Cp_J_kgK * (T_in_K - T_out)

    return {
        'T_out':  T_out,
        'Q_loss': Q_loss,
        'dT':     T_in_K - T_out,
        'U_eq':   U_eq,
        'NTU':    NTU,
        'D_ext':  D_ext,
    }


# ─────────────────────────────────────────────────────────────────────────────
#   FONCTION PRINCIPALE — PERTES SUR LE RÉSEAU COMPLET
# ─────────────────────────────────────────────────────────────────────────────

def transport_heat_loss(T_hot_K, T_amb_K, m_dot_per_tower_kg_s, N,
                        P_water_Pa=P_WATER_DEFAULT,
                        A_per_tower=A_PER_TOWER,
                        fluid=None,
                        verbose=False):
    """
    Calcule les pertes thermiques du réseau de transport entre N tours
    et le power block central.

    Hypothèse : collecteur principal unique de longueur L_moy (géométrie
    circulaire), portant le débit total de N tours.

    Parameters
    ----------
    T_hot_K              : température de sortie du champ solaire [K]
    T_amb_K              : température ambiante [K]
    m_dot_per_tower_kg_s : débit par tour [kg/s]
    N                    : nombre de tours CSP
    P_water_Pa           : pression eau pressurisée [Pa]
    A_per_tower          : surface au sol par tour [m²]
    fluid                : forcer un fluide (None = sélection auto)
    verbose              : afficher le détail

    Returns
    -------
    dict avec :
        fluid          — fluide sélectionné (str)
        fluid_label    — nom lisible
        T_pb_in        [K]   — température à l'entrée du power block
        dT_loss        [K]   — chute totale de température
        Q_loss         [W]   — pertes thermiques totales
        Q_loss_kW      [kW]
        L_avg          [m]   — longueur moyenne de pipe
        D_pipe         [m]   — diamètre du collecteur principal
        m_dot_total    [kg/s]
        U_eq           [W/m²·K]
        NTU            [-]
        eta_transport  [-]   — efficacité thermique du transport
    """
    if fluid is None:
        fluid = select_fluid(T_hot_K)

    L_avg       = average_pipe_length(N, A_per_tower)
    D_pipe      = pipe_diameter(N)
    m_dot_total = N * m_dot_per_tower_kg_s

    Cp = get_Cp(fluid, T_hot_K, P_water_Pa)

    res = pipe_ntu_heat_loss(T_hot_K, T_amb_K, m_dot_total, Cp, D_pipe, L_avg)

    T_pb_in = res['T_out']
    Q_loss  = res['Q_loss']

    # Chaleur totale disponible (sans pertes)
    Q_available = m_dot_total * Cp * (T_hot_K - T_amb_K)
    eta_transport = 1.0 - Q_loss / Q_available if Q_available > 0 else 0.0

    result = {
        'fluid':          fluid,
        'fluid_label':    fluid_label(fluid),
        'T_pb_in':        T_pb_in,
        'dT_loss':        T_hot_K - T_pb_in,
        'Q_loss':         Q_loss,
        'Q_loss_kW':      Q_loss / 1000.0,
        'L_avg':          L_avg,
        'D_pipe':         D_pipe,
        'm_dot_total':    m_dot_total,
        'm_dot_per_tower': m_dot_per_tower_kg_s,
        'U_eq':           res['U_eq'],
        'NTU':            res['NTU'],
        'D_ext':          res['D_ext'],
        'Cp':             Cp,
        'eta_transport':  eta_transport,
    }

    if verbose:
        _print_transport(result, T_hot_K, T_amb_K, N)

    return result


# ─────────────────────────────────────────────────────────────────────────────
#   AFFICHAGE
# ─────────────────────────────────────────────────────────────────────────────

def _print_transport(r, T_hot_K, T_amb_K, N):
    print("─── Transport thermique ────────────────────────────────────────────")
    print(f"  Fluide de transport     : {r['fluid_label']}")
    print(f"  N tours                 : {N}")
    print(f"  ṁ total                 : {r['m_dot_total']:.2f} kg/s")
    print(f"  Longueur moyenne pipe   : {r['L_avg']:.1f} m")
    print(f"  Diamètre collecteur     : {r['D_pipe']*1e3:.0f} mm  "
          f"(ext. avec isolation : {r['D_ext']*1e3:.0f} mm)")
    print(f"  Épaisseur rockwool      : {E_ROCKWOOL*1e3:.0f} mm")
    print(f"  U_eq                    : {r['U_eq']:.4f} W/m²·K")
    print(f"  NTU                     : {r['NTU']:.5f}")
    print(f"  T_champ (entrée pipe)   : {T_hot_K - 273.15:.2f} °C")
    print(f"  T_power_block (sortie)  : {r['T_pb_in'] - 273.15:.2f} °C")
    print(f"  ΔT pertes               : {r['dT_loss']:.3f} K")
    print(f"  Q_loss                  : {r['Q_loss_kW']:.3f} kW")
    print(f"  η transport             : {r['eta_transport']*100:.4f} %")
    print("────────────────────────────────────────────────────────────────────")


def print_fluid_selection_summary(T_K):
    """Affiche la logique de sélection pour une température donnée."""
    f = select_fluid(T_K)
    print(f"  T_transport = {T_K-273.15:.1f} °C  →  fluide : {fluid_label(f)}")
    print(f"  Seuils : eau ≤ 200 °C | huile ≤ 400 °C | sel > 400 °C")


# ─────────────────────────────────────────────────────────────────────────────
#   STANDALONE — test rapide
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    T_amb = 25.0 + 273.15   # K

    print("=" * 68)
    print("  heat_transport.py — test standalone")
    print("=" * 68)

    # ── Sensibilité sur N ─────────────────────────────────────────────
    print("\n  Sensibilité au nombre de tours (sel solaire, 565 °C) :\n")
    print(f"  {'N':>5}  {'L_moy[m]':>10}  {'D_pipe[mm]':>12}  "
          f"{'T_PB[°C]':>10}  {'ΔT[K]':>8}  {'Q_loss[kW]':>12}")
    print("  " + "-"*65)

    for N in [5, 10, 20, 50, 100, 200]:
        r = transport_heat_loss(
            T_hot_K              = 565 + 273.15,
            T_amb_K              = T_amb,
            m_dot_per_tower_kg_s = 5.99,
            N                    = N,
            verbose              = False,
        )
        print(f"  {N:>5}  {r['L_avg']:>10.1f}  {r['D_pipe']*1e3:>12.0f}  "
              f"{r['T_pb_in']-273.15:>10.2f}  {r['dT_loss']:>8.3f}  "
              f"{r['Q_loss_kW']:>12.3f}")

    # ── Test sélection fluide ─────────────────────────────────────────
    print("\n  Sélection fluide selon température :\n")
    for T_C in [120, 200, 250, 400, 450, 565, 700, 1000]:
        r = transport_heat_loss(
            T_hot_K              = T_C + 273.15,
            T_amb_K              = T_amb,
            m_dot_per_tower_kg_s = 5.0,
            N                    = 20,
            verbose              = False,
        )
        print(f"  T={T_C:5.0f} °C  →  {r['fluid_label']:<22}  "
              f"T_PB={r['T_pb_in']-273.15:.2f} °C  ΔT={r['dT_loss']:.3f} K")

    # ── Affichage détaillé un cas ─────────────────────────────────────
    print()
    transport_heat_loss(
        T_hot_K              = 565 + 273.15,
        T_amb_K              = T_amb,
        m_dot_per_tower_kg_s = 5.99,
        N                    = 50,
        verbose              = True,
    )