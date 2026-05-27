# -*- coding: utf-8 -*-
"""
csp_piping.py
=============
Code exact de l'optimiseur de réseau de piping CSP.
Source : code fourni par Grégoire Hendrix — conservé tel quel.

Ajout : fonction wrapper `network_heat_loss(N, T_hot_K, T_amb_K, m_dot_unit, Cp, fluid_label)`
qui remplace `transport_heat_loss` de heat_transport.py dans combined_progressive.py.
"""

import numpy as np
import math
from CoolProp.CoolProp import PropsSI

# =============================================================================
# PARAMETERS
# =============================================================================

INSULATION_MODE     = "rockwool_only"
E_INSULATION_TOTAL  = 0.08              # m  ← épaisseur rockwool (valeur arbitraire)
E_SILICA            = 0.020             # m  (non utilisé en rockwool_only)
E_ROCKWOOL          = E_INSULATION_TOTAL
LAMBDA_SILICA       = 0.040             # W/m·K
LAMBDA_ROCKWOOL     = 0.040             # W/m·K
H_EXT               = 10.0             # W/m²·K

DIST_CC = 200   # m — centre-to-centre distance entre tours

DIAM_TABLE = {1: '1 1/4"', 2: '2"', 3: '2 1/2"', 4: '3"', 5: '3"'}
DIAMETER_M = {
    '1 1/4"': 0.0422,
    '2"'    : 0.0603,
    '2 1/2"': 0.0730,
    '3"'    : 0.0889,
}

VERBOSE = False


# =============================================================================
# 1. DISTRIBUTION
# =============================================================================

def calculate_pb_distribution(n_units, n_pb):
    target    = n_units // n_pb
    remainder = n_units % n_pb
    return [target + 1 if i < remainder else target for i in range(n_pb)]


# =============================================================================
# 2. LAYOUT OPTIMISER
# =============================================================================

class CSPOptimizerRows:
    def __init__(self, DIST, CAP):
        self.DIST = DIST
        self.CAP  = CAP

    def coord(self, p):
        i, j = p
        return (j * self.DIST, i * self.DIST)

    def manhattan(self, a, b):
        return abs(a[0] - b[0]) + abs(a[1] - b[1])

    def generate_blue_towers(self, nb):
        n_side    = math.ceil(math.sqrt(nb))
        M         = n_side * 2 - 1
        positions = [(i, j) for i in range(0, M, 2)
                             for j in range(0, M, 2)]
        L  = (M - 1) * self.DIST
        cx, cy = L / 2, L / 2
        coords = np.array([self.coord(p) for p in positions])
        dist_c = np.sqrt((coords[:, 0] - cx)**2 + (coords[:, 1] - cy)**2)
        sel    = np.argsort(dist_c)[:nb]
        return [positions[i] for i in sel], M

    def group_by_rows(self, blue_pos):
        ordered = sorted(
            list(enumerate(blue_pos, start=1)),
            key=lambda x: (x[1][0], x[1][1])
        )
        clusters, row, buffer = [], None, []
        for idx, pos in ordered:
            r, c = pos
            if row is None:
                row = r
            if r != row:
                for k in range(0, len(buffer), self.CAP):
                    clusters.append(buffer[k:k + self.CAP])
                buffer, row = [], r
            buffer.append(idx)
        for k in range(0, len(buffer), self.CAP):
            clusters.append(buffer[k:k + self.CAP])
        return clusters

    def mst_with_hub(self, nodes, coords, hub):
        if len(nodes) == 1:
            return {}
        visited, unvisited, parent = {hub}, set(n for n in nodes if n != hub), {}
        while unvisited:
            best, best_d = None, math.inf
            for u in visited:
                for v in unvisited:
                    d = self.manhattan(coords[u], coords[v])
                    if d < best_d:
                        best_d, best = d, (u, v)
            u, v = best
            parent[v] = u
            visited.add(v)
            unvisited.remove(v)
        return parent

    def build_network_for_pb(self, blue_pos, clusters, pb):
        coords   = {0: self.coord(pb)}
        coords.update({i: self.coord(blue_pos[i - 1])
                       for i in range(1, len(blue_pos) + 1)})
        pb_coord = coords[0]
        parent, cluster_map = {}, {}
        for cluster in clusters:
            hub = min(cluster, key=lambda i: self.manhattan(coords[i], pb_coord))
            for v, u in self.mst_with_hub(cluster, coords, hub).items():
                parent[v] = u
            parent[hub]      = 0
            cluster_map[hub] = cluster
        total = sum(self.manhattan(coords[b], coords[p]) for b, p in parent.items())
        return parent, coords, cluster_map, total

    def find_best_pb(self, blue_pos, M):
        clusters = self.group_by_rows(blue_pos)
        blue_set = set(blue_pos)
        candidates = [(i, j) for i in range(M) for j in range(M)
                      if (i, j) not in blue_set]
        best = None
        for pb in candidates:
            parent, coords, info, cost = self.build_network_for_pb(
                blue_pos, clusters, pb)
            if best is None or cost < best["cost"]:
                best = dict(pb=pb, parent=parent, coords=coords,
                            clusters=info, cost=cost)
        return best


# =============================================================================
# 3. PIPE CHARGES & DIAMETERS
# =============================================================================

def compute_pipe_charges_and_diameters(parent, nb_units, diameter_table):
    children = {i: [] for i in range(nb_units + 1)}
    for child, par in parent.items():
        children[par].append(child)

    charge = {}

    def dfs(node):
        if node != 0 and not children[node]:
            charge[node] = 1
            return 1
        total = 1 if node != 0 else 0
        for ch in children[node]:
            total += dfs(ch)
        charge[node] = total
        return total

    dfs(0)

    pipe_info = {}
    for b, p in parent.items():
        c    = charge[b]
        diam = diameter_table.get(c, diameter_table[max(diameter_table)])
        pipe_info[(b, p)] = dict(charge=c, diameter=diam)

    return pipe_info, charge, children


# =============================================================================
# 4. THERMAL MODEL
# =============================================================================

def _resolve_insulation(D_ext_pipe):
    if INSULATION_MODE == "silica+rockwool":
        return E_SILICA, E_ROCKWOOL, LAMBDA_SILICA, LAMBDA_ROCKWOOL
    elif INSULATION_MODE == "rockwool_only":
        return 0.0, E_INSULATION_TOTAL, LAMBDA_SILICA, LAMBDA_ROCKWOOL
    else:
        raise ValueError(f"Unknown INSULATION_MODE: '{INSULATION_MODE}'.")


def compute_U_and_Dtotal(D_ext_pipe, e_silica, e_rockwool,
                          lambda_silica, lambda_rockwool, h_ext):
    D1 = D_ext_pipe
    D2 = D1 + 2 * e_silica
    D3 = D2 + 2 * e_rockwool

    R_silica   = math.log(D2 / D1) / (2 * math.pi * lambda_silica) if e_silica > 0 else 0.0
    R_rockwool = math.log(D3 / D2) / (2 * math.pi * lambda_rockwool)
    R_conv     = 1.0 / (h_ext * math.pi * D3)

    R_total = R_silica + R_rockwool + R_conv
    U_eq    = 1.0 / (R_total * math.pi * D1)
    return U_eq, D3


def pipe_heat_loss(T_in, T_amb, m_dot, Cp, D_ext_pipe, L,
                   e_silica, e_rockwool, lambda_silica, lambda_rockwool, h_ext):
    U_eq, D_total = compute_U_and_Dtotal(
        D_ext_pipe, e_silica, e_rockwool,
        lambda_silica, lambda_rockwool, h_ext)

    A       = math.pi * D_ext_pipe * L
    NTU_0   = U_eq * A / (m_dot * Cp)
    T_out_0 = T_amb + (T_in - T_amb) * math.exp(-NTU_0)
    T_mean  = (T_in + T_out_0) / 2
    Cp_mean = 1443.0 + 0.172 * (T_mean - 273.15)   # SolarSalt par défaut

    NTU   = U_eq * A / (m_dot * Cp_mean)
    T_out = T_amb + (T_in - T_amb) * math.exp(-NTU)
    Q_loss = m_dot * Cp_mean * (T_in - T_out)

    return T_out, Q_loss, U_eq, NTU, D_total


def pipe_heat_loss_fluid(T_in, T_amb, m_dot, Cp, D_ext_pipe, L,
                          e_silica, e_rockwool, lambda_silica, lambda_rockwool, h_ext):
    """
    Variante de pipe_heat_loss qui utilise Cp constant (huile, eau)
    au lieu de la corrélation SolarSalt.
    """
    U_eq, D_total = compute_U_and_Dtotal(
        D_ext_pipe, e_silica, e_rockwool,
        lambda_silica, lambda_rockwool, h_ext)

    A     = math.pi * D_ext_pipe * L
    NTU   = U_eq * A / (m_dot * Cp)
    T_out = T_amb + (T_in - T_amb) * math.exp(-NTU)
    Q_loss = m_dot * Cp * (T_in - T_out)

    return T_out, Q_loss, U_eq, NTU, D_total


def compute_cascaded_temperatures(parent, coords, pipe_info, children,
                                   T_hot, T_amb, m_dot_unit, Cp,
                                   use_salt_cp=False):
    """
    Calcul des températures en cascade dans le réseau.
    use_salt_cp=True  → corrélation SolarSalt (code original)
    use_salt_cp=False → Cp constant (huile, eau)
    """
    node_T_out   = {}
    pipe_results = {}

    def _segment(child, par, T_ch):
        m_ch = pipe_info[(child, par)]["charge"] * m_dot_unit
        diam = pipe_info[(child, par)]["diameter"]
        D_k  = DIAMETER_M[diam]
        x1, y1 = coords[child]
        x2, y2 = coords[par]
        L_k  = abs(x1 - x2) + abs(y1 - y2)
        Cp_in = 1443.0 + 0.172 * (T_ch - 273.15) if use_salt_cp else Cp
        e_s, e_r, ls, lr = _resolve_insulation(D_k)

        if use_salt_cp:
            T_out_k, Q_loss_k, U_k, NTU_k, D_total_k = pipe_heat_loss(
                T_ch, T_amb, m_ch, Cp_in, D_k, L_k, e_s, e_r, ls, lr, H_EXT)
        else:
            T_out_k, Q_loss_k, U_k, NTU_k, D_total_k = pipe_heat_loss_fluid(
                T_ch, T_amb, m_ch, Cp_in, D_k, L_k, e_s, e_r, ls, lr, H_EXT)

        pipe_results[(child, par)] = dict(
            T_in=T_ch, T_out=T_out_k, m_dot=m_ch, L=L_k,
            diam=diam, D_total=D_total_k, U_eq=U_k, NTU=NTU_k, Q_loss=Q_loss_k)
        return T_out_k, m_ch

    def dfs(node):
        if node != 0 and not children[node]:
            node_T_out[node] = T_hot
            return
        for ch in children[node]:
            dfs(ch)
        if node == 0:
            return
        m_total = sum(pipe_info[(ch, node)]["charge"] * m_dot_unit
                      for ch in children[node])
        for ch in children[node]:
            _segment(ch, node, node_T_out[ch])
        T_mix = sum(
            pipe_info[(ch, node)]["charge"] * m_dot_unit *
            pipe_results[(ch, node)]["T_out"]
            for ch in children[node]
        ) / m_total
        node_T_out[node] = T_mix

    for ch in children[0]:
        dfs(ch)

    m_pb_total = sum(pipe_info[(ch, 0)]["charge"] * m_dot_unit
                     for ch in children[0])
    for ch in children[0]:
        _segment(ch, 0, node_T_out[ch])

    T_PB = sum(
        pipe_info[(ch, 0)]["charge"] * m_dot_unit *
        pipe_results[(ch, 0)]["T_out"]
        for ch in children[0]
    ) / m_pb_total

    return node_T_out, pipe_results, T_PB


# =============================================================================
# WRAPPER — remplace transport_heat_loss de heat_transport.py
# =============================================================================

# Propriétés des fluides de transport
CP_OIL = 2200.0   # J/kg·K — huile thermique

def _Cp_salt(T_K):
    T_C = T_K - 273.15
    return 1443.0 + 0.172 * T_C

FLUID_WATER = "water"
FLUID_OIL   = "oil"
FLUID_SALT  = "salt"

def select_fluid(T_K):
    T_C = T_K - 273.15
    if T_C <= 200: return FLUID_WATER
    elif T_C <= 400: return FLUID_OIL
    else: return FLUID_SALT

def fluid_label(fluid):
    return {FLUID_WATER: "Eau pressurisée",
            FLUID_OIL:   "Huile thermique",
            FLUID_SALT:  "Sel solaire (SolarSalt)"}[fluid]


def network_heat_loss(T_hot_K, T_amb_K, m_dot_per_tower_kg_s, N,
                      CAP=5, n_pb=1,
                      P_water=30e5,
                      fluid=None,
                      verbose=False):
    """
    Calcule les pertes thermiques du réseau de transport CSP
    en utilisant l'optimiseur de layout exact (MST, distances Manhattan).

    Remplace transport_heat_loss(T_hot_K, T_amb_K, m_dot_per_tower, N)
    dans combined_progressive.py.

    Parameters
    ----------
    T_hot_K              : température de sortie du champ solaire [K]
    T_amb_K              : température ambiante [K]
    m_dot_per_tower_kg_s : débit par tour [kg/s]
    N                    : nombre de tours CSP
    CAP                  : nombre max de tours par groupe (layout)
    n_pb                 : nombre de power blocks (1 en général)
    P_water              : pression eau si fluide = eau [Pa]
    fluid                : forcer un fluide (None = auto)
    verbose              : affichage détaillé

    Returns
    -------
    dict compatible avec transport_heat_loss :
        fluid, fluid_label, T_pb_in, dT_loss, Q_loss, Q_loss_kW,
        L_avg, D_pipe, m_dot_total, eta_transport
    """
    if fluid is None:
        fluid = select_fluid(T_hot_K)

    # Cp du fluide
    if fluid == FLUID_SALT:
        Cp         = _Cp_salt(T_hot_K)
        use_salt   = True
    elif fluid == FLUID_OIL:
        Cp         = CP_OIL
        use_salt   = False
    else:  # water
        Cp         = PropsSI('C', 'T', T_hot_K, 'P', P_water, 'Water')
        use_salt   = False

    # ── Layout optimiser (code original) ──────────────────────────────
    DIST = DIST_CC / 2
    distribution = calculate_pb_distribution(N, n_pb)
    nb           = distribution[0]   # unités pour le PB principal

    opt         = CSPOptimizerRows(DIST, CAP)
    blue_pos, M = opt.generate_blue_towers(nb)
    sol         = opt.find_best_pb(blue_pos, M)

    pipe_info, charge, children = compute_pipe_charges_and_diameters(
        sol["parent"], nb, DIAM_TABLE)

    # ── Calcul thermique en cascade ────────────────────────────────────
    node_T_out, pipe_results, T_PB = compute_cascaded_temperatures(
        parent     = sol["parent"],
        coords     = sol["coords"],
        pipe_info  = pipe_info,
        children   = children,
        T_hot      = T_hot_K,
        T_amb      = T_amb_K,
        m_dot_unit = m_dot_per_tower_kg_s,
        Cp         = Cp,
        use_salt_cp = use_salt,
    )

    # ── Métriques agrégées ─────────────────────────────────────────────
    total_Q     = sum(info["Q_loss"] for info in pipe_results.values())  # W
    m_dot_total = N * m_dot_per_tower_kg_s

    # Longueur moyenne pondérée par débit
    L_total_weighted = sum(
        info["L"] * info["m_dot"] for info in pipe_results.values()
    )
    L_avg = L_total_weighted / m_dot_total if m_dot_total > 0 else 0

    # Diamètre du tronçon principal (vers le PB)
    D_pipe_main = DIAMETER_M.get(
        pipe_results[list(pipe_results.keys())[0]]["diam"],
        0.089
    )

    # Longueur totale du réseau
    L_network = sol["cost"]

    Q_available = m_dot_total * Cp * (T_hot_K - T_amb_K)
    eta_transport = 1.0 - total_Q / Q_available if Q_available > 0 else 0.0

    result = {
        'fluid':          fluid,
        'fluid_label':    fluid_label(fluid),
        'T_pb_in':        T_PB,
        'dT_loss':        T_hot_K - T_PB,
        'Q_loss':         total_Q,
        'Q_loss_kW':      total_Q / 1000.0,
        'L_avg':          L_avg,
        'L_network':      L_network,
        'D_pipe':         D_pipe_main,
        'm_dot_total':    m_dot_total,
        'm_dot_per_tower': m_dot_per_tower_kg_s,
        'Cp':             Cp,
        'eta_transport':  eta_transport,
        'N':              N,
        'nb_per_pb':      nb,
    }

    if verbose:
        print(f"─── Réseau piping CSP ({'exact MST'}) ───────────────────────────────")
        print(f"  Fluide           : {fluid_label(fluid)}")
        print(f"  N tours          : {N}  (CAP={CAP})")
        print(f"  Longueur réseau  : {L_network:.1f} m  (moy. pondérée = {L_avg:.1f} m)")
        print(f"  T_champ          : {T_hot_K-273.15:.2f} °C")
        print(f"  T_power_block    : {T_PB-273.15:.2f} °C")
        print(f"  ΔT pertes        : {T_hot_K-T_PB:.3f} K")
        print(f"  Q_loss           : {total_Q/1000:.3f} kW")
        print(f"  η transport      : {eta_transport*100:.4f} %")
        print("────────────────────────────────────────────────────────────────────")

    return result


# =============================================================================
# STANDALONE TEST
# =============================================================================

if __name__ == "__main__":
    print("Test network_heat_loss — N=50 tours, sel solaire 565°C\n")
    r = network_heat_loss(
        T_hot_K              = 565 + 273.15,
        T_amb_K              = 25  + 273.15,
        m_dot_per_tower_kg_s = 5.99,
        N                    = 50,
        verbose              = True,
    )
    print(f"\nComparaison formule approchée vs MST exact :")
    print(f"  L_moy approché  : {(2/3)*math.sqrt(50*10000/math.pi):.1f} m")
    print(f"  L_moy MST exact : {r['L_avg']:.1f} m  (réseau total = {r['L_network']:.1f} m)")