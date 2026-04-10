# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 10:48:28 2026

@author: gregoire.hendrix
"""

# -*- coding: utf-8 -*-
"""
Sweep : T_PB en fonction du nombre d'unités CSP
  - n_units varie de 2 à 500
  - 1 power block, CAP = 1 (pas de groupement)
"""

import numpy as np
import math
import matplotlib.pyplot as plt

# =============================================================================
# PARAMETERS
# =============================================================================

INSULATION_MODE     = "rockwool_only"
E_INSULATION_TOTAL  = 0.12
E_SILICA            = 0.020
E_ROCKWOOL          = E_INSULATION_TOTAL - E_SILICA
LAMBDA_SILICA       = 0.017
LAMBDA_ROCKWOOL     = 0.040
H_EXT               = 10.0

T_HOT               = 565 + 273.15
T_AMB               =  25 + 273.15
M_DOT_UNIT          = 5.99

DIAM_TABLE = {1: '1 1/4"', 2: '2"', 3: '2 1/2"', 4: '3"', 5: '3"'}
DIAMETER_M = {
    '1 1/4"': 0.0422,
    '2"'    : 0.0603,
    '2 1/2"': 0.0730,
    '3"'    : 0.0889,
}

DIST_CC = 200   # m centre-to-centre


# =============================================================================
# LAYOUT OPTIMISER
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
# PIPE CHARGES & DIAMETERS
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
# THERMAL MODEL
# =============================================================================

def _resolve_insulation(D_ext_pipe):
    if INSULATION_MODE == "silica+rockwool":
        return E_SILICA, E_ROCKWOOL, LAMBDA_SILICA, LAMBDA_ROCKWOOL
    else:  # rockwool_only
        return 0.0, E_INSULATION_TOTAL, LAMBDA_SILICA, LAMBDA_ROCKWOOL


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
    Cp_mean = 1443.0 + 0.172 * (T_mean - 273.15)

    NTU   = U_eq * A / (m_dot * Cp_mean)
    T_out = T_amb + (T_in - T_amb) * math.exp(-NTU)
    Q_loss = m_dot * Cp_mean * (T_in - T_out)

    return T_out, Q_loss, U_eq, NTU, D_total


def compute_cascaded_temperatures(parent, coords, pipe_info, children,
                                   T_hot, T_amb, m_dot_unit, Cp_salt):
    node_T_out   = {}
    pipe_results = {}

    def _segment(child, par, T_ch):
        m_ch = pipe_info[(child, par)]["charge"] * m_dot_unit
        diam = pipe_info[(child, par)]["diameter"]
        D_k  = DIAMETER_M[diam]
        x1, y1 = coords[child]
        x2, y2 = coords[par]
        L_k  = abs(x1 - x2) + abs(y1 - y2)
        Cp_in = 1443.0 + 0.172 * (T_ch - 273.15)
        e_s, e_r, ls, lr = _resolve_insulation(D_k)
        T_out_k, Q_loss_k, U_k, NTU_k, D_total_k = pipe_heat_loss(
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
# SWEEP
# =============================================================================

def run_single(nb):
    """Run the full pipeline for nb units, 1 PB, CAP=1. Returns T_PB in °C."""
    DIST     = DIST_CC / 2
    CAP      = 1
    T_C      = T_HOT - 273.15
    Cp_salt  = 1443.0 + 0.172 * T_C

    opt         = CSPOptimizerRows(DIST, CAP)
    blue_pos, M = opt.generate_blue_towers(nb)
    sol         = opt.find_best_pb(blue_pos, M)

    pipe_info, charge, children = compute_pipe_charges_and_diameters(
        sol["parent"], nb, DIAM_TABLE)

    _, pipe_results, T_PB = compute_cascaded_temperatures(
        parent     = sol["parent"],
        coords     = sol["coords"],
        pipe_info  = pipe_info,
        children   = children,
        T_hot      = T_HOT,
        T_amb      = T_AMB,
        m_dot_unit = M_DOT_UNIT,
        Cp_salt    = Cp_salt,
    )

    total_Q = sum(info["Q_loss"] for info in pipe_results.values())
    return T_PB - 273.15, total_Q / 1e6   # °C, MW


if __name__ == "__main__":
    n_values = list(range(2, 501))
    T_PB_list = []
    Q_list    = []
    for i, nb in enumerate(n_values):
        T_pb, Q = run_single(nb)
        T_PB_list.append(T_pb)
        Q_list.append(Q)
        if nb % 1 == 0:
            print(f"  n={nb:4d}  T_PB={T_pb:.3f}°C  Q_loss={Q:.4f} MW")

    # --- Plot ---
    fig, ax1 = plt.subplots(figsize=(10, 5))

    color_T = "#1f77b4"
    ax1.plot(n_values, T_PB_list, color=color_T, linewidth=1.5)
    ax1.set_xlabel("Number of CSP units", fontsize=12)
    ax1.set_ylabel("T_PB [°C]", color=color_T, fontsize=12)
    ax1.tick_params(axis="y", labelcolor=color_T)
    ax1.axhline(y=565, color="gray", linestyle="--", linewidth=0.8, label="T_hot = 565 °C")
    ax1.set_xlim(2, 500)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    color_Q = "#d62728"
    ax2.plot(n_values, Q_list, color=color_Q, linewidth=1.2, linestyle="--")
    ax2.set_ylabel("Total thermal losses [MW]", color=color_Q, fontsize=12)
    ax2.tick_params(axis="y", labelcolor=color_Q)

    ax1.set_title("Power block inlet temperature and thermal losses\n"
                  "as a function of the number of CSP units (1 PB, no grouping)",
                  fontsize=12)

    lines1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(lines1, labels1, loc="lower left", fontsize=9)

    plt.tight_layout()
    plt.savefig("sweep_T_PB.png", dpi=150)
    plt.show()
    print("\nPlot saved: sweep_T_PB.png")