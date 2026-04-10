# -*- coding: utf-8 -*-
"""
CSP Merged Pipeline
  Step 1 : distribute N units across P power blocks (from layout optimiser)
  Step 2 : for each unique unit count, run piping + thermal loss model
Author: Grégoire Hendrix
"""

import numpy as np
import math


# =============================================================================
# PARAMETERS — edit here
# =============================================================================

INSULATION_MODE     = "rockwool_only"   # "silica+rockwool" or "rockwool_only"
E_INSULATION_TOTAL  = 0.12              # m
E_SILICA            = 0.020             # m (only used in silica+rockwool mode)
E_ROCKWOOL          = E_INSULATION_TOTAL - E_SILICA
LAMBDA_SILICA       = 0.017             # W/m·K
LAMBDA_ROCKWOOL     = 0.040             # W/m·K
H_EXT               = 10.0             # W/m²·K

T_HOT               = 565 + 273.15     # K  — solar field outlet
T_AMB               =  25 + 273.15     # K  — ambient
M_DOT_UNIT          = 5.99             # kg/s per CSP unit

DIAM_TABLE = {1: '1 1/4"', 2: '2"', 3: '2 1/2"', 4: '3"', 5: '3"'}
DIAMETER_M = {
    '1 1/4"': 0.0422,
    '2"'    : 0.0603,
    '2 1/2"': 0.0730,
    '3"'    : 0.0889,
}

DIST_CC = 200   # m — centre-to-centre distance between CSP units

VERBOSE = False  # True → print per-hub pipe details ; False → summary only


# =============================================================================
# 1.  DISTRIBUTION  (from code 1)
# =============================================================================

def calculate_pb_distribution(n_units, n_pb):
    target    = n_units // n_pb
    remainder = n_units % n_pb
    return [target + 1 if i < remainder else target for i in range(n_pb)]


# =============================================================================
# 2.  LAYOUT OPTIMISER  (from code 2)
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
# 3.  PIPE CHARGES & DIAMETERS
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


def compute_length_per_diameter(parent, coords, pipe_info):
    length_per_diam = {}
    for (child, par), info in pipe_info.items():
        diam = info["diameter"]
        x1, y1 = coords[child]
        x2, y2 = coords[par]
        dist   = abs(x1 - x2) + abs(y1 - y2)
        length_per_diam[diam] = length_per_diam.get(diam, 0) + dist
    return length_per_diam


# =============================================================================
# 4.  THERMAL MODEL
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
# 5.  REPORTING
# =============================================================================

def print_thermal_results(pipe_results, children, nb):
    hubs = set(children[0])

    def _all_desc(node):
        r = [node]
        for ch in children[node]:
            r.extend(_all_desc(ch))
        return r

    member_to_hub = {}
    for hub in hubs:
        for m in _all_desc(hub):
            member_to_hub[m] = hub

    hub_pipes = {h: [] for h in hubs}
    for (child, par), info in pipe_results.items():
        h = member_to_hub.get(child)
        if h is not None:
            hub_pipes[h].append((child, par, info))

    total_Q = 0.0
    for hub in sorted(hubs):
        pipes  = sorted(hub_pipes[hub], key=lambda x: x[0])
        hub_Q  = 0.0
        print(f"\n  --- Hub {hub} ---")
        for child, par, info in pipes:
            par_label = "PB" if par == 0 else str(par)
            q_flux = info['Q_loss'] / (math.pi * info['D_total'] * info['L'])
            print(f"    [{child:3d} → {par_label:3s}]  "
                  f"diam={info['diam']:8s}  "
                  f"L={info['L']:6.1f}m  "
                  f"ṁ={info['m_dot']:5.2f}kg/s  "
                  f"T_in={info['T_in']-273.15:6.2f}°C  "
                  f"T_out={info['T_out']-273.15:6.2f}°C  "
                  f"Q={info['Q_loss']/1e3:6.3f}kW  "
                  f"q={q_flux:.1f}W/m²")
            hub_Q += info['Q_loss']
        total_Q += hub_Q
        print(f"    → Hub total loss : {hub_Q/1e3:.3f} kW")

    return total_Q


def print_insulation_summary():
    print(f"\n  Insulation mode : {INSULATION_MODE}")
    print(f"  Total thickness : {E_INSULATION_TOTAL*1e3:.0f} mm  |  h_ext = {H_EXT:.1f} W/m²·K")
    print(f"  {'Diam':10s}  {'D_ext[mm]':10s}  {'D_total[mm]':12s}  {'U_eq[W/m²K]':12s}")
    for diam, D_ext in DIAMETER_M.items():
        e_s, e_r, ls, lr = _resolve_insulation(D_ext)
        U_eq, D_total = compute_U_and_Dtotal(D_ext, e_s, e_r, ls, lr, H_EXT)
        print(f"  {diam:10s}  {D_ext*1e3:10.1f}  {D_total*1e3:12.1f}  {U_eq:12.4f}")


# =============================================================================
# 6.  MAIN
# =============================================================================

if __name__ == "__main__":

    n_units = int(input("Total CSP units        : "))
    n_pb    = int(input("Number of power blocks : "))
    CAP     = int(input("Max units per group    : "))

    # --- Step 1 : distribution ---
    distribution = calculate_pb_distribution(n_units, n_pb)
    unique_counts = sorted(set(distribution))

    print(f"\nDistribution across {n_pb} PBs : {distribution}")
    counts_summary = {u: distribution.count(u) for u in unique_counts}
    for u, cnt in counts_summary.items():
        print(f"  {cnt}x PB with {u} units")

    # --- Step 2 : for each unique unit count ---
    DIST = DIST_CC / 2
    T_C  = T_HOT - 273.15
    Cp_salt = 1443.0 + 0.172 * T_C

    if VERBOSE:
        print_insulation_summary()

    results_by_count = {}

    for nb in unique_counts:
        if VERBOSE:
            print(f"\n{'='*70}")
            print(f"  CASE : {nb} units per PB  "
                  f"(applies to {counts_summary[nb]} power block(s))")
            print(f"{'='*70}")

        opt         = CSPOptimizerRows(DIST, CAP)
        blue_pos, M = opt.generate_blue_towers(nb)
        sol         = opt.find_best_pb(blue_pos, M)

        pipe_info, charge, children = compute_pipe_charges_and_diameters(
            sol["parent"], nb, DIAM_TABLE)

        branch_lengths = [
            abs(sol["coords"][b][0] - sol["coords"][p][0]) +
            abs(sol["coords"][b][1] - sol["coords"][p][1])
            for b, p in sol["parent"].items()
        ]
        rms = math.sqrt(sum(L**2 for L in branch_lengths) / len(branch_lengths))

        if VERBOSE:
            print(f"\n  Total piping length          : {sol['cost']:.1f} m")
            print(f"  Average piping length/unit   : {sol['cost']/nb:.1f} m")
            print(f"  RMS piping length per branch : {rms:.1f} m")

            lengths = compute_length_per_diameter(sol["parent"], sol["coords"], pipe_info)
            print("\n  Length per diameter :")
            for diam, L in lengths.items():
                print(f"    {diam} : {L:.2f} m")

        node_T_out, pipe_results, T_PB = compute_cascaded_temperatures(
            parent     = sol["parent"],
            coords     = sol["coords"],
            pipe_info  = pipe_info,
            children   = children,
            T_hot      = T_HOT,
            T_amb      = T_AMB,
            m_dot_unit = M_DOT_UNIT,
            Cp_salt    = Cp_salt,
        )

        if VERBOSE:
            total_Q = print_thermal_results(pipe_results, children, nb)
            print(f"\n  >>> Salt temperature at PB inlet : {T_PB - 273.15:.3f} °C")
            print(f"      Drop vs solar field outlet   : {T_HOT - T_PB:.3f} K")
            print(f"      Total thermal losses         : {total_Q/1e6:.4f} MW")
        else:
            total_Q = sum(info["Q_loss"] for info in pipe_results.values())

        results_by_count[nb] = dict(T_PB=T_PB, total_Q=total_Q,
                                    piping_length=sol["cost"])

    # --- Final summary ---
    print(f"\n{'='*70}")
    print("  SUMMARY ACROSS ALL POWER BLOCKS")
    print(f"{'='*70}")
    print(f"  {'PB config':20s}  {'Count':6s}  {'T_PB [°C]':12s}  "
          f"{'ΔT [K]':8s}  {'Q_loss [MW]':12s}  {'Piping [m]':10s}")
    for nb, res in results_by_count.items():
        cnt = counts_summary[nb]
        print(f"  {f'{nb} units/PB':20s}  {cnt:6d}  "
              f"{res['T_PB']-273.15:12.3f}  "
              f"{T_HOT - res['T_PB']:8.3f}  "
              f"{res['total_Q']/1e6:12.4f}  "
              f"{res['piping_length']:10.1f}")