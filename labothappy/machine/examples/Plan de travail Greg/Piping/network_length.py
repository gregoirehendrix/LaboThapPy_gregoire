import math
import numpy as np

# =========================================================
# OPTIMIZER
# =========================================================

class CSPOptimizerRows:
    def __init__(self, DIST, CAP):
        self.DIST = DIST
        self.CAP = CAP

    def coord(self, p):
        i, j = p
        return (j * self.DIST, i * self.DIST)

    def manhattan(self, a, b):
        return abs(a[0] - b[0]) + abs(a[1] - b[1])

    def generate_blue_towers(self, nb):
        n_side = math.ceil(math.sqrt(nb))
        M = n_side * 2 - 1
        positions = [(i, j) for i in range(0, M, 2)
                           for j in range(0, M, 2)]

        L = (M - 1) * self.DIST
        cx, cy = L / 2, L / 2
        coords = np.array([self.coord(p) for p in positions])
        dist_center = np.sqrt((coords[:,0] - cx)**2 +
                              (coords[:,1] - cy)**2)

        sel = np.argsort(dist_center)[:nb]
        return [positions[i] for i in sel], M

    def group_by_rows(self, blue_pos):
        ordered = sorted(
            list(enumerate(blue_pos, start=1)),
            key=lambda x: (x[1][0], x[1][1])
        )

        clusters = []
        row = None
        buffer = []

        for idx, (r, _) in ordered:
            if row is None:
                row = r
            if r != row:
                for k in range(0, len(buffer), self.CAP):
                    clusters.append(buffer[k:k+self.CAP])
                buffer = []
                row = r
            buffer.append(idx)

        for k in range(0, len(buffer), self.CAP):
            clusters.append(buffer[k:k+self.CAP])

        return clusters

    def mst_with_hub(self, nodes, coords, hub):
        visited = {hub}
        unvisited = set(nodes) - {hub}
        parent = {}

        while unvisited:
            best = None
            best_d = math.inf
            for u in visited:
                for v in unvisited:
                    d = self.manhattan(coords[u], coords[v])
                    if d < best_d:
                        best_d = d
                        best = (u, v)
            u, v = best
            parent[v] = u
            visited.add(v)
            unvisited.remove(v)

        return parent

    def build_network_for_pb(self, blue_pos, clusters, pb):
        coords = {0: self.coord(pb)}
        for i, pos in enumerate(blue_pos, start=1):
            coords[i] = self.coord(pos)

        parent = {}
        hub_map = {}

        for cluster in clusters:
            hub = min(cluster, key=lambda i: self.manhattan(coords[i], coords[0]))
            mst = self.mst_with_hub(cluster, coords, hub)
            parent.update(mst)
            parent[hub] = 0
            hub_map[hub] = cluster

        return parent, coords, hub_map

    def find_best_pb(self, blue_pos, M):
        clusters = self.group_by_rows(blue_pos)
        blue_set = set(blue_pos)

        best = None
        for i in range(M):
            for j in range(M):
                if (i, j) in blue_set:
                    continue
                parent, coords, hub_map = self.build_network_for_pb(
                    blue_pos, clusters, (i, j)
                )
                cost = sum(self.manhattan(coords[b], coords[p])
                           for b, p in parent.items())
                if best is None or cost < best["cost"]:
                    best = {
                        "parent": parent,
                        "coords": coords,
                        "hubs": hub_map,
                        "cost": cost
                    }
        return best


# =========================================================
# DIAMETER LOGIC
# =========================================================

DIAM_TABLE = {
    1: '1 1/4"',
    2: '2"',
    3: '2 1/2"',
    4: '3"',
    5: '3"'
}

def compute_pipe_info(parent):
    children = {}
    for c, p in parent.items():
        children.setdefault(p, []).append(c)

    charge = {}

    def dfs(n):
        total = 1 if n != 0 else 0
        for ch in children.get(n, []):
            total += dfs(ch)
        charge[n] = total
        return total

    dfs(0)

    pipe_info = {}
    for c, p in parent.items():
        pipe_info[(c, p)] = DIAM_TABLE.get(charge[c], '3"')

    return pipe_info


# =========================================================
# LENGTH PER HUB AND PER DIAMETER (✓ what you asked)
# =========================================================

def compute_lengths_per_hub_and_diameter(sol, pipe_info):
    parent = sol["parent"]
    coords = sol["coords"]
    hubs   = sol["hubs"]

    result = {}

    for hub, nodes in hubs.items():
        per_diam = {}

        # --- inside cluster ---
        for node in nodes:
            cur = node
            while cur != hub:
                par = parent[cur]
                L = abs(coords[cur][0] - coords[par][0]) + \
                    abs(coords[cur][1] - coords[par][1])
                d = pipe_info[(cur, par)]
                per_diam[d] = per_diam.get(d, 0.0) + L
                cur = par

        # --- hub -> PB ---
        par = parent[hub]  # = 0
        L = abs(coords[hub][0] - coords[par][0]) + \
            abs(coords[hub][1] - coords[par][1])
        d = pipe_info[(hub, par)]
        per_diam[d] = per_diam.get(d, 0.0) + L

        result[hub] = per_diam

    return result


# =========================================================
# MAIN
# =========================================================

if __name__ == "__main__":

    nb_units = int(input("Number of CSP units: "))
    CAP      = int(input("Max units per hub: "))
    DIST     = 100.0

    opt = CSPOptimizerRows(DIST, CAP)
    blue_pos, M = opt.generate_blue_towers(nb_units)
    sol = opt.find_best_pb(blue_pos, M)

    pipe_info = compute_pipe_info(sol["parent"])
    hub_lengths = compute_lengths_per_hub_and_diameter(sol, pipe_info)

    print("\n=== PIPING LENGTH PER HUB AND DIAMETER ===")
    for hub in sorted(hub_lengths):
        print(f"\nHub {hub}")
        for diam, L in sorted(hub_lengths[hub].items()):
            print(f"  {diam} : {L:.1f} m")