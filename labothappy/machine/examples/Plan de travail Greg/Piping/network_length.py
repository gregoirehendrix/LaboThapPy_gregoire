import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

# CLASSE OPTIMISATEUR
plt.close("all")
plt.rcParams.update({
    'font.size'        : 28,
    'axes.labelsize'   : 28,
    'xtick.labelsize'  : 24,
    'ytick.labelsize'  : 24,
    'legend.fontsize'  : 22,
})

SAVE_DIR = (r"C:\Users\gregoire.hendrix@johncockerill.com"
            r"\OneDrive - John Cockerill\Documents\Cockerill"
            r"\Images\3. Thermodynamic analysis\1. Piping")
os.makedirs(SAVE_DIR, exist_ok=True)

# Chemin complet du fichier PDF
save_path = os.path.join(SAVE_DIR, "layout.pdf")


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
        cx, cy = L/2, L/2
        coords = np.array([self.coord(p) for p in positions])
        dist_center = np.sqrt((coords[:,0]-cx)**2 + (coords[:,1]-cy)**2)

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

        for idx, pos in ordered:
            r, c = pos

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
        if len(nodes) == 1:
            return {}

        visited = {hub}
        unvisited = set(n for n in nodes if n != hub)
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

        coords = {}
        coords[0] = self.coord(pb)

        for i in range(1, len(blue_pos)+1):
            coords[i] = self.coord(blue_pos[i-1])

        pb_coord = coords[0]

        parent = {}
        cluster_map = {}

        for cluster in clusters:
            hub = min(cluster, key=lambda i: self.manhattan(coords[i], pb_coord))
            mst_parent = self.mst_with_hub(cluster, coords, hub)

            for v, u in mst_parent.items():
                parent[v] = u

            parent[hub] = 0
            cluster_map[hub] = cluster

        total = 0
        for b, p in parent.items():
            total += self.manhattan(coords[b], coords[p])

        return parent, coords, cluster_map, total

    def find_best_pb(self, blue_pos, M):

        clusters = self.group_by_rows(blue_pos)
        blue_set = set(blue_pos)

        candidates = [(i,j) for i in range(M) for j in range(M)
                      if (i,j) not in blue_set]

        best = None

        for pb in candidates:
            parent, coords, info, cost = self.build_network_for_pb(
                blue_pos, clusters, pb
            )
            if best is None or cost < best["cost"]:
                best = dict(pb=pb, parent=parent, coords=coords,
                            clusters=info, cost=cost)

        return best


# CALCUL CHARGES + DIAMETRES + LONGUEURS

DIAM_TABLE = {
    1: '1 1/4"',
    2: '2"',
    3: '2 1/2"',
    4: '3"',
    5: '3"'
}

def compute_pipe_charges_and_diameters(parent, nb_units, diameter_table):
    children = {i: [] for i in range(nb_units+1)}

    for child, par in parent.items():
        children[par].append(child)

    charge = {}

    def dfs(node):
        if node != 0 and len(children[node]) == 0:
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
        c = charge[b]
        if c in diameter_table:
            diam = diameter_table[c]
        else:
            diam = diameter_table[max(diameter_table.keys())]
        pipe_info[(b, p)] = dict(charge=c, diameter=diam)

    return pipe_info, charge


def compute_length_per_diameter(parent, coords, pipe_info):
    length_per_diam = {}

    for (child, par), info in pipe_info.items():
        diam = info["diameter"]
        x1, y1 = coords[child]
        x2, y2 = coords[par]

        dist = abs(x1 - x2) + abs(y1 - y2)
        length_per_diam[diam] = length_per_diam.get(diam, 0) + dist

    return length_per_diam


# LONGUEURS PAR HUB

def compute_hub_pipe_details(parent, coords, pipe_info, clusters, nb_units):
    """
    Pour chaque hub, retourne un dict {diamètre: longueur totale}.
    Inclut toutes les branches internes du cluster + tronçon hub→PB.
    """
    # Construire membre → hub
    member_to_hub = {}
    for hub, members in clusters.items():
        for m in members:
            member_to_hub[m] = hub

    hub_pipes = {hub: {} for hub in clusters}

    for (child, par), info in pipe_info.items():
        diam = info["diameter"]
        x1, y1 = coords[child]
        x2, y2 = coords[par]
        dist = abs(x1 - x2) + abs(y1 - y2)

        if par == 0:
            # Tronçon hub → PB : appartient au hub lui-même
            hub = child
        else:
            # Tronçon interne : appartient au hub du child
            hub = member_to_hub.get(child)

        if hub is not None and hub in hub_pipes:
            hub_pipes[hub][diam] = hub_pipes[hub].get(diam, 0) + dist

    return hub_pipes


def print_hub_pipe_details(hub_pipes, internal_to_visual):
    print("\n=== Longueurs de tuyauterie par hub ===")
    result = []
    for hub_int, diam_dict in hub_pipes.items():
        hub_vis = internal_to_visual[hub_int]
        result.append((hub_vis, diam_dict))

    for hub_vis, diam_dict in sorted(result):
        print(f"\nHub {hub_vis} :")
        for diam, L in sorted(diam_dict.items()):
            print(f"  {diam} : {L:.2f} m")


# AFFICHAGE VISUEL

def make_visual_mapping(blue_pos):
    blue_for_display = sorted(
        [(i, pos) for i, pos in enumerate(blue_pos, start=1)],
        key=lambda x: (-x[1][0], x[1][1])
    )
    internal_to_visual = {
        internal: display
        for display, (internal, _) in enumerate(blue_for_display, start=1)
    }
    return blue_for_display, internal_to_visual


def plot_network(blue_pos, sol, DIST, M):

    parent  = sol["parent"]
    coords  = sol["coords"]
    clusters = sol["clusters"]
    hubs    = set(clusters.keys())

    # --- Translation : PB devient l'origine (0, 0) ---
    xpb_orig, ypb_orig = coords[0]
    def translate(x, y):
        return x - xpb_orig, y - ypb_orig

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect("equal")

    blue_for_display, internal_to_visual = make_visual_mapping(blue_pos)

    for display_id, (internal_id, (r, c)) in enumerate(blue_for_display, start=1):
        x_raw = c * DIST
        y_raw = r * DIST
        x, y  = translate(x_raw, y_raw)

        if internal_id in hubs:
            edge = "#FFD700"
            lw   = 2
        else:
            edge = "white"
            lw   = 0.5

        ax.add_patch(Rectangle((x, y), DIST, DIST,
                                facecolor="#0B3D91", edgecolor=edge, linewidth=lw))
        # Pas de label texte

    # Power block à (0, 0)
    PB_SIZE = 0.3 * DIST
    ax.add_patch(Rectangle(
        ((DIST - PB_SIZE) / 2, (DIST - PB_SIZE) / 2),
        PB_SIZE, PB_SIZE,
        facecolor="red", edgecolor="black", linewidth=2
    ))

    # Connexions
    for b, p in parent.items():
        bx, by = translate(*coords[b])
        px, py = translate(*coords[p])

        bx += DIST / 2;  by += DIST / 2
        px += DIST / 2;  py += DIST / 2

        ax.plot([bx, px], [by, by], color="gray", linewidth=1.2)
        ax.plot([px, px], [by, py], color="gray", linewidth=1.2)

    # Légende
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#0B3D91", edgecolor="white",      label="Solar unit"),
        Patch(facecolor="#0B3D91", edgecolor="#FFD700", linewidth=2, label="Hub unit"),
        Patch(facecolor="red",     edgecolor="black",       label="Power block"),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    ax.set_xlabel("Distance from power block [m]")
    ax.set_ylabel("Distance from power block [m]")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, format="pdf", bbox_inches="tight", dpi=300)
    plt.show()

    return internal_to_visual


def print_clusters_visual(clusters, mapping):

    #print("\n=== Groupings ===")

    result = []
    for hub_int, members_int in clusters.items():
        hub_vis = mapping[hub_int]
        members_vis = [mapping[m] for m in members_int]
        result.append((hub_vis, sorted(members_vis)))

    #for hub, members in sorted(result):
        #print(f"Hub {hub} → {members}")


# MAIN PROGRAMME

if __name__ == "__main__":

    nb = int(input("Amount of CSP units per power block : "))
    CAP = int(input("Maximal amount of units per group : "))
    DIST_2 = 200
    DIST = DIST_2 / 2

    opt = CSPOptimizerRows(DIST, CAP)

    blue_pos, M = opt.generate_blue_towers(nb)
    sol = opt.find_best_pb(blue_pos, M)

    # Longueurs individuelles par branche
    branch_lengths = []
    for b, p in sol["parent"].items():
        bx, by = sol["coords"][b]
        px, py = sol["coords"][p]
        branch_lengths.append(abs(bx - px) + abs(by - py))

    rms = math.sqrt(sum(L**2 for L in branch_lengths) / len(branch_lengths))

    mapping = plot_network(blue_pos, sol, DIST, M)

    # Calcul diamètres et charges
    pipe_info, charge = compute_pipe_charges_and_diameters(sol["parent"], nb, DIAM_TABLE)

    # Détail par hub
    hub_pipes = compute_hub_pipe_details(
        sol["parent"], sol["coords"], pipe_info, sol["clusters"], nb
    )
    print_hub_pipe_details(hub_pipes, mapping)

    print("\nTotal piping length :", sol["cost"], 'm')
    print("Average piping length per unit :", sol["cost"] / nb, 'm')
    #print("RMS piping length per branch   :", rms, 'm')

    print_clusters_visual(sol["clusters"], mapping)

    lengths = compute_length_per_diameter(sol["parent"], sol["coords"], pipe_info)

    print("\n=== Total length per diameter ===")
    for diam, L in lengths.items():
        print(f"{diam} : {L:.2f} m")