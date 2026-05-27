# -*- coding: utf-8 -*-
"""
CSP Piping — 2D Sweep  nb_units × CAP
  x-axis : CAP        (2 → 10, step 1)   →  9 values
  y-axis : nb_units   (5 → 150, step 5)  → 30 values

Outputs : two heatmaps
  1. Total piping length   [m]
  2. Temperature drop ΔT   [K]

Run from the same directory as full_piping_losses.py, or add its
path to sys.path before importing.

Author: Grégoire Hendrix
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless — change to "Qt5Agg" in Spyder
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

# ── locate full_piping_losses if needed ─────────────────────────────────────
# Uncomment and adapt if the two scripts are not in the same folder:
# sys.path.insert(0, r"C:\path\to\your\scripts")

from full_piping_losses import (
    CSPOptimizerRows,
    compute_pipe_charges_and_diameters,
    compute_cascaded_temperatures,
    DIAM_TABLE,
    T_HOT, T_AMB, M_DOT_UNIT,
    DIST_CC,
)

# ── graphic style ────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"      : "serif",
    "font.size"        : 28,
    "axes.labelsize"   : 28,
    "xtick.labelsize"  : 24,
    "ytick.labelsize"  : 24,
    "legend.fontsize"  : 24,
    "figure.dpi"       : 150,
})

SAVE_DIR = (
    r"C:\Users\gregoire.hendrix@johncockerill.com"
    r"\OneDrive - John Cockerill\Documents\Cockerill"
    r"\Images\3. Thermodynamic analysis\1. Piping"
)
os.makedirs(SAVE_DIR, exist_ok=True)

# ── sweep parameters ─────────────────────────────────────────────────────────
NB_UNITS_RANGE = range(5, 151, 5)   # 5, 10, …, 150  → 30 values
CAP_RANGE      = range(2, 11, 1)    # 2, 3,  …, 10   →  9 values

NB_VALS  = list(NB_UNITS_RANGE)
CAP_VALS = list(CAP_RANGE)

n_nb  = len(NB_VALS)
n_cap = len(CAP_VALS)

# ── result arrays  (rows = nb_units, cols = CAP) ─────────────────────────────
Z_length = np.full((n_nb, n_cap), np.nan)
Z_dT     = np.full((n_nb, n_cap), np.nan)

DIST    = DIST_CC / 2          # m  half the centre-to-centre distance
T_C     = T_HOT - 273.15       # °C  for initial Cp estimate
Cp_salt = 1443.0 + 0.172 * T_C

# ── 2-D sweep ────────────────────────────────────────────────────────────────
t0 = time.time()
print(f"Starting 2-D sweep : {n_nb} × {n_cap} = {n_nb * n_cap} evaluations …\n")

for i, nb in enumerate(NB_VALS):
    for j, cap in enumerate(CAP_VALS):
        try:
            opt         = CSPOptimizerRows(DIST, cap)
            blue_pos, M = opt.generate_blue_towers(nb)
            sol         = opt.find_best_pb(blue_pos, M)

            pipe_info, charge, children = compute_pipe_charges_and_diameters(
                sol["parent"], nb, DIAM_TABLE
            )

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

            Z_length[i, j] = sol["cost"]
            Z_dT[i, j]     = T_HOT - T_PB

        except Exception as exc:
            print(f"  ✗  nb={nb:3d}, CAP={cap}: {exc}")

    elapsed = time.time() - t0
    avg     = elapsed / (i + 1)
    eta     = avg * (n_nb - i - 1)
    print(f"  nb={nb:3d} done   [{i+1:2d}/{n_nb}]  "
          f"elapsed={elapsed:5.1f}s  ETA≈{eta:5.1f}s")

print(f"\nSweep complete in {time.time()-t0:.1f} s\n")


# ── heatmap helper ───────────────────────────────────────────────────────────
def plot_heatmap(Z, nb_vals, cap_vals, cbar_label, fname, cmap="viridis"):
    """
    Render a heatmap with:
      x-axis : CAP  (columns)
      y-axis : nb_units (rows, low → high bottom-up)
    """
    # Step sizes for pixel extents
    dx = (cap_vals[-1] - cap_vals[0]) / (len(cap_vals) - 1)    # = 1
    dy = (nb_vals[-1]  - nb_vals[0])  / (len(nb_vals)  - 1)    # = 5

    extent = [
        cap_vals[0]  - dx / 2,   # left
        cap_vals[-1] + dx / 2,   # right
        nb_vals[0]   - dy / 2,   # bottom
        nb_vals[-1]  + dy / 2,   # top
    ]

    fig, ax = plt.subplots(figsize=(14, 10))

    im = ax.imshow(
        Z,
        origin  = "lower",     # row 0 (nb=5) at bottom
        aspect  = "auto",
        cmap    = cmap,
        extent  = extent,
        interpolation = "nearest",
    )

    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.035)
    cbar.set_label(cbar_label, fontsize=28)
    cbar.ax.tick_params(labelsize=24)

    ax.set_xlabel(r"Maximum units per group  $\mathrm{CAP}$  [–]")
    ax.set_ylabel(r"Number of CSP units  $N$  [–]")

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))

    fig.tight_layout()

    path = os.path.join(SAVE_DIR, fname)
    fig.savefig(path, bbox_inches="tight")
    print(f"  Saved → {path}")
    plt.close(fig)


# ── figure 1 : piping length ─────────────────────────────────────────────────
plot_heatmap(
    Z          = Z_length,
    nb_vals    = NB_VALS,
    cap_vals   = CAP_VALS,
    cbar_label = r"Total piping length  [m]",
    fname      = "heatmap_piping_length.pdf",
    cmap       = "viridis",
)

# ── figure 2 : ΔT ────────────────────────────────────────────────────────────
plot_heatmap(
    Z          = Z_dT,
    nb_vals    = NB_VALS,
    cap_vals   = CAP_VALS,
    cbar_label = r"Temperature drop  $\Delta T$  [K]",
    fname      = "heatmap_dT.pdf",
    cmap       = "plasma",
)

print("\nDone — both heatmaps saved.")