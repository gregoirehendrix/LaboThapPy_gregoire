# -*- coding: utf-8 -*-
"""
3-point curve extrapolation — PCHIP
Author: Grégoire Hendrix
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import PchipInterpolator

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"      : "serif",
    "font.size"        : 28,
    "axes.labelsize"   : 28,
    "xtick.labelsize"  : 24,
    "ytick.labelsize"  : 24,
    "legend.fontsize"  : 22,
    "figure.dpi"       : 150,
    "axes.grid"        : True,
    "grid.alpha"       : 0.3,
    "grid.linestyle"   : "--",
})

# ════════════════════════════════════════════════════════════════════════════
# ▶▶  ENCODE YOUR 3 POINTS HERE
# ════════════════════════════════════════════════════════════════════════════
points = np.array([
    (250,  3608),
    (800,  7605),
    (1200, 24381),
])

X_MIN = 250
X_MAX = 1250

X_LABEL = r"Air temperature [°$C$]"
Y_LABEL  = r"Unitary cost [$EUR$]"

# ════════════════════════════════════════════════════════════════════════════
# FIT — PCHIP through the 3 points
# ════════════════════════════════════════════════════════════════════════════
x_pts = points[:, 0]
y_pts = points[:, 1]

pchip = PchipInterpolator(x_pts, y_pts)

x_fit = np.linspace(X_MIN, X_MAX, 500)
y_fit = pchip(x_fit)

# ── Key points every 50°C from 600 to 1000 ───────────────────────────────────
x_key = np.arange(600, 1001, 50)
y_key = pchip(x_key)

print("\nKey points (600–1000 °C, step 50):")
for x, y in zip(x_key, y_key):
    print(f"  T = {x:.0f} °C  →  η = {y:.2f} %")

# ════════════════════════════════════════════════════════════════════════════
# PLOT
# ════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 10))

ax.plot(x_fit, y_fit, color="#1f77b4", lw=2.5, label="Sodeco Valves DN800 (Interpolated)")

ax.scatter(x_key, y_key, color="darkred", s=100, zorder=5,
           marker="D", label=r"Key points (600–1000 °C, $\Delta$T = 50 °C)")

for x, y in zip(x_key, y_key):
    ax.annotate(f"{y:.1f}",
                xy=(x-20, y-1), xytext=(0, 10), textcoords="offset points",
                fontsize=16, color="darkred", ha="center")

ax.scatter(x_pts, y_pts, color="darkolivegreen", s=140, zorder=6,
           marker="o", label="Input points")

for x, y in zip(x_pts, y_pts):
    ax.annotate(f"({x:.0f}, {y:.0f})",
                xy=(x, y), xytext=(12, -12), textcoords="offset points",
                fontsize=16, color="darkolivegreen")
ax.set_xlabel(X_LABEL)
ax.set_ylabel(Y_LABEL)
ax.set_xlim(X_MIN, X_MAX)
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
ax.grid(True, which="major", alpha=1, linestyle="--")
ax.grid(True, which="minor", alpha=1, linestyle=":")
ax.legend(loc="best")
fig.tight_layout()
plt.show()