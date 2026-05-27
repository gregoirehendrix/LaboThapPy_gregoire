# -*- coding: utf-8 -*-
"""
3-curve extrapolation + average envelope
3 points per curve → quadratic fit → 4th curve = average of the 3.
Author: Grégoire Hendrix
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"      : "serif",
    "font.size"        : 28,
    "axes.labelsize"   : 28,
    "xtick.labelsize"  : 24,
    "ytick.labelsize"  : 24,
    "legend.fontsize"  : 20,
    "figure.dpi"       : 150,
    "axes.grid"        : True,
    "grid.alpha"       : 0.3,
    "grid.linestyle"   : "--",
})

# ════════════════════════════════════════════════════════════════════════════
# ▶▶  ENCODE YOUR 3 CURVES HERE  (3 points each)
# ════════════════════════════════════════════════════════════════════════════
CURVES = [
    {
        "label"  : "Data Set 1",
        "color"  : "#1f77b4",          # blue
        "points" : np.array([
            (100, 97.5),
            (400, 96),
            (700, 90),
        ]),
    },
    {
        "label"  : "Data Set 2",
        "color"  : "sandybrown",          # red
        "points" : np.array([
            (100, 97.4),
            (400, 95.8),
            (550, 93.3),
        ]),
    },
    {
        "label"  : "Data Set 3",
        "color"  : "#2ca02c",          # green
        "points" : np.array([
            (100, 97.5),
            (400, 95.9),
            (700, 88.5),
        ]),
    },
]

# Extrapolation range
X_MIN = 50
X_MAX = 1050

# Key points to annotate
X_KEY_START = 600
X_KEY_STOP  = 1000
X_KEY_STEP  = 50

# Axis labels
X_LABEL = r"Air exit temperature  $T_\mathrm{exit}$  [°C]"
Y_LABEL  = r"Thermal efficiency  $\eta_\mathrm{rec}$  [%]"

# ════════════════════════════════════════════════════════════════════════════
# FIT — quadratic through each set of 3 points
# ════════════════════════════════════════════════════════════════════════════
x_fit = np.linspace(X_MIN, X_MAX, 800)
all_coeffs = []

for curve in CURVES:
    x_pts = curve["points"][:, 0]
    y_pts = curve["points"][:, 1]
    coeffs = np.polyfit(x_pts, y_pts, deg=2)
    curve["poly"]   = np.poly1d(coeffs)
    curve["y_fit"]  = curve["poly"](x_fit)
    all_coeffs.append(coeffs)
    print(f"\n{curve['label']} — polynomial: {curve['poly']}")

# ── 4th curve : average of the 3 quadratics ──────────────────────────────────
# Averaging coefficients is exact (linear operation on polynomials)
mean_coeffs = np.mean(all_coeffs, axis=0)
poly_mean   = np.poly1d(mean_coeffs)
y_mean      = poly_mean(x_fit)
print(f"\nAverage curve — polynomial: {poly_mean}")

# ── Key points on average curve ──────────────────────────────────────────────
x_key = np.arange(X_KEY_START, X_KEY_STOP + 1, X_KEY_STEP)
y_key = poly_mean(x_key)

print(f"\nKey points on average curve ({X_KEY_START}–{X_KEY_STOP} °C, step {X_KEY_STEP}):")
for x, y in zip(x_key, y_key):
    print(f"  T = {x:.0f} °C  →  η = {y:.3f} %")

# ════════════════════════════════════════════════════════════════════════════
# PLOT
# ════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 10))

# Individual curves + their input points
for curve in CURVES:
    ax.plot(x_fit, curve["y_fit"],
            color=curve["color"], lw=2, ls="--", alpha=0.75,
            label=curve["label"])
    #ax.scatter(curve["points"][:, 0], curve["points"][:, 1], color=curve["color"], s=100, zorder=5, marker="o")

# Average (4th) curve
ax.plot(x_fit, y_mean,
        color="darkred", lw=3, ls="-",
        label="Average fit")

# Key points on average curve
ax.scatter(x_key, y_key,
           color="darkred", s=110, zorder=6, marker="D")

for x, y in zip(x_key, y_key):
    ax.annotate(f"{y:.1f}",
                xy=(x-40, y-1), xytext=(0, 11), textcoords="offset points",
                fontsize=16, color="darkred", ha="center")

# Grid
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
ax.grid(True, which="major", alpha=1, linestyle="--")
ax.grid(True, which="minor", alpha=1, linestyle=":")

ax.set_xlabel(X_LABEL)
ax.set_ylabel(Y_LABEL)
ax.set_xlim(X_MIN, X_MAX)
ax.legend(loc="upper right")
fig.tight_layout()
plt.show()