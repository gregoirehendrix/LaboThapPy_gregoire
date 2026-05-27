# -*- coding: utf-8 -*-
"""
Corrosion rate comparison – chloride salts vs Solar Salt
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# ── Style (identical to reference scripts) ────────────────────────────────
plt.rcParams.update({
    'font.family'      : 'serif',
    'font.size'        : 28,
    'axes.labelsize'   : 28,
    'xtick.labelsize'  : 20,          # overridden locally: long rotated labels
    'ytick.labelsize'  : 24,
    'legend.fontsize'  : 18,
    'figure.dpi'       : 150,
    'axes.grid'        : True,
    'grid.alpha'       : 0.25,
    'grid.linestyle'   : '--',
    'text.usetex'      : True,
    'text.latex.preamble': r'\usepackage{amsmath}',
})

SAVE_DIR = r"C:\Users\gregoire.hendrix@johncockerill.com\OneDrive - John Cockerill\Documents\Cockerill\Images\1. State of the Art\2. HTF"

# ── Palette (from capex breakdown) ────────────────────────────────────────
C_SOLAR   = '#adb5bd'   # grey   – reference
C_LiNaK   = '#2a9d8f'   # teal   – LiCl-NaCl-KCl
C_NaKCa   = '#457b9d'   # blue   – NaCl-KCl-CaCl2
C_MgNaK   = '#e63946'   # red    – MgCl2-NaCl-KCl
C_NaKZn   = '#6d6875'   # purple – NaCl-KCl-ZnCl2

# ── Data (grouped by salt family, Solar Salt first) ───────────────────────
labels = [
    r"Solar Salt (316SS)" + "\n" + r"565$^\circ$C",
    r"Solar Salt (321SS)" + "\n" + r"565$^\circ$C",
    r"LiCl--NaCl--KCl (20G)" + "\n" + r"700$^\circ$C",
    r"NaCl--KCl--CaCl$_2$ (Haynes230)" + "\n" + r"600$^\circ$C",
    r"NaCl--KCl--CaCl$_2$ (TP347H)" + "\n" + r"600$^\circ$C",
    r"NaCl--KCl--CaCl$_2$ (In625)" + "\n" + r"600$^\circ$C",
    r"MgCl$_2$--NaCl--KCl (SS310)" + "\n" + r"700$^\circ$C",
    r"NaCl--KCl--ZnCl$_2$ (SS316)" + "\n" + r"700$^\circ$C",
]

rates = [
    7.3,        # Solar Salt 316SS
    9.0,        # Solar Salt 321SS
    2350.0,     # LiCl-NaCl-KCl  – Xie et al. 2024
    487.639,    # NaCl-KCl-CaCl2 – Xiao et al. 2025
    2383.628,   # NaCl-KCl-CaCl2 – Xiao et al. 2025
    5437.520,   # NaCl-KCl-CaCl2 – Xiao et al. 2025
    1581.0,     # MgCl2-NaCl-KCl – Ding et al. 2018
    1700.0,     # NaCl-KCl-ZnCl2 – Liu et al. 2019
]

colors = [
    C_SOLAR, C_SOLAR,
    C_LiNaK,
    C_NaKCa, C_NaKCa, C_NaKCa,
    C_MgNaK,
    C_NaKZn,
]

# ── Figure ────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 10))

bars = ax.bar(labels, rates, color=colors, edgecolor='white', linewidth=0.8,
              width=0.6)

# ── Bar annotations ───────────────────────────────────────────────────────
for bar, rate in zip(bars, rates):
    # Solar Salt bars are near-zero: annotate below the top to keep visible
    offset = max(rate * 0.02, 30)
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + offset,
        f'{rate:.0f}',
        ha='center', va='bottom', fontsize=18
    )

# ── Separator line between Solar Salt and chloride salts ──────────────────
ax.axvline(x=1.5, color='black', linewidth=1.2, linestyle=':', alpha=0.5)

# ── Axis labels ───────────────────────────────────────────────────────────
ax.set_ylabel(r'Corrosion rate ($\mu$m/year)')
ax.set_ylim(0, 6400)
ax.tick_params(axis='x', rotation=45)
ax.grid(axis='x', alpha=0)   # vertical gridlines off, keep horizontal

# ── Legend (salt family) ──────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color=C_SOLAR, label='Solar Salt (reference)'),
    mpatches.Patch(color=C_LiNaK, label=r'LiCl--NaCl--KCl'),
    mpatches.Patch(color=C_NaKCa, label=r'NaCl--KCl--CaCl$_2$'),
    mpatches.Patch(color=C_MgNaK, label=r'MgCl$_2$--NaCl--KCl'),
    mpatches.Patch(color=C_NaKZn, label=r'NaCl--KCl--ZnCl$_2$'),
]
ax.legend(handles=legend_handles, loc='upper left', framealpha=0.9)

plt.tight_layout()

fname = 'corrosion_rates_chloride_vs_solar.pdf'
try:
    plt.savefig(os.path.join(SAVE_DIR, fname), format='pdf', bbox_inches='tight')
    print(f"Saved -> {os.path.join(SAVE_DIR, fname)}")
except Exception:
    plt.savefig(fname, format='pdf', bbox_inches='tight')
    print(f"Saved locally as {fname}")

plt.show()