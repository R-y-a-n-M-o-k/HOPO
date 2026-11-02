#!/usr/bin/env python3
"""
Compare vibrational energies (<= EMAX cm^-1) from 2D_vib_energies.xlsx for the
sheets l=1, coupled and cis+trans, separately for A' and A''.

Sheet layout (no header row):
    col 0 = symmetry (A' / A'')   col 1 = state index   col 2 = energy (cm^-1)

Produces:
  1) vib_levels_compare.png  - energy-level diagram: left panel A', right panel A'',
                               one column per sheet, side by side
  3) vib_barcode_all.eps     - single "barcode" diagram: A' group (l=1, coupled, cis+trans)
                               followed by A'' group (same order), E <= EMAX
  2) vib_ordered_compare.png - energy vs. level number (sorted), sheets overlaid,
                               left panel A', right panel A''

Usage:  python plot_vib_energies.py [path/to/2D_vib_energies.xlsx] [output_dir]
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["ps.fonttype"] = 42   # embed TrueType fonts in EPS

XLSX = sys.argv[1] if len(sys.argv) > 1 else "2D_vib_energies.xlsx"
OUTDIR = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
OUTDIR.mkdir(parents=True, exist_ok=True)

EMAX = 8000.0                                   # energy cutoff (cm^-1)
SHEETS = ["full", "l=1", "cis+trans"]        # sheets to compare
SYMS = ["A'", 'A"']
COLORS = dict(zip(SHEETS, ["red", "blue", "black"]))


def load(xlsx):
    """Return {sheet: {sym: sorted energies <= EMAX}}."""
    raw = pd.read_excel(xlsx, sheet_name=SHEETS, header=None)
    data = {}
    for name, df in raw.items():
        sym = df[0].astype(str).str.strip().str.replace("''", '"', regex=False)
        en = pd.to_numeric(df[2], errors="coerce")
        data[name] = {s: en[(sym == s) & (en <= EMAX)].dropna().sort_values().to_numpy()
                      for s in SYMS}
    return data

def plot_barcode(data):
    """Barcode diagram: A' group (l=1, coupled, cis+trans) then A'' group (same order)."""
    n = len(SHEETS)
    gap = 0.8                       # space between the A' and A'' groups
    fig, ax = plt.subplots(figsize=(11, 9))
    ticks, labels = [], []
    for g, s in enumerate(SYMS):
        for i, name in enumerate(SHEETS):
            x0 = g * (n + gap) + i + 0.1
            for e in data[name][s]:
                ax.hlines(e, x0, x0 + 0.8, color=COLORS[name], lw=1.6)
            ticks.append(x0 + 0.4)
            labels.append(name)
        ax.text(g * (n + gap) + n / 2, -0.055, s, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=22, fontweight="bold")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.set_ylim(-100, EMAX + 100)
    ax.set_xlim(0, 2 * n + gap)
    ax.set_ylabel(r"Energy (cm$^{-1}$)", fontsize=20)
    ax.set_title(f"A' and A'' levels, E ≤ {EMAX:.0f} cm$^{{-1}}$", fontsize=20)
    fig.tight_layout()
    fig.savefig(OUTDIR / "vib_barcode_all.eps", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTDIR / "vib_levels_compare.png", dpi=200)


def plot_ordered(data):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, s in zip(axes, SYMS):
        for name in SHEETS:
            e = data[name][s]
            ax.plot(range(1, len(e) + 1), e, "o", ms=3, color=COLORS[name], label=name)
        ax.set_title(f"{s}  (E ≤ {EMAX:.0f} cm$^{{-1}}$)")
        ax.set_xlabel("Level number (in order of energy)")
        ax.grid(alpha=0.3)
        ax.legend()
    axes[0].set_ylabel(r"Energy (cm$^{-1}$)")
    fig.tight_layout()
    fig.savefig(OUTDIR / "vib_ordered_compare.png", dpi=200)


if __name__ == "__main__":
    data = load(XLSX)
    for name, d in data.items():
        print(name, {k: len(v) for k, v in d.items()})
    plot_ordered(data)
    plot_barcode(data)
    print(f"Saved plots to {OUTDIR.resolve()}")