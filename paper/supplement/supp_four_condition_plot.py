"""Supplementary figure S11: four-condition (paragraph-scrambled) sensitivity analysis.

Layout of the original supplement figure (Supplement_four_condition_paragraph_validation): 2 x 2 panels,
one per condition (A Intact narrative, B Paragraph-scrambled, C Word-scrambled, D Rest), each with its
own y-axis. Colour = condition; line style and marker = method (PCA dash-dot triangles, spatial MS-AA
solid squares, temporal MS-AA dashed circles). Raw HTFA is a dotted line in the condition colour and
chance (1/T; T = 272 for paragraph, 300 otherwise) a light grey dotted line. Error bars are +-1 SEM over
the 40 decoder repetitions (stated in the caption, not on the figure). Full-model decoding only, so no
same-data top-m selection enters.

Usage:
    python paper/supplement/supp_four_condition_plot.py [--source DIR] [--out DIR]
DIR holds raw_htfa_summary.csv, pca_summary_with_retention.csv, aa_full_summary_with_retention.csv
(analysis/06_four_condition/rerun_four_condition.py writes them).
"""
import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(os.environ.get("MSAA_ROOT", Path(__file__).resolve().parents[2]))
ap = argparse.ArgumentParser()
ap.add_argument("--source", default=str(ROOT / "four_condition_ratio_decoding_outputs"))
ap.add_argument("--out", default=str(ROOT / "outputs" / "paper_supplement"))
a = ap.parse_args()
SRC, OUT = Path(a.source), Path(a.out)
OUT.mkdir(parents=True, exist_ok=True)

CONDS = [("intact", "Intact narrative", "#7B3294"), ("paragraph", "Paragraph-scrambled", "#E69F00"),
         ("word", "Word-scrambled", "#2CA25F"), ("rest", "Rest", "#222222")]
METHODS = [("PCA", "-.", "^"), ("Spatial MS-AA", "-", "s"), ("Temporal MS-AA", "--", "o")]

raw = pd.read_csv(SRC / "raw_htfa_summary.csv").set_index("condition")
pca = pd.read_csv(SRC / "pca_summary_with_retention.csv")
aa = pd.read_csv(SRC / "aa_full_summary_with_retention.csv")
curves = {"PCA": pca, "Spatial MS-AA": aa[aa.analysis_type == "spatial"],
          "Temporal MS-AA": aa[aa.analysis_type == "temporal"]}

# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes are
# the printed sizes.
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
                     "ytick.major.size": 2.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
                     "axes.titlepad": 3, "axes.labelpad": 1.5})
fig, axes = plt.subplots(2, 2, figsize=(5.0, 4.4), sharex=True, constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.03, h_pad=0.03, wspace=0.06, hspace=0.06)
for ax, (cond, title, color), letter in zip(axes.ravel(), CONDS, "ABCD"):
    ax.axhline(float(raw.loc[cond, "mean"]), color=color, ls=":", lw=1.1)
    ax.axhline(float(raw.loc[cond, "chance"]), color="0.7", ls=":", lw=0.8)
    for name, ls, marker in METHODS:
        x = curves[name][curves[name].condition == cond].sort_values("actual_ratio")
        ax.errorbar(x.actual_ratio, x["mean"], yerr=x.err, color=color, ls=ls, marker=marker, ms=2.4,
                    capsize=1.2, lw=1.0, elinewidth=0.7, capthick=0.7)
    ax.set_title(f"{letter}   {title}", fontweight="bold", loc="left")   # letter in the title: no clash with labels
    if ax in axes[1]:
        ax.set_xlabel("Normalized component ratio")
    if ax in axes[:, 0]:
        ax.set_ylabel("Timepoint-identification accuracy")
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_tick_params(labelbottom=True)

handles = [Line2D([], [], color="0.25", ls=ls, marker=mk, ms=2.4, lw=1.0, label=n) for n, ls, mk in METHODS] + [
    Line2D([], [], color="0.25", ls=":", lw=1.1, label="Raw HTFA"),
    Line2D([], [], color="0.7", ls=":", lw=0.8, label="Chance")]
fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, handlelength=2.2, columnspacing=1.2)
fig.suptitle("Four-condition sensitivity analysis", fontsize=8.5, fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Supplement_four_condition_sensitivity.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
print("wrote", OUT / "Supplement_four_condition_sensitivity.png")
