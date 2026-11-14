"""Figure 2: raw HTFA, PCA reconstruction, spatial and temporal MS-AA at matched normalized ratio.

Reads the source table that paper_figure02_empirical_PCA_MSAA_benchmark_main.ipynb writes
(outputs/paper_figures/Figure2_source_data.csv) and redraws it with colour = method (the panel
already names the condition), shaded bands (mean +- 1 SEM over decoder repetitions) instead of
error bars, and one shared legend.

Usage (from the repo root, or with MSAA_ROOT set to the run root):
    python paper/main_figures/figure02_plot.py [--source CSV] [--out DIR] [--sharey] [--errorbars]
--sharey puts all three panels on one y-axis, side by side (absolute accuracies comparable across conditions),
with the legend in the empty upper part of the rest panel.
"""
import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(os.environ.get("MSAA_ROOT", Path(__file__).resolve().parents[2]))
ap = argparse.ArgumentParser()
ap.add_argument("--source", default=str(ROOT / "outputs" / "paper_figures" / "Figure2_source_data.csv"))
ap.add_argument("--out", default=str(ROOT / "outputs" / "paper_figures"))
ap.add_argument("--sharey", action="store_true", help="same y-axis in all three panels")
ap.add_argument("--errorbars", action="store_true", help="+-1 SEM (decoder repetitions) error bars instead of shaded bands")
a = ap.parse_args()
OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
SUFFIX = ("_sharey" if a.sharey else "") + ("_errorbars" if a.errorbars else "")

CONDS = [("intact", "Intact narrative"), ("word", "Word-scrambled"), ("rest", "Rest")]
# colour = method; colourblind-checked categorical slots 1-3 (blue, orange, aqua)
METHODS = [("PCA reconstruction", "#2a78d6"), ("Spatial MS-AA", "#eb6834"), ("Temporal MS-AA", "#1baf7a")]
CHANCE = 1 / 300

s = pd.read_csv(a.source)
pca = s[(s.source == "pca") & (s.decode_space == "reconstruction")]
aa = s[s.source == "msaa"]
raw = s[s.source == "raw_htfa"].set_index("condition")
curves = {"PCA reconstruction": pca, "Spatial MS-AA": aa[aa.analysis_type == "spatial"],
          "Temporal MS-AA": aa[aa.analysis_type == "temporal"]}

# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes
# are the printed sizes.
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})
# shared y: three panels side by side, legend in the empty upper part of the rest panel;
# own y per panel: stacked vertically, legend above
if a.sharey:
    fig, axes = plt.subplots(1, 3, figsize=(5.0, 2.6), sharex=True, sharey=True)
else:
    fig, axes = plt.subplots(3, 1, figsize=(5.0, 6.4), sharex=True)
for ax, (cond, title), letter in zip(axes, CONDS, "ABC"):
    for name, color in METHODS:
        x = curves[name][curves[name].condition == cond].sort_values("component_ratio")
        if a.errorbars:
            ax.errorbar(x.component_ratio, x["mean"], yerr=x.err, color=color, lw=1.0, elinewidth=0.7,
                        capsize=1.5, capthick=0.7)
        else:
            ax.fill_between(x.component_ratio, x["mean"] - x.err, x["mean"] + x.err, color=color, alpha=0.25, lw=0)
            ax.plot(x.component_ratio, x["mean"], color=color, lw=1.1)
    r = float(raw.loc[cond, "mean"])
    ax.axhline(r, color="#222222", ls=":", lw=1.0)
    ax.axhline(CHANCE, color="0.6", ls="--", lw=0.6)
    ax.set_title(title, fontsize=7.5, loc="left", fontweight="bold", pad=2)
    ax.text(-0.14 if not a.sharey else -0.10, 1.03, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", lw=0.5); ax.set_axisbelow(True)
    ax.tick_params(labelsize=6.5)
    if not a.sharey or ax is axes[0]:
        ax.set_ylabel("Timepoint-identification\naccuracy", fontsize=7, labelpad=1)
    if not a.sharey and ax is axes[-1]:
        ax.set_xlabel("Normalized representational ratio (components / dimension)", fontsize=7, labelpad=1)
axes[-1].set_xlim(0, 1.0)
# uncertainty (+-1 SEM over decoder repetitions) is stated in the caption, not on the figure
from matplotlib.lines import Line2D
handles = [Line2D([], [], color=c, lw=1.5, label=n) for n, c in METHODS] + [
    Line2D([], [], color="#222222", ls=":", lw=1.0, label="Raw HTFA (700 nodes)"),
    Line2D([], [], color="0.6", ls="--", lw=0.6, label="Chance (1/300)")]
if a.sharey:
    # legend below the panels, clear of every axis (ben, 14 Nov)
    axes[1].set_xlabel("Normalized representational ratio (components / dimension)", fontsize=7, labelpad=2)
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=6.5, handlelength=1.8,
               columnspacing=1.4, labelspacing=0.3, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(pad=0.3, w_pad=0.4, rect=(0, 0.11, 1, 1))
else:
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, 1.035))
    fig.tight_layout(rect=(0, 0.01, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure2_empirical_PCA_MSAA_benchmark{SUFFIX}.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
print("wrote", OUT / f"Figure2_empirical_PCA_MSAA_benchmark{SUFFIX}.png")
