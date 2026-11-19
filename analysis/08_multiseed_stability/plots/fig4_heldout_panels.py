"""Figure 4 panels A, B and D and the full-model supplement, from the 10-seed held-out outputs.

  A  spatial held-out top-m (m = 3, 5, 10) vs K, intact / word / rest
  B  temporal held-out top-m (m = 1, 2, 3, 5, 10) vs K, intact / word / rest
  D  temporal matched retained fraction (25, 50, 70 % of K) vs K, intact / word
     (100 % is 1 by definition and is not drawn)
  S  full-model held-out accuracy vs K, spatial and temporal, intact / word / rest
A and B are also drawn in the submitted Fig 4 layout (*_compact: one axis, colour = condition, shade = m;
temporal compact shows m = 1, 3, 10 to stay readable).

Lines are the mean over the 10 fit seeds; ribbons are +-1 SD across seeds (consistency across
initializations, not a sampling error). Titles are kept short; these details go in the captions. Each seed's value is
the mean over 20 outer participant splits, each split averaging its two directions (rank on one
half, score on the other). A cell is drawn only where m < K (at m >= K top-m is the full model).

Usage: python fig4_heldout_panels.py [--mf AGG_DIR] [--out OUTDIR]
Defaults (paths.py): heldout-everyk-<analysis>-all-results/, the matched-fraction aggregate in
outputs/multiseed_stability/matched_fraction/aggregate, OUTDIR outputs/multiseed_stability/fig4_heldout_panels.
Writes one PNG + PDF per panel and the source table behind each.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402

# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes
# below are the printed sizes.
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 1.5,
                     "ytick.major.pad": 1.5, "legend.handlelength": 1.4, "legend.borderaxespad": 0.2,
                     "legend.labelspacing": 0.25})

ap = argparse.ArgumentParser()
ap.add_argument("--mf", default=None)
ap.add_argument("--out", default=None)
ap.add_argument("--mark-k", type=int, default=None,
                help="draw a thin dashed line at this spatial K on panel A (the representative model, e.g. 50)")
a = ap.parse_args()
MF = Path(a.mf or paths.OUT / "matched_fraction" / "aggregate")
OUT = Path(a.out or paths.OUT / "fig4_heldout_panels"); OUT.mkdir(parents=True, exist_ok=True)

CONDS = ["intact", "word", "rest"]
LABEL = {"intact": "Intact narrative", "word": "Word-scrambled", "rest": "Rest"}
COLOR = {"intact": "#7B3294", "word": "#2CA25F", "rest": "#222222"}     # paper condition colours
CHANCE = 1 / 300
XMAX = {"spatial": 700, "temporal": 300}


def shades(hex_color, n):
    """n steps of one hue, light -> dark (sequential: larger m = darker)."""
    c = np.array(matplotlib.colors.to_rgb(hex_color))
    return [tuple(1 - (1 - c) * f) for f in np.linspace(0.35, 1.0, n)]


def style(ax, at):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=6.5)
    ax.grid(axis="y", color="0.9", lw=0.5); ax.set_axisbelow(True)
    ax.set_xlim(0, XMAX[at] * 1.02)
    ax.set_xlabel(f"{at.capitalize()} model order (K)", fontsize=7, labelpad=1)


def save(fig, name, table):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
    table.to_csv(OUT / f"{name}_source.csv", index=False)
    plt.close(fig)
    print("wrote", OUT / name)


def topm_panel(at, ms, name, letter):
    p = pd.read_csv(paths.heldout_everyk(at, "all") / "heldout_per_K.csv")
    p = p[p.top_m.isin(ms) & (p.K > p.top_m)].sort_values(["condition", "top_m", "K"])
    fig, axes = plt.subplots(1, 3, figsize=(5.0, 2.1), sharey=True, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.03)
    for ax, c in zip(axes, CONDS):
        for m, col in zip(ms, shades(COLOR[c], len(ms))):
            x = p[(p.condition == c) & (p.top_m == m)]
            ax.fill_between(x.K, x["mean"] - x.sd, x["mean"] + x.sd, color=col, alpha=0.25, lw=0)
            ax.plot(x.K, x["mean"], color=col, lw=1.0, label=f"top-{m}")
        ax.axhline(CHANCE, color="0.5", ls="--", lw=0.6)
        if at == "spatial" and a.mark_k:   # representative model used for the motif analyses
            ax.axvline(a.mark_k, color="0.6", ls=(0, (3, 2)), lw=0.6, zorder=0.5)
        ax.set_title(LABEL[c], fontsize=7.5, pad=2)
        style(ax, at)
        ax.legend(frameon=False, fontsize=6, loc="upper right")
    axes[0].set_ylabel("Held-out top-m\ntimepoint-identification accuracy", fontsize=7, labelpad=1)
    axes[0].text(-0.42, 1.08, letter, transform=axes[0].transAxes, fontsize=10, fontweight="bold", va="bottom")
    fig.suptitle(f"{at.capitalize()} MS-AA: held-out top-m accuracy", fontsize=8.5, fontweight="bold")
    save(fig, name, p[["analysis_type", "condition", "K", "top_m", "n_seeds", "mean", "sd", "min", "max", "seed0"]])


def topm_compact(at, ms, name, letter):
    """Same data as topm_panel, in the submitted Fig 4 layout: one axis, colour = condition, shade = m."""
    p = pd.read_csv(paths.heldout_everyk(at, "all") / "heldout_per_K.csv")
    p = p[p.top_m.isin(ms) & (p.K > p.top_m)].sort_values(["condition", "top_m", "K"])
    fig, ax = plt.subplots(figsize=(5.0, 3.3))
    for c in CONDS:
        for m, col in zip(ms, shades(COLOR[c], len(ms))):
            x = p[(p.condition == c) & (p.top_m == m)]
            ax.fill_between(x.K, x["mean"] - x.sd, x["mean"] + x.sd, color=col, alpha=0.25, lw=0)
            ax.plot(x.K, x["mean"], color=col, lw=1.0, label=f"{LABEL[c]}, top-{m}")
    ax.axhline(CHANCE, color="0.5", ls="--", lw=0.6)
    style(ax, at)
    ax.set_ylabel("Held-out top-m\ntimepoint-identification accuracy", fontsize=7)
    ax.legend(frameon=False, fontsize=6, ncol=len(CONDS), loc="upper center", bbox_to_anchor=(0.5, -0.16))
    ax.text(-0.12, 1.03, letter, transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax.set_title(f"{at.capitalize()} MS-AA: held-out top-m accuracy", fontsize=8.5, fontweight="bold")
    save(fig, name, p[["analysis_type", "condition", "K", "top_m", "n_seeds", "mean", "sd", "min", "max", "seed0"]])


def full_model_panel():
    rows = []
    for at in ("spatial", "temporal"):
        r = pd.read_csv(paths.heldout_everyk(at, "all") / "heldout_topm_runs_all_seeds.csv.gz",
                        usecols=["condition", "K", "top_m", "outer_rep", "direction", "fit_seed", "full_mean_same_eval"])
        r = r[r.top_m == r.groupby(["condition", "K"]).top_m.transform("min")]   # one row per split x direction
        s = (r.groupby(["condition", "K", "fit_seed", "outer_rep"]).full_mean_same_eval.mean()
             .groupby(["condition", "K", "fit_seed"]).mean())
        g = s.groupby(["condition", "K"]).agg(["mean", "std", "count"]).reset_index()
        g.columns = ["condition", "K", "mean", "sd", "n_seeds"]
        g.insert(0, "analysis_type", at)
        rows.append(g)
    t = pd.concat(rows, ignore_index=True)
    fig, axes = plt.subplots(1, 2, figsize=(5.0, 2.3), sharey=True, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.03)
    for ax, at in zip(axes, ("spatial", "temporal")):
        for c in CONDS:
            x = t[(t.analysis_type == at) & (t.condition == c)].sort_values("K")
            ax.fill_between(x.K, x["mean"] - x.sd, x["mean"] + x.sd, color=COLOR[c], alpha=0.25, lw=0)
            ax.plot(x.K, x["mean"], color=COLOR[c], lw=1.0, label=LABEL[c])
        ax.axhline(CHANCE, color="0.5", ls="--", lw=0.6)
        ax.set_title(f"{at.capitalize()} MS-AA", fontsize=7.5, pad=2)
        style(ax, at)
    axes[0].set_ylabel("Held-out full-model\ntimepoint-identification accuracy", fontsize=7, labelpad=1)
    h, l = axes[0].get_legend_handles_labels()
    h.append(plt.Line2D([], [], color="0.5", ls="--", lw=0.6)); l.append("Chance (1/300)")
    fig.suptitle("Full-model held-out accuracy", fontsize=8.5, fontweight="bold")
    fig.legend(h, l, frameon=False, fontsize=6.5, loc="outside lower center", ncol=4, handlelength=1.8)
    save(fig, "FigS_full_model_accuracy_vs_K", t)


def matched_fraction_panel():
    s = pd.read_csv(MF / "matched_fraction_seed_summary.csv")
    s = s[s.fraction_requested < 1.0]
    g = (s.groupby(["condition", "K", "fraction_requested"])
         .agg(mean_retention=("mean_retention", "mean"), sd_retention=("mean_retention", "std"),
              n_seeds=("fit_seed", "nunique"), top_m=("top_m", "first"), actual_fraction=("actual_fraction", "first"))
         .reset_index())
    fracs = sorted(g.fraction_requested.unique())
    fig, axes = plt.subplots(1, 2, figsize=(5.0, 2.3), sharey=True, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.03)
    for ax, c in zip(axes, ["intact", "word"]):
        for f, col in zip(fracs, shades(COLOR[c], len(fracs))):
            x = g[(g.condition == c) & (g.fraction_requested == f)].sort_values("K")
            # the lightest shade (25 %) gets a denser band so its spread stays visible (ben, 14 Nov)
            ax.fill_between(x.K, x.mean_retention - x.sd_retention, x.mean_retention + x.sd_retention,
                            color=col, alpha=0.45 if f == min(fracs) else 0.25, lw=0)
            ax.plot(x.K, x.mean_retention, "-o", color=col, lw=1.0, ms=2, label=f"{f:.0%} of K")
        ax.set_title(LABEL[c], fontsize=7.5, pad=2)
        ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=6.5)
        ax.grid(axis="y", color="0.9", lw=0.5); ax.set_axisbelow(True)
        ax.set_xlabel("Temporal model order (K)", fontsize=7, labelpad=1)
        ax.legend(frameon=False, fontsize=6, loc="lower right")
    axes[0].set_ylabel("Fraction of above-chance\nfull-model accuracy retained", fontsize=7, labelpad=1)
    axes[0].text(-0.3, 1.06, "B", transform=axes[0].transAxes, fontsize=10, fontweight="bold", va="bottom")
    fig.suptitle("Temporal MS-AA: retained accuracy at matched fraction of K", fontsize=8.5, fontweight="bold")
    save(fig, "Fig4D_temporal_matched_fraction", g)


topm_panel("spatial", [3, 5, 10], "Fig4A_spatial_heldout_topm", "A")
topm_panel("temporal", [1, 2, 3, 5, 10], "Fig4B_temporal_heldout_topm", "B")
topm_compact("spatial", [3, 5, 10], "Fig4A_spatial_heldout_topm_compact", "A")
topm_compact("temporal", [1, 3, 10], "Fig4B_temporal_heldout_topm_compact", "B")
matched_fraction_panel()
full_model_panel()
