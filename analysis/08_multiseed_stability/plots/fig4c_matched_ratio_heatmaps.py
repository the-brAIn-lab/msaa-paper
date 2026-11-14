"""Figure 4C / supplement: spatial minus temporal held-out top-m at exactly matched component ratios, 10 seeds.

Spatial K/700 = temporal K/300 at 15 pairs (EXACT_PAIRS below), for top-m in
1, 2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100 (a cell is blank when m exceeds either K of its pair).
m = 1, 2, 3, 5, 10 come from the every-K held-out run, the other m from the high-m run (same splits,
rankings and decoder seeds; that run reproduces m = 10 and the full model exactly).

Each fit seed's value per cell: average the two directions within each outer split, then the 20 outer
splits. Left: mean over the 10 seeds of (spatial - temporal).
Right: number of seeds (of 10) with spatial > temporal, annotated, on a diverging scale centred at 5
(0 = temporal higher in all seeds, 5 = mixed, 10 = spatial higher in all). The summary table also has
consistency = n_spatial - n_temporal (-10..+10; an exact tie counts for neither; there are none). Seed s spatial is paired with seed s temporal. Fit seeds are initializations, not
independent samples, so this is shown as consistency, without a significance test.

Usage: python fig4c_matched_ratio_heatmaps.py [--highm DIR] [--out OUTDIR]
Defaults (paths.py): heldout-everyk-<analysis>-all-results/, heldout-highm-results/ in the run root,
OUTDIR outputs/multiseed_stability/fig4c_matched_ratio_heatmaps. Writes the figure for all three
conditions, an intact-only panel for the main figure, and the per-seed and summary tables.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.colorbar as _mcb
_mcb.Colorbar.n_rasterize = 10**9  # keep colorbars vector in PDFs
_add_solids = _mcb.Colorbar._add_solids


def _add_solids_no_seams(self, *a, **k):
    # stroke each colour band in its own colour: no hairline gaps in PDF viewers (matplotlib colorbar docs)
    _add_solids(self, *a, **k)
    if getattr(self, "solids", None) is not None:
        self.solids.set_edgecolor("face")
        self.solids.set_linewidth(0.3)


_mcb.Colorbar._add_solids = _add_solids_no_seams


def vheat(ax, arr, origin="upper", **kw):
    """imshow look-alike drawn as vector cells (pcolormesh), so the PDF holds no raster image."""
    a = np.ma.masked_invalid(np.asarray(arr, dtype=float))
    n, m = a.shape
    im = ax.pcolormesh(np.arange(m + 1) - 0.5, np.arange(n + 1) - 0.5, a, shading="flat",
                       edgecolors="face", linewidth=0.3, **kw)  # face-coloured edges: no seams
    ax.set_xlim(-0.5, m - 0.5)
    ax.set_ylim((-0.5, n - 0.5) if origin == "lower" else (n - 0.5, -0.5))
    return im


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--highm", default=None)
ap.add_argument("--out", default=None)
a = ap.parse_args()
HIGHM = Path(a.highm or paths.ROOT / "heldout-highm-results")
OUT = Path(a.out or paths.OUT / "fig4c_matched_ratio_heatmaps"); OUT.mkdir(parents=True, exist_ok=True)

EXACT_PAIRS = [(7, 3), (14, 6), (21, 9), (35, 15), (42, 18), (63, 27), (70, 30), (105, 45), (126, 54),
               (140, 60), (175, 75), (189, 81), (210, 90), (252, 108), (700, 300)]
TOP_M = [1, 2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100]
CONDS = [("intact", "Intact narrative"), ("word", "Word-scrambled"), ("rest", "Rest")]
COLS = ["condition", "K", "top_m", "fit_seed", "outer_rep", "direction", "mean"]


def seed_values(at):
    ks = {p[0] if at == "spatial" else p[1] for p in EXACT_PAIRS}
    r1 = pd.read_csv(paths.heldout_everyk(at, "all") / "heldout_topm_runs_all_seeds.csv.gz", usecols=COLS)
    r2 = pd.read_csv(HIGHM / "heldout_highm_runs_all_seeds.csv.gz", usecols=COLS + ["analysis_type"])
    r2 = r2[r2.analysis_type == at].drop(columns="analysis_type")
    r = pd.concat([r1[r1.K.isin(ks)], r2], ignore_index=True)
    r = r[r.top_m.isin(TOP_M) & (r.top_m <= r.K)]
    # m = 10 is in both runs (identical); at K < requested m the every-K run stores top_m = K (also a duplicate)
    r = r.drop_duplicates(["condition", "K", "top_m", "fit_seed", "outer_rep", "direction"])
    s = (r.groupby(["condition", "K", "top_m", "fit_seed", "outer_rep"])["mean"].mean()
         .groupby(["condition", "K", "top_m", "fit_seed"]).agg(["mean", "count"]).reset_index())
    assert (s["count"] == 20).all(), s[s["count"] != 20].head()
    return s.drop(columns="count")


sp, tp = seed_values("spatial"), seed_values("temporal")
rows = []
for sK, tK in EXACT_PAIRS:
    for c, _ in CONDS:
        for m in TOP_M:
            if m > sK or m > tK:
                continue
            a_ = sp[(sp.condition == c) & (sp.K == sK) & (sp.top_m == m)].set_index("fit_seed")["mean"]
            b_ = tp[(tp.condition == c) & (tp.K == tK) & (tp.top_m == m)].set_index("fit_seed")["mean"]
            j = a_.index.intersection(b_.index)
            assert len(j) == 10, (sK, tK, c, m, len(j))
            for s in j:
                rows.append(dict(condition=c, spatial_K=sK, temporal_K=tK, ratio=sK / 700, top_m=m, fit_seed=int(s),
                                 spatial=a_[s], temporal=b_[s], diff=a_[s] - b_[s]))
seed_tab = pd.DataFrame(rows)
seed_tab.to_csv(OUT / "matched_ratio_per_seed.csv", index=False)
summ = (seed_tab.groupby(["condition", "spatial_K", "temporal_K", "ratio", "top_m"])
        .agg(mean_spatial=("spatial", "mean"), mean_temporal=("temporal", "mean"), mean_diff=("diff", "mean"),
             sd_diff=("diff", "std"), n_seeds_spatial_gt_temporal=("diff", lambda d: int((d > 0).sum())),
             n_seeds_temporal_gt_spatial=("diff", lambda d: int((d < 0).sum())),
             n_seeds_tied=("diff", lambda d: int((d == 0).sum())),
             n_seeds=("diff", "size"), seed0_diff=("diff", "first"))
        .reset_index())
# initialization consistency: seeds with spatial > temporal minus seeds with temporal > spatial (-10..+10).
# An exactly equal seed counts for neither side; with no ties this is 2 * n_spatial - 10.
summ["consistency"] = summ.n_seeds_spatial_gt_temporal - summ.n_seeds_temporal_gt_spatial
summ.to_csv(OUT / "matched_ratio_summary.csv", index=False)
print("exact seed-level ties:", int(summ.n_seeds_tied.sum()), "of", int(summ.n_seeds.sum()))


def grid(c, col):
    g = np.full((len(TOP_M), len(EXACT_PAIRS)), np.nan)
    x = summ[summ.condition == c]
    for _, r in x.iterrows():
        g[TOP_M.index(int(r.top_m)), [p[0] for p in EXACT_PAIRS].index(int(r.spatial_K))] = r[col]
    return g


CONS_LABEL = "Seeds with spatial > temporal (of 10)"
# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so these
# point sizes are the printed sizes. Every other ratio is labelled on x (ratio only); the K pairs are
# listed in the caption.
FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 7.5, 7, 6.5, 6
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2, "ytick.major.size": 2,
                     "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})


def draw(ax_d, ax_n, c, vmax, annotate=True):
    cm_d = plt.get_cmap("RdBu_r").copy(); cm_d.set_bad("0.93")
    im = vheat(ax_d, grid(c, "mean_diff"), cmap=cm_d, vmin=-vmax, vmax=vmax, origin="lower")
    # seed count 0..10 on a diverging scale centred at 5 (mixed): 0 = temporal higher in every seed,
    # 10 = spatial higher in every seed (Lucy, 12 Nov). No exact ties, so this carries the same
    # information as the consistency score in the summary table.
    n = grid(c, "n_seeds_spatial_gt_temporal")
    im2 = vheat(ax_n, n, cmap=cm_d, vmin=0, vmax=10, origin="lower")
    for i in range(n.shape[0] if annotate else 0):
        for j in range(n.shape[1]):
            if np.isfinite(n[i, j]):
                v = int(n[i, j])
                ax_n.text(j, i, str(v), ha="center", va="center", fontsize=FS_ANNOT,
                          color="white" if abs(v - 5) >= 4 else "#222222")
    idx = list(range(0, len(EXACT_PAIRS), 2))
    for ax in (ax_d, ax_n):
        ax.set_yticks(range(len(TOP_M))); ax.set_yticklabels(TOP_M, fontsize=FS_TICK - 0.5)
        ax.set_xticks(idx)
        ax.set_xticklabels([f"{EXACT_PAIRS[j][0] / 700:.2f}".lstrip("0") for j in idx], fontsize=FS_TICK)
    return im, im2


def cbar(fig, im, ax, label, ticks=None, **kw):
    cb = fig.colorbar(im, ax=ax, ticks=ticks, **kw)
    cb.set_label(label, fontsize=FS_LABEL, labelpad=2)
    cb.ax.tick_params(labelsize=FS_TICK, length=2, width=0.6)
    cb.outline.set_linewidth(0.6)
    return cb


vmax = float(np.nanmax(np.abs(summ.mean_diff)))
fig, axes = plt.subplots(3, 2, figsize=(5.0, 6.6), constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.03, wspace=0.04, hspace=0.04)
for (c, title), row in zip(CONDS, axes):
    im, im2 = draw(row[0], row[1], c, vmax)
    row[0].set_ylabel(f"{title}\nTop-m", fontsize=FS_LABEL, labelpad=2)
    row[1].set_ylabel("Top-m", fontsize=FS_LABEL, labelpad=2)
axes[0, 0].set_title("Mean spatial \u2212 temporal", fontsize=FS_TITLE, fontweight="bold")
axes[0, 1].set_title("Seeds with spatial > temporal", fontsize=FS_TITLE, fontweight="bold")
for ax in axes[-1]:
    ax.set_xlabel("Matched component ratio", fontsize=FS_LABEL, labelpad=1)
cbar(fig, im, axes[:, 0], "Held-out accuracy, spatial \u2212 temporal", ticks=[-0.02, 0, 0.02], location="bottom", shrink=0.85, aspect=25, pad=0.01)
cbar(fig, im2, axes[:, 1], CONS_LABEL, ticks=[0, 5, 10], location="bottom", shrink=0.85, aspect=25, pad=0.01)
fig.suptitle("Spatial vs temporal MS-AA at matched component ratios", fontsize=9, fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"FigS_matched_ratio_heatmaps_all_conditions.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

# main-figure panel: intact pair, counts annotated
vmax_i = float(np.nanmax(np.abs(summ[summ.condition == "intact"].mean_diff)))
fig, axes = plt.subplots(1, 2, figsize=(5.0, 2.75), constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.05)
im, im2 = draw(axes[0], axes[1], "intact", vmax_i)
for ax in axes:
    ax.set_xlabel("Matched component ratio", fontsize=FS_LABEL, labelpad=1)
    ax.set_ylabel("Top-m", fontsize=FS_LABEL, labelpad=2)
axes[0].set_title("Intact: mean spatial \u2212 temporal", fontsize=FS_TITLE, fontweight="bold")
axes[1].set_title("Intact: seeds with spatial > temporal", fontsize=FS_TITLE, fontweight="bold")
cbar(fig, im, axes[0], "Held-out accuracy,\nspatial \u2212 temporal", ticks=[-0.02, 0, 0.02], location="bottom", shrink=0.9, aspect=22, pad=0.02)
cbar(fig, im2, axes[1], CONS_LABEL, ticks=[0, 5, 10], location="bottom", shrink=0.9, aspect=22, pad=0.02)
axes[0].text(-0.2, 1.04, "C", transform=axes[0].transAxes, fontsize=10, fontweight="bold", va="bottom")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig4C_matched_ratio_heatmaps_intact.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
print("wrote", OUT, "|", len(summ), "cells")
