"""Supplement: complete participant-held-out spatial top-m landscape over K x m, 10 seeds.

Rows: held-out top-m accuracy, and the fraction of above-chance full-model accuracy retained.
Columns: intact / word / rest. x = spatial K (10, 20, ..., 700), y = m = 1, 2, 3, 5, 7, 10, 15, 20, 25,
30, 40, 50 (blank where m > K). Across-condition spatial MS-AA basis. Each seed's value averages the two
directions within each outer split, then the 20 splits; cells show the mean over the 10 fit seeds.

m = 1, 2, 3, 5, 10 come from the every-K held-out run; m = 7, 15, ..., 50 from the high-m run with
HM_GRID=landscape (same rankings, splits and decoder seeds; it reproduces m = 10 exactly). Without the
landscape results the figure is drawn with the m <= 10 rows only (marked DRAFT).

Usage: python spatial_topm_landscape.py [--landscape DIR] [--out OUTDIR]
Defaults (paths.py): heldout-everyk-spatial-all-results/, heldout-highm-landscape-results/ in the run root,
OUTDIR outputs/multiseed_stability/spatial_topm_landscape.
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
ap.add_argument("--landscape", default=None)
ap.add_argument("--out", default=None)
a = ap.parse_args()
LAND = Path(a.landscape or paths.ROOT / "heldout-highm-landscape-results")
OUT = Path(a.out or paths.OUT / "spatial_topm_landscape"); OUT.mkdir(parents=True, exist_ok=True)
KS = list(range(10, 701, 10))
MS = [1, 2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50]
CONDS = [("intact", "Intact narrative"), ("word", "Word-scrambled"), ("rest", "Rest")]
COLS = ["condition", "K", "top_m", "fit_seed", "outer_rep", "direction", "mean", "retention_above_chance"]

r = pd.read_csv(paths.heldout_everyk("spatial", "all") / "heldout_topm_runs_all_seeds.csv.gz", usecols=COLS)
r = r[r.K.isin(KS)]
draft = True
f = LAND / "heldout_highm_runs_all_seeds.csv.gz"
if f.exists():
    h = pd.read_csv(f, usecols=COLS + ["analysis_type"])
    h = h[(h.analysis_type == "spatial") & h.K.isin(KS)].drop(columns="analysis_type")
    r = pd.concat([r, h], ignore_index=True); draft = False
r = r[r.top_m.isin(MS) & (r.top_m <= r.K)]
r = r.drop_duplicates(["condition", "K", "top_m", "fit_seed", "outer_rep", "direction"])
seed = (r.groupby(["condition", "K", "top_m", "fit_seed", "outer_rep"])[["mean", "retention_above_chance"]].mean()
        .groupby(["condition", "K", "top_m", "fit_seed"]).mean().reset_index())
summ = (seed.groupby(["condition", "K", "top_m"])
        .agg(accuracy=("mean", "mean"), accuracy_sd=("mean", "std"),
             retention=("retention_above_chance", "mean"), retention_sd=("retention_above_chance", "std"),
             n_seeds=("fit_seed", "nunique")).reset_index())
summ.to_csv(OUT / "spatial_topm_landscape.csv", index=False)


def grid(c, col):
    g = np.full((len(MS), len(KS)), np.nan)
    x = summ[summ.condition == c]
    for row in x.itertuples():
        g[MS.index(int(row.top_m)), KS.index(int(row.K))] = getattr(row, col)
    return g


# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes are
# the printed sizes; K is labelled every 200.
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2, "ytick.major.size": 2, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})
fig, axes = plt.subplots(2, 3, figsize=(5.0, 4.0), sharex=True, sharey=True, constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.03, hspace=0.04)
acc_max = float(np.nanmax(summ.accuracy))
for j, (c, title) in enumerate(CONDS):
    for i, (col, cmap, vmin, vmax, lab) in enumerate([
            ("accuracy", "viridis", 1 / 300, acc_max, "Held-out top-m\ntimepoint-ident. accuracy"),
            ("retention", "magma", 0, 1, "Fraction of above-chance\nfull-model accuracy\nretained")]):
        ax = axes[i, j]
        cm = plt.get_cmap(cmap).copy(); cm.set_bad("0.93")
        im = vheat(ax, grid(c, col), cmap=cm, vmin=vmin, vmax=vmax, origin="lower")
        ax.set_yticks(range(len(MS))); ax.set_yticklabels(MS, fontsize=6)
        xt = [KS.index(k) for k in (200, 400, 600)]
        ax.set_xticks(xt); ax.set_xticklabels([KS[k] for k in xt], fontsize=6.5)
        if j == 0:
            ax.set_ylabel("Top-m", fontsize=7, labelpad=1)
        if i == 1:
            ax.set_xlabel("Spatial model order (K)", fontsize=7, labelpad=1)
        if i == 0:
            ax.set_title(title, fontsize=7.5, fontweight="bold", pad=2)
        if j == 2:
            cb = fig.colorbar(im, ax=axes[i, :], shrink=0.95, aspect=14, pad=0.02)
            cb.set_label(lab, fontsize=6.5, labelpad=2, linespacing=1.0)
            cb.ax.tick_params(labelsize=6, length=2, width=0.6); cb.outline.set_linewidth(0.6)
# details (10 fit seeds, grey = m > K or not yet run) go in the caption
fig.suptitle(("DRAFT, m \u2264 10 only: " if draft else "") + "Spatial MS-AA held-out top-m landscape", fontsize=8.5,
             fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"FigS_spatial_topm_landscape.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
print("wrote", OUT, "| draft" if draft else "")
