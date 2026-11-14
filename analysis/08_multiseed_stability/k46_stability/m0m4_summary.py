"""Compact M0-M4 stability summary across the nine further fit seeds (plan: K=46 paper motifs).

For each paper motif M0-M4 and each seed 1-9 (Hungarian match to the seed-0 motif on the 700-node
group-mean maps): spatial map r, temporal profile r (per-scan time-course r, Fisher-z averaged over
the 108 scans) and Yeo-7 composition r, plus recovered = map r >= 0.9. Only the M0-M4 names are
shown; the raw component IDs are in the supplementary crosswalk (make_supp_tables.py).

Reads k46_cross_seed_stability.py's K46_paper_archetypes_alignment.csv.
Usage: python m0m4_summary.py [--alignment CSV] [--out OUTDIR]
Writes M0M4_stability_summary.csv (one row per motif: median and minimum of each r, recovered n/9)
and M0M4_stability_summary.png/.pdf (three motif x seed heatmaps; x = not recovered).
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
ap.add_argument("--alignment", default=str(paths.OUT / "k46_cross_seed_stability" / "K46_paper_archetypes_alignment.csv"))
ap.add_argument("--out", default=str(paths.OUT / "m0m4_summary"))
a = ap.parse_args()
OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
REC = 0.9

al = pd.read_csv(a.alignment)
al = al[al.seed > 0]
MOT = sorted(al.paper_label.unique())
SEEDS = sorted(al.seed.unique())
METRICS = [("map_r", "Spatial map r"), ("timecourse_r", "Temporal profile r"), ("yeo7_profile_r", "Yeo-7 composition r")]

summ = al.groupby("paper_label").agg(
    map_r_median=("map_r", "median"), map_r_min=("map_r", "min"),
    temporal_r_median=("timecourse_r", "median"), temporal_r_min=("timecourse_r", "min"),
    yeo7_r_median=("yeo7_profile_r", "median"), yeo7_r_min=("yeo7_profile_r", "min"),
    n_recovered=("recovered", "sum"), n_seeds=("seed", "nunique"),
    dominant_network=("dominant_net_seed0", "first")).reset_index().rename(columns={"paper_label": "motif"})
summ["recovered"] = summ.n_recovered.astype(int).astype(str) + "/" + summ.n_seeds.astype(str)
summ.to_csv(OUT / "M0M4_stability_summary.csv", index=False)

# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes are
# the printed sizes; the three metrics are stacked so the cells stay wide enough for their values.
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})
fig, axes = plt.subplots(3, 1, figsize=(5.0, 4.1), sharex=True, constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, hspace=0.04)
cmap = plt.get_cmap("Purples")
for ax, (col, title) in zip(axes, METRICS):
    g = al.pivot(index="paper_label", columns="seed", values=col).loc[MOT, SEEDS]
    rec = al.pivot(index="paper_label", columns="seed", values="recovered").loc[MOT, SEEDS].astype(bool)
    im = vheat(ax, g.clip(lower=0).to_numpy(), cmap=cmap, vmin=0, vmax=1)
    for i in range(len(MOT)):
        for j in range(len(SEEDS)):
            v = g.iat[i, j]
            ax.text(j, i, f"{v:.2f}".replace("0.", ".").replace("1.00", "1"), ha="center", va="center", fontsize=6.5,
                    color="white" if v > 0.6 else "#222222")
            if not rec.iat[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#e34a33", lw=1.1))
    ax.set_xticks(range(len(SEEDS))); ax.set_xticklabels(SEEDS, fontsize=6.5)
    ax.set_yticks(range(len(MOT))); ax.set_yticklabels(MOT, fontsize=6.5)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=7.5, fontweight="bold", loc="left", pad=2)
axes[-1].set_xlabel("Fit seed", fontsize=7, labelpad=1)
axes[0].set_yticklabels([f"{m}  ({r})" for m, r in zip(summ.motif, summ.recovered)], fontsize=6.5)
cb = fig.colorbar(im, ax=axes, shrink=0.6, aspect=22, pad=0.02, ticks=[0, 0.5, 1])
cb.set_label("Pearson r with the seed-0 motif", fontsize=7, labelpad=2)
cb.ax.tick_params(labelsize=6.5, length=2, width=0.6); cb.outline.set_linewidth(0.6)
fig.suptitle("Spatial MS-AA, K = 46: motifs M0\u2013M4 across further fit seeds", fontsize=8.5, fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"M0M4_stability_summary.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
print(summ.round(3).to_string(index=False))
