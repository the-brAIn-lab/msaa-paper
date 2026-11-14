"""Supplementary figure for the temporal K=14 analyses (notebooks revision_10a and revision_10c), paper style.

  A  K=14 accumulation: held-out retention of above-chance full-model accuracy for top-m, m = 1..14
     (reconstruction space), intact / word; band = SEM over the 20 outer participant splits. Rest is drawn
     greyed out: its full-model accuracy is near chance, so its retention ratio is not interpreted
  B  generator-anchor masking (+/-2 TR around each archetype's generator TR removed, frozen unmasked
     ranking): above-chance held-out accuracy, unmasked vs masked, for m = 1..14; chance is 1/300
     unmasked and 1/(remaining TRs) masked, so the two curves are compared above their own chance
  C  K=14 motifs at higher K: best map r of each K=14 archetype with any archetype of the fit at K

Reads the outputs of revision_10a / revision_10c (seed-0 fit, every-K held-out splits, 1e-8 tie rule).
Usage: python k14_supplement.py [--in DIR] [--out OUTDIR]
Defaults: outputs/temporal_compact_dictionary_K14_local in the run root, OUTDIR outputs/multiseed_stability/k14_supplement.
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
ap.add_argument("--in", dest="inp", default=str(paths.ROOT / "outputs" / "temporal_compact_dictionary_K14_local"))
ap.add_argument("--out", default=str(paths.OUT / "k14_supplement"))
a = ap.parse_args()
IN, OUT = Path(a.inp), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
CONDS = [("intact", "Intact narrative"), ("word", "Word-scrambled"), ("rest", "Rest")]
COLOR = {"intact": "#7B3294", "word": "#2CA25F", "rest": "#222222"}

acc = pd.read_csv(IN / "K14_accumulation_summary.csv")
acc = acc[(acc.space == "reconstruction")].sort_values("top_m")
cmp_ = pd.read_csv(IN / "K14_anchor_masked_vs_unmasked.csv").sort_values("top_m")
frag = pd.read_csv(IN / "temporal_K14_fragmentation_by_K.csv")


# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes
# are the printed sizes. A and B side by side, C below at full width.
FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 7.5, 7, 6.5, 6
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 1.5,
                     "ytick.major.pad": 1.5, "legend.handlelength": 1.6, "legend.labelspacing": 0.25})


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", lw=0.5); ax.set_axisbelow(True)
    ax.tick_params(labelsize=FS_TICK)


fig = plt.figure(figsize=(5.0, 5.0), constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.03, wspace=0.05, hspace=0.06)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

for c, lab in CONDS:
    x = acc[acc.condition == c]
    if c == "rest":
        # rest full-model accuracy is near chance, so the retention ratio has a near-floor denominator and
        # can exceed 1; it is greyed out and not interpreted (Lucy, 12 Nov)
        axA.plot(x.top_m, x.mean_retention, ":", color="0.72", lw=0.9, marker="o", ms=1.8,
                 label="Rest (near-floor\ndenominator,\nnot interpreted)")
        continue
    axA.fill_between(x.top_m, x.mean_retention - x.mean_outer_sem_retention,
                     x.mean_retention + x.mean_outer_sem_retention, color=COLOR[c], alpha=0.2, lw=0)
    axA.plot(x.top_m, x.mean_retention, "-o", color=COLOR[c], ms=2.2, lw=1.0, label=lab)
axA.axhline(1, color="0.6", ls="--", lw=0.6)
axA.set_xlabel("Top-m (temporal, K = 14)", fontsize=FS_LABEL, labelpad=1)
axA.set_ylabel("Fraction of above-chance\nfull-model accuracy retained", fontsize=FS_LABEL, labelpad=1)
axA.set_xticks(range(1, 15, 2)); axA.set_title("Accumulation over top-m", fontsize=FS_TITLE, loc="left", fontweight="bold")
axA.legend(frameon=False, fontsize=FS_ANNOT, loc="lower right"); style(axA)

for c, lab in CONDS[:2]:
    x = cmp_[cmp_.condition == c]
    um = x.mean_accuracy_unmasked - x.mean_chance_unmasked
    mm = x.mean_accuracy_masked - x.mean_chance_masked
    axB.fill_between(x.top_m, um - x.outer_sem_accuracy_unmasked, um + x.outer_sem_accuracy_unmasked,
                     color=COLOR[c], alpha=0.18, lw=0)
    axB.plot(x.top_m, um, "-", color=COLOR[c], lw=1.0, label=f"{lab.split()[0]}, unmasked")
    axB.fill_between(x.top_m, mm - x.outer_sem_accuracy_masked, mm + x.outer_sem_accuracy_masked,
                     color=COLOR[c], alpha=0.18, lw=0)
    axB.plot(x.top_m, mm, "--", color=COLOR[c], lw=1.0, label=f"{lab.split()[0]}, anchors masked")
    full = x[x.top_m == x.top_m.max()]
    ratio = float(full.above_chance_ratio_masked_to_unmasked.iloc[0])
    axB.annotate(f"{ratio:.0%}\nretained", (14, float(mm.iloc[-1])), xytext=(3, 0), textcoords="offset points",
                 fontsize=FS_ANNOT, color="#222222", va="center", linespacing=0.95)
axB.axhline(0, color="0.6", ls="--", lw=0.6)
axB.set_xlabel("Top-m (temporal, K = 14)", fontsize=FS_LABEL, labelpad=1)
axB.set_ylabel("Held-out timepoint-identification\naccuracy above chance", fontsize=FS_LABEL, labelpad=1)
axB.set_xticks(range(1, 15, 2)); axB.set_xlim(0.5, 18.5)
axB.set_title("Generator-anchor masking", fontsize=FS_TITLE, loc="left", fontweight="bold")
# colour = condition as in A; line style = masking
from matplotlib.lines import Line2D  # noqa: E402
axB.legend(handles=[Line2D([], [], color="0.25", lw=1.0, ls="-", label="Unmasked"),
                    Line2D([], [], color="0.25", lw=1.0, ls="--", label="Anchors masked")],
           frameon=False, fontsize=FS_ANNOT, loc="lower right", bbox_to_anchor=(1.0, 0.06)); style(axB)

g = frag.pivot(index="reference_archetype", columns="target_K", values="best_map_correlation")
im = vheat(axC, g.to_numpy(), cmap="Purples", vmin=0, vmax=1)
for i in range(g.shape[0]):
    for j in range(g.shape[1]):
        v = g.iat[i, j]
        axC.text(j, i, f"{v:.2f}".replace("0.", "."), ha="center", va="center", fontsize=FS_ANNOT,
                 color="white" if v > 0.6 else "#222222")
axC.set_xticks(range(g.shape[1])); axC.set_xticklabels(g.columns, fontsize=FS_TICK)
axC.set_yticks(range(g.shape[0])); axC.set_yticklabels([f"T{k}" for k in g.index], fontsize=FS_TICK)
axC.tick_params(length=0)
axC.set_xlabel("Temporal model order (K)", fontsize=FS_LABEL, labelpad=1); axC.set_ylabel("K = 14 motif", fontsize=FS_LABEL, labelpad=1)
axC.set_title("K = 14 motifs at higher model order", fontsize=FS_TITLE, loc="left", fontweight="bold")
cb = fig.colorbar(im, ax=axC, shrink=0.9, aspect=18, pad=0.015, ticks=[0, 0.5, 1])
cb.set_label("Best map r", fontsize=FS_LABEL, labelpad=2); cb.ax.tick_params(labelsize=FS_TICK, length=2, width=0.6)
cb.outline.set_linewidth(0.6)
for ax, l, dx in zip((axA, axB, axC), "ABC", (-0.30, -0.30, -0.07)):
    ax.text(dx, 1.03, l, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"FigS_temporal_K14_accumulation_anchor_fragmentation.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
print("wrote", OUT)
