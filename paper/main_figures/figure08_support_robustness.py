"""Main Figure 8: support and robustness of motifs M0-M4 at the representative K (panels A-F).

A-D as in paper/main_figures/paper_figure06_support_and_robustness.ipynb (cross-seed recovery, scan-level support,
condition deployment, leave-one-condition-out stability); E: fraction of above-chance full-model accuracy kept
by the held-out top-m subsets at the representative K (10 fit seeds); F: top-5 signal kept after trimming the
first and last 10 TRs (edge-trim run, S10 source).

Plot-only: reads existing CSV outputs, fits and evaluates nothing. Run from the results root (the cwd
build_figures uses):
    python paper/main_figures/figure08_support_robustness.py
"""
import glob
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from msaa_pipeline.paper_archetypes import PAPER_K, paper_archetypes  # noqa: E402

ARCHETYPES, ARCHETYPE_LABELS = paper_archetypes()
LABEL = dict(zip(ARCHETYPES, ARCHETYPE_LABELS))
CONDITIONS = ["intact", "word", "rest"]
COND_LABEL = {"intact": "Intact narrative", "word": "Word-scrambled", "rest": "Rest"}
COLORS = {"intact": "#7B3294", "word": "#2CA25F", "rest": "#222222"}
RECOVERED_R = 0.9
K_REP = PAPER_K
TOP_M = [3, 5, 10]

SUPPORT_DIR = Path("revision_05_archetype_interpretation_and_support")
CI_DIR = Path(f"outputs/interpretation/condition_influence_and_seed_stability_K{K_REP}")
SEED_ALIGN_CSV = Path(os.environ.get(
    "MSAA_SEED_ALIGNMENT",
    f"outputs/multiseed_stability/k{K_REP}_cross_seed_stability/K{K_REP}_paper_archetypes_alignment.csv"))
SUPPORT_CSV = SUPPORT_DIR / "selected_archetype_support_by_independent_participant.csv"
COND_CSV = SUPPORT_DIR / "selected_archetype_condition_summary_independent.csv"
LOCO_CSV = CI_DIR / "leave_one_condition_out_archetype_comparison.csv"
# E: per-seed held-out top-m summaries of the every-K 10-seed run (mean_retention = per-split
# (top-m - chance)/(full K on the same held-out participants - chance), averaged over 20 splits x 2 directions)
HELDOUT_GLOB = os.environ.get(
    "MSAA_HELDOUT_SUMMARIES",
    "heldout-everyk-spatial-all-results/seed*_heldout_topm_summary.csv")
# F: edge-trim summary (S10 source; seed-0 fit at K_REP, held-out top-5, 20 outer splits)
EDGE_CSV = Path(f"msaa_edge_trim_sensitivity_K{K_REP}_top5/edge_trim_accuracy_summary.csv")

OUT = Path(os.environ.get("FIG8_OUT", "outputs/paper_figures"))
OUT.mkdir(parents=True, exist_ok=True)

support = pd.read_csv(SUPPORT_CSV)
cond_summary = pd.read_csv(COND_CSV)
loco = pd.read_csv(LOCO_CSV)
align = pd.read_csv(SEED_ALIGN_CSV)
edge = pd.read_csv(EDGE_CSV).set_index("condition")
seed_files = sorted(glob.glob(HELDOUT_GLOB))
assert len(seed_files) == 10, seed_files
ho = pd.concat([pd.read_csv(f).assign(fit_seed=int(Path(f).name[4:6])) for f in seed_files])

# ---- derived tables (A-D identical to the notebook)
support = support[support["archetype"].isin(ARCHETYPES)].copy()
assert support.groupby("archetype").size().eq(108).all()
align = align[align["seed"] > 0]
SEEDS = sorted(align["seed"].unique())
seed_r = align.pivot(index="paper_label", columns="seed", values="map_r").loc[ARCHETYPE_LABELS, SEEDS]
loco = loco[loco["reference_archetype"].isin(ARCHETYPES)]
loco_r = (loco.assign(dropped=loco["comparison"].str.replace("drop_", "", regex=False))
          .pivot(index="reference_archetype", columns="dropped", values="mean_spatial_map_correlation")
          .reindex(index=ARCHETYPES, columns=CONDITIONS))

ret = (ho[(ho["analysis_type"] == "spatial") & (ho["K"] == K_REP) & (ho["top_m"].isin(TOP_M))]
       .groupby(["condition", "top_m"])["mean_retention"].agg(["mean", "std", "count"]).reset_index())
assert ret["count"].eq(10).all()
ret.to_csv(OUT / f"Figure8E_K{K_REP}_retention_source.csv", index=False)
edge.loc[["intact", "word", "rest"], ["baseline_accuracy", "trim_accuracy", "fraction_preserved",
                                      "fraction_preserved_sem", "n_outer_reps"]].to_csv(OUT / f"Figure8F_K{K_REP}_edge_trim_source.csv")
print(ret.round(4).to_string(index=False))

# ---- figure
import matplotlib.colorbar as _mcb
_mcb.Colorbar.n_rasterize = 10**9
_add_solids = _mcb.Colorbar._add_solids


def _add_solids_no_seams(self, *a, **k):
    _add_solids(self, *a, **k)
    if getattr(self, "solids", None) is not None:
        self.solids.set_edgecolor("face")
        self.solids.set_linewidth(0.3)


_mcb.Colorbar._add_solids = _add_solids_no_seams


def vheat(ax, arr, origin="upper", **kw):
    a = np.ma.masked_invalid(np.asarray(arr, dtype=float))
    n, m = a.shape
    im = ax.pcolormesh(np.arange(m + 1) - 0.5, np.arange(n + 1) - 0.5, a, shading="flat",
                       edgecolors="face", linewidth=0.3, **kw)
    ax.set_xlim(-0.5, m - 0.5)
    ax.set_ylim((-0.5, n - 0.5) if origin == "lower" else (n - 0.5, -0.5))
    return im


FS_SUPTITLE, FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT, FS_LETTER = 9, 8, 7.5, 7, 6, 10
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})
fig = plt.figure(figsize=(4.25, 6.75), constrained_layout=True)  # same width as Fig 8; one extra row
fig.get_layout_engine().set(w_pad=0.03, h_pad=0.04, wspace=0.06, hspace=0.06)
gs = fig.add_gridspec(3, 2, width_ratios=[1.4, 1])
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[2, 0])
axF = fig.add_subplot(gs[2, 1])
x = np.arange(len(ARCHETYPES))
offsets = {"intact": -0.24, "word": 0.0, "rest": 0.24}


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", lw=0.5)
    ax.set_axisbelow(True)


# A
g = seed_r.to_numpy(dtype=float)
im = vheat(axA, np.clip(g, 0, 1), cmap="Purples", vmin=0, vmax=1)
for i in range(g.shape[0]):
    for j in range(g.shape[1]):
        v = g[i, j]
        axA.text(j, i, f"{v:.2f}".replace("0.", ".").replace("1.00", "1"), ha="center", va="center",
                 fontsize=FS_ANNOT, color="white" if v > 0.6 else "#222222")
        if v < RECOVERED_R:
            axA.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#e34a33", lw=1.2))
axA.set_xticks(range(len(SEEDS)))
axA.set_xticklabels(SEEDS)
axA.set_yticks(range(len(ARCHETYPE_LABELS)))
axA.set_yticklabels(ARCHETYPE_LABELS)
axA.tick_params(length=0)
axA.set_xlabel("Fit seed", fontsize=FS_LABEL, labelpad=1)
axA.set_title("Cross-seed motif recovery", fontweight="bold", loc="left", fontsize=FS_TITLE)
cb = fig.colorbar(im, ax=axA, shrink=0.95, aspect=18, pad=0.015, fraction=0.06)
cb.set_label("Spatial-map r\nwith seed-0 motif", fontsize=FS_LABEL, labelpad=2)
cb.set_ticks([0, 0.5, 1])
cb.ax.tick_params(labelsize=FS_TICK)
cb.outline.set_linewidth(0.6)

# B
rng = np.random.default_rng(42)
box_data = []
for pos, k in enumerate(ARCHETYPES):
    cur = support[support["archetype"] == k]
    box_data.append(cur["support_l2"].to_numpy(dtype=float))
    for condition in CONDITIONS:
        vals = cur.loc[cur["condition"] == condition, "support_l2"].to_numpy(dtype=float)
        axB.scatter(np.full(len(vals), pos) + rng.normal(0, 0.06, size=len(vals)), vals, s=3, alpha=0.5,
                    color=COLORS[condition], lw=0)
axB.boxplot(box_data, positions=x, widths=0.5, showfliers=False,
            medianprops=dict(color="#222222", lw=1.0), boxprops=dict(color="0.35", lw=0.7),
            whiskerprops=dict(color="0.35", lw=0.7), capprops=dict(color="0.35", lw=0.7))
axB.set_xticks(x)
axB.set_xticklabels(ARCHETYPE_LABELS)
axB.set_ylabel(r"Spatial support $\|W^{(s)}[k,:]\|_2$", fontsize=FS_LABEL, labelpad=1)
axB.set_title("Broad scan-level support", fontweight="bold", loc="left", fontsize=FS_TITLE - 0.5)
style(axB)

# C
for condition in CONDITIONS:
    cur = cond_summary[cond_summary["condition"] == condition].set_index("archetype").reindex(ARCHETYPES)
    axC.errorbar(x + offsets[condition], cur["mean_support"], yerr=cur["sem_support"], marker="o", ms=3.2,
                 linestyle="none", capsize=1.5, elinewidth=0.8, capthick=0.8, color=COLORS[condition],
                 label=COND_LABEL[condition])
axC.set_xticks(x)
axC.set_xticklabels(ARCHETYPE_LABELS)
axC.set_ylabel(r"Mean spatial support $\|W^{(s)}[k,:]\|_2$", fontsize=FS_LABEL, labelpad=1)
axC.set_title("Condition-specific deployment", fontweight="bold", loc="left", fontsize=FS_TITLE)
style(axC)

# D
for condition in CONDITIONS:
    axD.bar(x + offsets[condition], loco_r[condition], width=0.23, color=COLORS[condition])
axD.axhline(RECOVERED_R, color="#E69F00", ls="--", lw=1.0)
axD.set_xticks(x)
axD.set_xticklabels(ARCHETYPE_LABELS)
axD.set_ylim(0, 1.02)
axD.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axD.set_ylabel("Spatial-map r with\nfull-data motif", fontsize=FS_LABEL, labelpad=1)
axD.set_title("Leave-one-condition-out stability", fontweight="bold", loc="left", fontsize=FS_TITLE - 0.5)
style(axD)

# E - fraction of above-chance full-model accuracy retained by held-out top-m subsets at K_REP
# (mean +- SD over the 10 fit seeds; rest omitted: near-floor denominator)
xe = np.arange(len(TOP_M))
for condition, off in [("intact", -0.18), ("word", 0.18)]:
    cur = ret[ret["condition"] == condition].set_index("top_m").reindex(TOP_M)
    axE.bar(xe + off, cur["mean"], width=0.34, yerr=cur["std"], color=COLORS[condition],
            error_kw=dict(elinewidth=0.8, capsize=1.5, capthick=0.8, ecolor="0.25"))
    for xi, v in zip(xe + off, cur["mean"]):
        axE.text(xi, 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=FS_ANNOT, color="white", rotation=90)
axE.set_xticks(xe)
axE.set_xticklabels([f"{m}\n({m / K_REP:.0%})" for m in TOP_M])
axE.set_xlabel(f"Top-m archetypes kept (% of K = {K_REP})", fontsize=FS_LABEL, labelpad=1)
axE.set_ylim(0, 1.0)
axE.set_ylabel("Fraction of above-chance\nfull-model accuracy", fontsize=FS_LABEL, labelpad=1)
axE.set_title(f"Information retained at K = {K_REP}", fontweight="bold", loc="left", fontsize=FS_TITLE)
style(axE)

# F - top-5 above-chance accuracy kept after removing the first and last 10 TRs (S10B)
xf = np.arange(len(CONDITIONS))
for i, condition in enumerate(CONDITIONS):
    if condition == "rest":
        axF.text(i, 0.03, "n/a", ha="center", va="bottom", fontsize=FS_ANNOT, color="0.35")
        continue
    v = float(edge.loc[condition, "fraction_preserved"])
    s = float(edge.loc[condition, "fraction_preserved_sem"])
    axF.bar(i, v, width=0.6, yerr=s, color=COLORS[condition],
            error_kw=dict(elinewidth=0.8, capsize=1.5, capthick=0.8, ecolor="0.25"))
    axF.text(i, 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=FS_ANNOT, color="white", rotation=90)
axF.axhline(1.0, color="0.5", ls=":", lw=0.7)
axF.set_xticks(xf)
axF.set_xticklabels(["Intact", "Word", "Rest"])
axF.set_xlim(-0.6, len(CONDITIONS) - 0.4)
axF.set_ylim(0, 1.1)
axF.set_ylabel("Top-5 signal kept\nafter edge trimming", fontsize=FS_LABEL, labelpad=1)
axF.set_title("Robust to run-boundary removal", fontweight="bold", loc="left", fontsize=FS_TITLE - 0.5)
style(axF)

h = [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in CONDITIONS] + [
    plt.Line2D([], [], color="#E69F00", ls="--", lw=1.0)]
fig.legend(h, [COND_LABEL[c] for c in CONDITIONS] + ["Recovery threshold (r = 0.9)"], frameon=False, ncol=4,
           loc="upper center", bbox_to_anchor=(0.5, 0.0), fontsize=FS_LABEL, handlelength=1.0, handleheight=0.8,
           handletextpad=0.4, columnspacing=0.9,
           title="Condition (B, C, E, F) / condition left out (D)", title_fontsize=FS_TICK)

for ax, letter in zip([axA, axB, axC, axD, axE, axF], "ABCDEF"):
    ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontsize=FS_LETTER, fontweight="bold", va="bottom")
    ax.tick_params(labelsize=FS_TICK)
axA.tick_params(axis="y", labelsize=FS_LABEL)

fig.suptitle("Support and robustness of motifs M0–M4", fontsize=FS_SUPTITLE, fontweight="bold")
png = OUT / "Figure8_motif_support_and_robustness.png"
pdf = OUT / "Figure8_motif_support_and_robustness.pdf"
fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.02)
fig.savefig(pdf, dpi=600, bbox_inches="tight", pad_inches=0.02)
print("Saved:", png)
print("Saved:", pdf)
