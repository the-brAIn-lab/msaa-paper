"""Supplementary figures for the K=46 paper motifs M0-M4, in paper style.

  atlas audit  Yeo-7 composition of each motif with all 700 HTFA nodes, only nodes whose centre lies
               inside the cortical atlas, and only nodes assigned with confidence >= 0.5, plus the
               node-level audit numbers (Figure 5 notebook outputs).
  conditions   how strongly each motif is deployed in each condition (per-scan rank-1 contribution
               ||X_i c s_i||_F, 36 scans per condition) and how well each motif is recovered when a
               whole condition is left out of the fit (map r with the full-data motif; revision_09 outputs).

Usage: python paper/supplement/supp_atlas_condition_plot.py [--root RUN_ROOT] [--out DIR]
RUN_ROOT holds paper_figure_outputs/ and outputs/interpretation/ (defaults to MSAA_ROOT or the repo root).
"""
import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("MSAA_ROOT", str(Path(__file__).resolve().parents[2])))
ap.add_argument("--out", default=None)
a = ap.parse_args()
ROOT = Path(a.root)
OUT = Path(a.out or ROOT / "outputs" / "paper_supplement"); OUT.mkdir(parents=True, exist_ok=True)
FIG5 = ROOT / "paper_figure_outputs"
CI = ROOT / "outputs" / "interpretation" / "condition_influence_and_seed_stability_K46"

M = {28: "M0", 35: "M1", 16: "M2", 39: "M3", 45: "M4"}      # seed-0 id -> paper name (see the crosswalk)
CONDS = [("intact", "Intact narrative"), ("word", "Word-scrambled"), ("rest", "Rest")]
COLOR = {"intact": "#7B3294", "word": "#2CA25F", "rest": "#222222"}
NETS = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode"]
NET_SHORT = ["Vis", "SomMot", "DorsAttn", "VentAttn", "Limbic", "FP", "DMN"]
MODES = [("all_nodes", "All 700 nodes", "#2a78d6"), ("center_inside_atlas", "Center inside atlas", "#eb6834"),
         ("confidence_ge_0.50", "Assignment confidence \u2265 0.5", "#1baf7a")]


# Drawn at the journal's print width (Network Neuroscience reproduces figures ~5 in across), so the point sizes
# are the printed sizes.
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})


def style(ax, axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=axis, color="0.9", lw=0.5); ax.set_axisbelow(True)
    ax.tick_params(labelsize=6.5)


# ---------------- atlas audit
prof = pd.read_csv(FIG5 / "Figure5_Yeo7_composition_sensitivity_profiles.csv")
modes_present = [m for m in MODES if m[0] in set(prof["mode"])]
assert len(modes_present) == len(MODES), sorted(prof["mode"].unique())
aud = pd.read_csv(FIG5 / "Figure5_HTFA_Yeo_mapping_audit_summary.csv").set_index("metric")["value"]
# horizontal bars, so the full Yeo-7 network names can be written out once, unrotated (Lucy, 12 Nov)
fig, axes = plt.subplots(1, 5, figsize=(5.0, 2.9), sharey=True, sharex=True, constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02)
w = 0.8 / len(modes_present)
for ax, (aid, name) in zip(axes, M.items()):
    for k, (mode, lab, col) in enumerate(modes_present):
        x = prof[(prof.archetype == aid) & (prof["mode"] == mode)].set_index("network").reindex(NETS)
        ax.barh(np.arange(7) + (k - (len(modes_present) - 1) / 2) * w, x.mass_fraction, height=w * 0.92, color=col,
                label=lab)
    ax.set_yticks(range(7)); ax.set_yticklabels(NETS, fontsize=6.5)
    ax.set_ylim(6.6, -0.6)
    ax.set_title(name, fontsize=7.5, fontweight="bold", pad=2); style(ax, axis="x")
    ax.tick_params(axis="y", length=0)
axes[0].set_xlim(0, None)
axes[0].xaxis.set_major_locator(plt.MaxNLocator(2))
axes[2].set_xlabel("Share of |weight| (Yeo-7)", fontsize=7, labelpad=1)
fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, ncol=3, loc="outside lower center",
           handlelength=1.0, handleheight=0.8, columnspacing=1.0, handletextpad=0.4)
fig.suptitle("Atlas mapping audit", fontsize=8.5, fontweight="bold")
# the node-level audit numbers go in the caption; print them so the caption can be checked against the data
print(f"atlas audit: {int(aud['n_centers_inside_cortical_atlas'])}/{int(aud['n_nodes_total'])} centres inside atlas; "
      f"centre vs integrated agreement {aud['center_vs_integrated_assignment_agreement_pct_inside']:.1f}%; "
      f"{int(aud['n_assignment_confidence_below_0.50'])} nodes with confidence < 0.5")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Supplement_M0M4_atlas_mapping_audit.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

# ---------------- condition contributions
obs = pd.read_csv(CI / "selected_archetype_deployment_by_observation.csv")
obs = obs[obs.archetype.isin(M)]
loco = pd.read_csv(CI / "leave_one_condition_out_archetype_comparison.csv")
loco = loco[loco.reference_archetype.isin(M)]
fig, (axA, axB) = plt.subplots(1, 2, figsize=(4.52, 2.4), gridspec_kw=dict(width_ratios=[1.5, 1]), constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.05)
for i, aid in enumerate(M):
    for k, (c, lab) in enumerate(CONDS):
        v = obs[(obs.archetype == aid) & (obs.condition == c)].rank1_frobenius.to_numpy()
        pos = i + (k - 1) * 0.26
        bp = axA.boxplot(v, positions=[pos], widths=0.22, patch_artist=True, showfliers=False,
                         medianprops=dict(color="white", lw=1.0), whiskerprops=dict(color=COLOR[c], lw=0.7),
                         capprops=dict(color=COLOR[c], lw=0.7), boxprops=dict(lw=0.7))
        bp["boxes"][0].set(facecolor=COLOR[c], edgecolor=COLOR[c], alpha=0.85)
        axA.scatter(np.full(len(v), pos) + np.linspace(-0.07, 0.07, len(v)), v, s=2, color=COLOR[c], alpha=0.35, lw=0)
axA.set_xticks(range(len(M))); axA.set_xticklabels(M.values(), fontsize=7)
axA.set_ylabel("Rank-1 contribution per\nscan (Frobenius norm)", fontsize=7, labelpad=1)
fig.legend([plt.Rectangle((0, 0), 1, 1, color=COLOR[c]) for c, _ in CONDS] + [plt.Line2D([], [], color="#E69F00", ls="--", lw=1.0)],
           [l for _, l in CONDS] + ["Recovery threshold (r = 0.9)"], frameon=False,
           fontsize=6.5, ncol=4, loc="outside lower center", handlelength=1.0, handleheight=0.8, columnspacing=1.2,
           handletextpad=0.4, title="Condition (A) / condition left out (B)", title_fontsize=6.5)
axA.set_title("Deployment by condition", fontsize=7.5, loc="left", fontweight="bold", pad=2); style(axA)
axA.tick_params(axis="x", labelsize=7)
drops = [("drop_intact", "intact"), ("drop_word", "word"), ("drop_rest", "rest")]
for k, (d, c) in enumerate(drops):
    x = loco[loco.comparison == d].set_index("reference_archetype").reindex(list(M))
    axB.bar(np.arange(len(M)) + (k - 1) * 0.26, x.mean_spatial_map_correlation, width=0.24, color=COLOR[c])
axB.axhline(0.9, color="#E69F00", ls="--", lw=1.0)   # r = 0.90 recovery threshold, as in Fig 6D
# the same table and metric as Figure 6D: Pearson r of the group-mean spatial maps, Hungarian-matched on the
# generators, averaged over the two conditions present in both fits
axB.set_xticks(range(len(M))); axB.set_xticklabels(M.values(), fontsize=7)
axB.set_ylim(0, 1.02); axB.set_ylabel("Spatial-map r with\nfull-data motif", fontsize=7, labelpad=1)
axB.set_title("Leave-one-condition-out stability", fontsize=7.5, loc="left", fontweight="bold", pad=2); style(axB)
axB.tick_params(axis="x", labelsize=7)
for ax, l, dx in zip((axA, axB), "AB", (-0.2, -0.32)):
    ax.text(dx, 1.02, l, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Supplement_M0M4_condition_contributions.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
print("wrote", OUT)
