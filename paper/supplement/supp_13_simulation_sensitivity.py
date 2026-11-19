"""Supplementary Figure S13: simulation sensitivity of direct motif recovery.

Reads the one-at-a-time sweep written by
``analysis/07_simulations/s13_sensitivity_sweep.py`` and draws a 2 x 4 figure
sized for a ~5 inch print width:

* top row: observation noise, between-subject heterogeneity, fitted K,
  spatial motif overlap;
* bottom row: the four components of dynamical mixing (mean dwell, transition
  width, secondary-motif weight, AR(1) rho of the secondary background).

Visual language follows Figure 3A: MS-AA blue / PCA orange, mean with a
bootstrap 95% CI (5000 resamples), methods offset side by side. Because the
settings are ordered, thin lines connect the means (as in Figure 3C). Instead
of 20 replicate dots, a thin faint bar spans the replicate min to max. All
panels share a 0.4-1.0 y-axis. A grey band marks the shared Figure 3 baseline
setting. Numbers above each setting count replicates (of 20) in which MS-AA
recovery exceeds PCA recovery. The fitted-K panel uses coverage-adjusted
recovery (unmatched true motifs scored 0); short dashed segments mark its
ceiling K/6 for K < 6.

Usage (from the paper run directory):
    python paper/supplement/supp_13_simulation_sensitivity.py \
        --results simulation_pca_msaa_extended_outputs/s13_sensitivity_results.csv \
        --outdir outputs/paper_supplement
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METHOD_COLORS = {"MS-AA": "#eb6834", "PCA": "#2a78d6"}  # same as Fig 3 (PCA blue, MS-AA orange)
N_TRUE = 6
# Smallest win count out of 20 with two-sided exact sign-test p < 0.05 (p = 0.041).
SIGN_TEST_MIN = 15
DEFAULTS = dict(noise_sd=0.25, heterogeneity_scale=1.0, K_fit=6, distinctiveness_scale=1.0,
                mean_dwell=12, transition_width=4, secondary_weight=0.25, ar1_rho=0.85)

# (param, title, x label, level order easiest -> hardest)
PANELS = [
    ("noise_sd", "Noise", "Noise SD", [0.05, 0.10, 0.25, 0.40, 0.60]),
    ("heterogeneity_scale", "Heterogeneity", "Heterogeneity scale", [0.0, 0.5, 1.0, 1.5, 2.0]),
    ("K_fit", "Fitted K", "Fitted K (true K = 6)", [4, 5, 6, 7, 8]),
    ("distinctiveness_scale", "Spatial motif\noverlap", "Mean spatial-motif cosine", [1.15, 1.0, 0.85, 0.70, 0.55]),
    ("mean_dwell", "Dwell time", "Mean dwell (TRs)", [16, 14, 12, 10, 8]),
    ("transition_width", "Transition width", "Transition width (TRs)", [2, 3, 4, 6, 8]),
    ("secondary_weight", "Secondary weight", "Secondary-motif weight", [0.10, 0.175, 0.25, 0.35, 0.45]),
    ("ar1_rho", "Persistence", "AR(1) \u03c1", [0.75, 0.80, 0.85, 0.875, 0.90]),
]

STYLE = {
    "font.size": 6.5,
    "axes.titlesize": 7,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def mean_ci(x, n_boot=5000, seed=42):
    """Mean and percentile-bootstrap 95% CI (same procedure as Figure 3)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    rng = np.random.default_rng(seed)
    means = rng.choice(x, size=(n_boot, len(x)), replace=True).mean(axis=1)
    return float(x.mean()), np.quantile(means, [0.025, 0.975])


def fmt(v):
    # drop the leading zero so five ticks fit in a ~1 inch panel
    return f"{v:g}".replace("0.", ".", 1) if 0 < abs(v) < 1 else f"{v:g}"


def panel_table(res: pd.DataFrame, param: str, levels: list) -> pd.DataFrame:
    """Rows for one panel (baseline rows reused at the default level)."""
    use_cov = param == "K_fit"
    aa_col = "aa_motif_recovery_coverage" if use_cov else "aa_motif_recovery"
    pca_col = "pca_motif_recovery_coverage" if use_cov else "pca_motif_recovery"
    rows = []
    for i, lv in enumerate(levels):
        if np.isclose(lv, DEFAULTS[param]):
            cur = res[res["param"] == "baseline"]
        else:
            cur = res[(res["param"] == param) & np.isclose(res["level"].astype(float), lv)]
        cur = cur.sort_values("rep")
        assert len(cur) == 20, (param, lv, len(cur))
        aa, pca = cur[aa_col].to_numpy(float), cur[pca_col].to_numpy(float)
        seed_base = 1000 * (PANELS.index(next(p for p in PANELS if p[0] == param)) + 1) + 10 * i
        for method, vals, off in [("MS-AA", aa, 0), ("PCA", pca, 1)]:
            m, ci = mean_ci(vals, seed=seed_base + off)
            rows.append(dict(param=param, level=lv, x=i, method=method, mean=m, ci_lo=ci[0], ci_hi=ci[1],
                             min=vals.min(), max=vals.max(), n=len(vals),
                             n_msaa_gt_pca=int((aa > pca).sum()),
                             mean_diff=float((aa - pca).mean()),
                             mean_pairwise_cosine=float(cur["mean_pairwise_cosine"].iloc[0]),
                             metric=aa_col.replace("aa_", "")))
    return pd.DataFrame(rows)


def draw_panel(ax, tab, param, title, xlabel, letter, show_ylabel):
    offset = 0.14
    base_x = int(tab.loc[np.isclose(tab["level"], DEFAULTS[param]), "x"].iloc[0])
    ax.axvspan(base_x - 0.42, base_x + 0.42, color="0.92", zorder=0, lw=0)

    if param == "K_fit":
        for x, lv in enumerate([4, 5]):
            ax.plot([x - 0.36, x + 0.36], [lv / N_TRUE] * 2, ls=(0, (2, 1.5)), color="0.45",
                    lw=0.7, zorder=1)

    for method, sign in [("MS-AA", -1), ("PCA", 1)]:
        cur = tab[tab["method"] == method].sort_values("x")
        xs = cur["x"].to_numpy() + sign * offset
        c = METHOD_COLORS[method]
        ax.vlines(xs, cur["min"], cur["max"], color=c, alpha=0.35, lw=0.6, zorder=2)
        ax.plot(xs, cur["mean"], color=c, lw=0.8, zorder=3)
        yerr = np.vstack([cur["mean"] - cur["ci_lo"], cur["ci_hi"] - cur["mean"]])
        ax.errorbar(xs, cur["mean"], yerr=yerr, fmt="o", color=c, ms=2.8, elinewidth=1.1,
                    capsize=1.6, capthick=0.8, zorder=4, label=method)

    wins = tab[tab["method"] == "MS-AA"].sort_values("x")
    for x, n in zip(wins["x"], wins["n_msaa_gt_pca"]):
        ax.text(x, 0.975, f"{n}", ha="center", va="top", fontsize=5.5,
                color="0.15" if n >= SIGN_TEST_MIN else "#b2182b",
                fontweight="normal" if n >= SIGN_TEST_MIN else "bold")

    if param == "distinctiveness_scale":
        labels = [f"{v:.2f}".lstrip("0") for v in wins["mean_pairwise_cosine"]]
    else:
        labels = [fmt(v) for v in wins["level"]]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(0.4, 1.0)
    ax.set_yticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel(xlabel, labelpad=1.5)
    ax.set_title(title, fontweight="bold", pad=3)
    ax.text(-0.10, 1.03, letter, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    if show_ylabel:
        ax.set_ylabel("Direct spatial\nrecovery (r)", labelpad=2)
    else:
        ax.tick_params(labelleft=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="simulation_pca_msaa_extended_outputs/s13_sensitivity_results.csv")
    ap.add_argument("--outdir", default="outputs/paper_supplement")
    ap.add_argument("--width", type=float, default=5.0)
    ap.add_argument("--height", type=float, default=3.4)
    a = ap.parse_args()

    res = pd.read_csv(a.results)
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tables = [panel_table(res, p, lv) for p, _, _, lv in PANELS]
    summary = pd.concat(tables, ignore_index=True)

    plt.rcParams.update(STYLE)
    fig, axes = plt.subplots(2, 4, figsize=(a.width, a.height), sharey=True)
    fig.subplots_adjust(left=0.105, right=0.995, top=0.88, bottom=0.115, wspace=0.12, hspace=0.80)
    for i, (tab, (p, title, xlabel, _)) in enumerate(zip(tables, PANELS)):
        ax = axes.flat[i]
        draw_panel(ax, tab, p, title, xlabel, "ABCDEFGH"[i], show_ylabel=(i % 4 == 0))

    handles, labels = axes.flat[0].get_legend_handles_labels()
    # say in the legend that the grey band is the shared baseline setting (ben, 14 Nov)
    from matplotlib.patches import Patch
    handles.append(Patch(facecolor="0.92", edgecolor="none"))
    labels.append("Baseline setting")
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.55, 0.985), handletextpad=0.3, columnspacing=1.5)
    # header over the bottom row
    top_b = axes[1, 0].get_position().y1
    fig.text(0.55, top_b + 0.075, "Dynamical mixing components (others held at baseline)",
             ha="center", va="bottom", fontsize=7, style="italic")

    stem = outdir / "S13_simulation_sensitivity"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.02)
    summary.to_csv(outdir / "S13_simulation_sensitivity_source.csv", index=False)
    print("wrote", stem.with_suffix(".pdf"), stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
