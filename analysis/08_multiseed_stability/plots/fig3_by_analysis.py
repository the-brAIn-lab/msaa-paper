"""Figure 3 split by analysis: one 2 x 2 figure per analysis type (spatial or temporal).

Left column = paper fits (seed 0); right column = mean of the 10 fits per K (seeds 0-9).
Top row: across-condition timepoint decoding accuracy (full model) vs K, paper K grid.
Bottom row: HELD-OUT top-m accuracy, intact only, every K (rank on half A, score on half B; 20 outer splits).
Same y-limits in both columns so seed 0 and the 10-fit mean can be compared directly.
Optional peak-K bootstrap CI strips (peakk_ci.py output, raw W=1 versions).

Usage:
  python fig3_by_analysis.py ANALYSIS [--summary0 full_summary.csv] [--perk per_K.csv] [--heldout0 DIR] [--heldoutA DIR]
         [--out OUTDIR] [--peakci0 seed0_ci.csv] [--peakciA all_ci.csv]
  Defaults (paths.py): seed-0 decoding summary, seeds/results/per_K.csv, heldout-everyk-ANALYSIS-{seed0,all}-results/,
  OUTDIR outputs/multiseed_stability/fig3_by_analysis; the peak-K CI strips use the peakk_ci.py W=1 outputs if present.
  ANALYSIS = spatial | temporal
"""
import argparse, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
try:
    from scipy import stats
    def t975(df): return stats.t.ppf(0.975, df)
except ImportError:
    _T = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228}
    def t975(df): return df.map(lambda v: _T.get(int(v), 1.96))

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("analysis", choices=["spatial", "temporal"])
for o in ("summary0", "perk", "heldout0", "heldoutA", "out", "peakci0", "peakciA"): ap.add_argument(f"--{o}", default=None)
a = ap.parse_args(); AT = a.analysis
SUMM0 = a.summary0 or paths.decoding_summary(AT); PERK = a.perk or paths.SEEDS_RESULTS / "per_K.csv"
HO0 = a.heldout0 or paths.heldout_everyk(AT, "seed0"); HOA = a.heldoutA or paths.heldout_everyk(AT, "all")
OUT = Path(a.out or paths.OUT / "fig3_by_analysis")
CI = {}
for given, which, key in ((a.peakci0, "seed0", "seed0"), (a.peakciA, "all", "seedavg")):
    f = Path(given) if given else paths.OUT / "peakk_ci" / f"{AT}_{which}_W1.csv"
    if f.exists():
        c_ = pd.read_csv(f); CI[key] = c_[c_.condition == "intact"].set_index("top_m")

COND = {"intact": "#7b2d8e", "word": "#2ca02c", "rest": "#222222"}
TOPM = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c", 5: "#d62728", 10: "#8c564b"}
SHADE_M = (5, 10)
REP_M = 10
XMAX = {"spatial": 700, "temporal": 300}[AT]
GRID = {"spatial": [5, 7, 10, 14, 21, 25, 35, 40, 42, 44, 45, 46, 48, 50, 52, 53, 54, 55, 56, 58, 60, 63, 65, 70,
                    75, 88, 100, 105, 126, 140, 175, 189, 200, 210, 252, 300, 400, 500, 600, 700],
        "temporal": [2, 3, 5, 6, 8, 9, 10, 11, 14, 15, 17, 18, 20, 23, 25, 27, 30, 35, 38, 40, 45, 50, 54, 55, 60,
                     65, 75, 81, 90, 100, 108, 129, 171, 200, 214, 257, 300]}[AT]

# --- decoding (top row)
d0 = pd.read_csv(SUMM0)
d0["lo"], d0["hi"] = d0["mean"] - 1.96 * d0["err"], d0["mean"] + 1.96 * d0["err"]
pk = pd.read_csv(PERK)
pk = pk[(pk.fit_scope == "across") & (pk.metric == "full_mean") & (pk.analysis_type == AT) & pk.K.isin(GRID)].copy()
h = pk["sd"] / pk["n_seeds"] ** 0.5 * t975(pk["n_seeds"] - 1)
pk["lo"], pk["hi"] = pk["mean"] - h, pk["mean"] + h
DEC = {"seed0": d0, "seedavg": pk}


def load_heldout(dirname):
    d = Path(dirname)
    p = pd.read_csv(d / "heldout_per_K.csv"); p = p[p.condition == "intact"]
    R = pd.read_csv(d / "heldout_topm_runs_all_seeds.csv.gz")
    R = R[R.condition == "intact"].groupby(["fit_seed", "K", "top_m", "outer_rep"])["mean"].mean()
    sem = R.groupby(["fit_seed", "K", "top_m"]).sem().groupby(["K", "top_m"]).mean().rename("e").reset_index()
    return p.merge(sem, on=["K", "top_m"])


HO = {"seed0": load_heldout(HO0), "seedavg": load_heldout(HOA)}
COLS = ["seed0", "seedavg"]
TITLE = {"seed0": "Paper fits (seed 0)", "seedavg": "Mean of 10 fits (seeds 0–9)"}

fig, axes = plt.subplots(2, 2, figsize=(10, 8), gridspec_kw=dict(hspace=0.32, wspace=0.2))
peaks = {}
for col, mode in enumerate(COLS):
    ax = axes[0, col]; d = DEC[mode]
    for c, color in COND.items():
        x = d[d.condition == c].sort_values("K")
        ax.fill_between(x.K, x.lo, x.hi, color=color, alpha=0.2, lw=0)
        ax.plot(x.K, x["mean"], "-o", color=color, ms=2.2, lw=1.6)
    ax.set_title(f"{AT.capitalize()} AA: {TITLE[mode]}", fontsize=10.5)
    ax.set_ylabel("Across-condition decoding accuracy", fontsize=8.5)
    ax.text(-0.08, 1.06, "AB"[col], transform=ax.transAxes, fontsize=12)

    ax = axes[1, col]; ho = HO[mode]
    for m, color in TOPM.items():
        x = ho[(ho.top_m == m) & (ho.K > m)].sort_values("K")
        ax.fill_between(x.K, x["mean"] - x.e, x["mean"] + x.e, color=color, alpha=0.25, lw=0)
        ax.plot(x.K, x["mean"], "-", color=color, lw=1.1, label=f"top-{m}")
    r = ho[(ho.top_m == REP_M) & (ho.K > REP_M)]
    peaks[mode] = int(r.loc[r["mean"].idxmax(), "K"])
    ax.axhline(1 / 300, color="gray", ls="--", lw=0.8)
    ax.set_ylabel("Held-out top-m decoding accuracy\nfor intact condition", fontsize=8.5)
    ax.legend(title="Top-m", frameon=False, fontsize=7.5, title_fontsize=8.5, loc="upper right", ncol=2)
    ax.text(-0.08, 1.06, "CD"[col], transform=ax.transAxes, fontsize=12)

# shared y-limits per row
for row in (0, 1):
    top = max(a.get_ylim()[1] for a in axes[row])
    if row == 1 and CI:
        top *= 1.45
    elif row == 1:
        top *= 1.08
    for a in axes[row]:
        a.set_ylim(0, top)

# peak-K CI strips (after y-limits are fixed)
for col, mode in enumerate(COLS):
    ax = axes[1, col]
    if mode in CI:
        ci = CI[mode]
        for j, (m, color) in enumerate(TOPM.items()):
            if m not in ci.index: continue
            lo_, hi_, pk_ = ci.loc[m, "ci_lo"], ci.loc[m, "ci_hi"], ci.loc[m, "peak_K"]
            if m in SHADE_M:
                ax.axvspan(lo_ - 0.5, hi_ + 0.5, color=color, alpha=0.13, lw=0, zorder=0)
            y = 0.97 - 0.035 * j
            tr = ax.get_xaxis_transform()
            ax.plot([lo_, hi_], [y, y], "-", color=color, lw=3, solid_capstyle="butt", transform=tr)
            ax.plot([pk_], [y], "o", color="white", mec=color, ms=3.5, mew=1.2, transform=tr)
            ax.text(hi_ + XMAX * 0.012, y, f"{lo_}–{hi_}" if lo_ != hi_ else f"{lo_}",
                    transform=tr, va="center", fontsize=6.5, color=color)
    else:
        ax.axvline(peaks[mode], color="0.6", ls=":", lw=0.9)
        ax.text(peaks[mode], 0.99, f" K={peaks[mode]}", transform=ax.get_xaxis_transform(), ha="left", va="top", fontsize=9)

axes[0, 1].legend([Patch(color=v) for v in COND.values()], [k.capitalize() for k in COND], title="Condition",
                  frameon=False, fontsize=8, title_fontsize=9, loc="upper left", bbox_to_anchor=(0.03, 1.0))
for a in axes.flat:
    a.tick_params(labelsize=8); a.spines[["top", "right"]].set_visible(False)
    a.set_xlabel("Number of archetypes (K)", fontsize=9); a.set_xlim(-XMAX * 0.015, XMAX * 1.03)
fig.suptitle(f"Figure 3 ({AT}), held-out top-m at every K. Left: paper fits (seed 0); right: mean of 10 fits per K.\n"
             "Decoding bands: 95% CI (right: of the 10-fit mean)."
             "Held-out (C/D) bands: SEM over 20 outer splits (right: mean over fits). Dashed: chance 1/300."
             + ("\nStrips (top): bootstrap 95% CI of the K where each top-m curve peaks (o = peak); shaded: top-5/top-10 CIs."
                if CI else ""),
             fontsize=8.5, y=0.995)
OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure3_{AT}_seed0_vs_seedavg.{ext}", dpi=300, bbox_inches="tight")
print(f"wrote {OUT}/Figure3_{AT}_seed0_vs_seedavg; intact held-out top-{REP_M} peaks:", peaks)
