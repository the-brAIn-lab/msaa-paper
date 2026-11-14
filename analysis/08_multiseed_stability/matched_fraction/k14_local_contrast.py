"""Is the K=14 matched-fraction retention bump driven by the subset numerator, the full-model denominator,
or both? Uses only the aggregate outputs of temporal_matched_fraction.py (no new decoding).

Writes seed_level_table.csv (seed, condition, K, requested fraction, selected m, actual fraction, raw top-m
and full-model accuracy, above-chance retention), delta14_contrast.csv (paired local contrast
Delta14 = R(14) - [R(12) + R(16)] / 2 over the 10 fit seeds, bootstrap 95% CI, seeds with Delta14 > 0;
at 50% the fraction is exactly matched: 6/12 = 7/14 = 8/16) and three plots.
    python k14_local_contrast.py AGG_DIR OUT_DIR"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

agg, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
CH = 1 / 300
s = pd.read_csv(agg / "matched_fraction_seed_summary.csv")
t = pd.DataFrame({
    "fit_seed": s.fit_seed.astype(int), "condition": s.condition, "K": s.K, "requested_fraction": s.fraction_requested,
    "selected_m": s.top_m, "actual_fraction": s.actual_fraction,
    "topm_raw_accuracy": s.mean_accuracy, "full_model_raw_accuracy": s.mean_full_accuracy,
    "retention_above_chance": s.mean_retention,
    "retention_of_means": (s.mean_accuracy - CH) / (s.mean_full_accuracy - CH),
    "n_outer_splits": s.n_outer})
t = t.sort_values(["condition", "requested_fraction", "K", "fit_seed"])
t.to_csv(out / "seed_level_table.csv", index=False)

rng = np.random.default_rng(20261111)
rows = []
for (c, f), g in t.groupby(["condition", "requested_fraction"]):
    if f == 1.0:
        continue
    p = g.pivot(index="fit_seed", columns="K")
    for metric in ["retention_above_chance", "topm_raw_accuracy", "full_model_raw_accuracy"]:
        x = p[metric]
        d = (x[14] - (x[12] + x[16]) / 2).to_numpy()
        boots = rng.choice(d, (10000, len(d)), replace=True).mean(1)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append(dict(condition=c, requested_fraction=f, metric=metric,
                         m_at_12_14_16="/".join(str(int(p["selected_m"][k].iloc[0])) for k in (12, 14, 16)),
                         n_seeds=len(d), mean_delta14=d.mean(), ci95_lo=lo, ci95_hi=hi, n_positive=int((d > 0).sum()),
                         per_seed=" ".join(f"{v:+.4f}" for v in d)))
con = pd.DataFrame(rows)
con.to_csv(out / "delta14_contrast.csv", index=False)

# summaries for plots: mean +- SD over seeds
m = t.groupby(["condition", "requested_fraction", "K"]).agg(
    acc=("topm_raw_accuracy", "mean"), acc_sd=("topm_raw_accuracy", "std"),
    full=("full_model_raw_accuracy", "mean"), full_sd=("full_model_raw_accuracy", "std"),
    n=("fit_seed", "nunique")).reset_index()
assert (m.n == 10).all()
conds = ["intact", "word"]

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=False)
for a, c in zip(ax, conds):
    for f in [0.25, 0.5, 0.7, 1.0]:
        g = m[(m.condition == c) & np.isclose(m.requested_fraction, f)]
        lab = "full model (100%)" if f == 1.0 else f"{f:.0%} retained"
        a.plot(g.K, g.acc, marker="o", lw=2.4 if f == 0.5 else 1.4, label=lab)
        a.fill_between(g.K, g.acc - g.acc_sd, g.acc + g.acc_sd, alpha=.15)
    a.axhline(CH, ls=":", c="grey", lw=1); a.axvspan(11.5, 16.5, color="k", alpha=.05)
    a.set_title(f"{c}: raw held-out accuracy of matched subset"); a.set_xlabel("Temporal K"); a.set_ylabel("Accuracy (chance = 1/300, dotted)")
    a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(out / "raw_subset_accuracy_vs_K.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for a, c in zip(ax, conds):
    g = m[(m.condition == c) & np.isclose(m.requested_fraction, 1.0)]
    a.plot(g.K, g.full, marker="o", c="k"); a.fill_between(g.K, g.full - g.full_sd, g.full + g.full_sd, color="k", alpha=.12)
    for k, sub in t[(t.condition == c) & (t.requested_fraction == 1.0)].groupby("K"):
        a.scatter([k] * len(sub), sub.full_model_raw_accuracy, s=8, c="grey", alpha=.6)
    a.axvspan(11.5, 16.5, color="k", alpha=.05)
    a.set_title(f"{c}: full-model raw held-out accuracy (denominator)"); a.set_xlabel("Temporal K"); a.set_ylabel("Accuracy")
fig.tight_layout(); fig.savefig(out / "full_model_accuracy_vs_K.png", dpi=200); plt.close(fig)

# supplement figure: per-seed retention at K = 12, 14, 16. The bootstrap CI over fit seeds stays in the CSV only;
# seeds are initializations, not independent samples, so the figure shows the seed count with Delta14 > 0.
LABEL = {"intact": "Intact narrative", "word": "Word-scrambled"}
COLOR = {"intact": "#7B3294", "word": "#2CA25F"}
# drawn at the journal's print width (~5 in), so these point sizes are the printed sizes
plt.rcParams.update({"axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5})
fig, ax = plt.subplots(2, 3, figsize=(5.0, 3.9), constrained_layout=True)
fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.04, hspace=0.06)
for i, c in enumerate(conds):
    for j, f in enumerate([0.25, 0.5, 0.7]):
        a = ax[i, j]
        g = t[(t.condition == c) & np.isclose(t.requested_fraction, f) & t.K.isin([12, 14, 16])]
        p = g.pivot(index="fit_seed", columns="K", values="retention_above_chance")
        for k_, (sd, r) in enumerate(p.iterrows()):
            a.plot([12, 14, 16], r.values, marker="o", ms=1.6, lw=0.5, color="0.65",
                   label="individual fit seed" if k_ == 0 else None)
        a.plot([12, 14, 16], p.mean().values, c=COLOR[c], lw=1.4, marker="s", ms=3, label="mean of 10 seeds")
        r0 = con[(con.condition == c) & np.isclose(con.requested_fraction, f) & (con.metric == "retention_above_chance")].iloc[0]
        ms = g.groupby("K").selected_m.first()
        a.set_xticks([12, 14, 16], [f"{k}\n(m = {ms[k]})" for k in (12, 14, 16)], fontsize=6.5)
        a.tick_params(axis="y", labelsize=6.5)
        a.set_xlim(11.3, 16.7)
        a.set_title(f"{LABEL[c]}, {f:.0%} of K\n", fontsize=7, fontweight="bold", pad=1)
        a.text(0.5, 1.02, f"\u039414 = {r0.mean_delta14:+.3f}; {r0.n_positive}/10 seeds > 0", transform=a.transAxes,
               fontsize=6, color="0.3", ha="center", va="bottom")
        a.spines[["top", "right"]].set_visible(False)
        if j == 0:
            a.set_ylabel("Fraction of above-chance\nfull-model accuracy retained", fontsize=6.5, labelpad=1)
        if i == 1:
            a.set_xlabel("Temporal model order (K)", fontsize=7, labelpad=1)
# legend below the panels, clear of the data (ben, 14 Nov); grey = one fit seed, colour = condition mean
from matplotlib.lines import Line2D  # noqa: E402
fig.legend(handles=[Line2D([], [], color="0.65", lw=0.5, marker="o", ms=1.6, label="Individual fit seed"),
                    Line2D([], [], color=COLOR["intact"], lw=1.4, marker="s", ms=3, label="Mean, Intact narrative"),
                    Line2D([], [], color=COLOR["word"], lw=1.4, marker="s", ms=3, label="Mean, Word-scrambled")],
           fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.015), ncol=3, handlelength=1.8,
           columnspacing=1.2)   # placed just below the figure, with a gap under the x labels
fig.suptitle("Temporal MS-AA: local contrast at K = 12, 14, 16", fontsize=8.5, fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig(out / f"per_seed_K12_14_16.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

pd.set_option("display.width", 200)
print(con.drop(columns="per_seed").round(4).to_string(index=False))
print(m[m.K.isin([12, 14, 16])].round(4).to_string(index=False))
