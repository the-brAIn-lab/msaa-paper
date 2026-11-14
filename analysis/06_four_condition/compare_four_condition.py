"""Compare the four-condition tables of two runs (e.g. August 2026 vs the current-pipeline rerun).

For each condition x method (raw HTFA, PCA, spatial MS-AA, temporal MS-AA) at the shared K values:
max and mean |accuracy difference|, Pearson r between the old and new accuracy curves, and the
largest difference in SEM units (|diff| / sqrt(err_old^2 + err_new^2)). Also checks the
Intact > Paragraph-scrambled > Word-scrambled > Rest ordering of raw HTFA accuracy and, per
method, at every shared ratio. Writes comparison.csv, raw_htfa_old_vs_new.csv, ordering.csv and an
old-vs-new overlay figure (PNG + PDF).

Usage: python compare_four_condition.py --old DIR --new DIR --out DIR
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser()
ap.add_argument("--old", required=True)
ap.add_argument("--new", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
OLD, NEW, OUT = Path(a.old), Path(a.new), Path(a.out)
OUT.mkdir(parents=True, exist_ok=True)

CONDS = [("intact", "Intact narrative"), ("paragraph", "Paragraph-scrambled"),
         ("word", "Word-scrambled"), ("rest", "Rest")]
METHODS = [("PCA reconstruction", "#2a78d6"), ("Spatial MS-AA", "#eb6834"), ("Temporal MS-AA", "#1baf7a")]


def load(d):
    raw = pd.read_csv(d / "raw_htfa_summary.csv")
    pca = pd.read_csv(d / "pca_summary_with_retention.csv").assign(method="PCA reconstruction")
    aa = pd.read_csv(d / "aa_full_summary_with_retention.csv")
    aa["method"] = aa.analysis_type.map({"spatial": "Spatial MS-AA", "temporal": "Temporal MS-AA"})
    raw = raw.assign(method="Raw HTFA")
    cols = ["condition", "method", "K", "target_ratio", "actual_ratio", "mean", "err"]
    return pd.concat([raw[cols], pca[cols], aa[cols]], ignore_index=True)


old, new = load(OLD), load(NEW)
m = old.merge(new, on=["condition", "method", "K"], suffixes=("_old", "_new"))
m["diff"] = m.mean_new - m.mean_old
m["z"] = m["diff"].abs() / np.sqrt(m.err_old ** 2 + m.err_new ** 2)
rows = []
for (c, meth), g in m.groupby(["condition", "method"], sort=False):
    rows.append(dict(condition=c, method=meth, n_shared=len(g),
                     n_old=int(((old.condition == c) & (old.method == meth)).sum()),
                     n_new=int(((new.condition == c) & (new.method == meth)).sum()),
                     max_abs_diff=g["diff"].abs().max(), mean_abs_diff=g["diff"].abs().mean(),
                     mean_signed_diff=g["diff"].mean(),
                     r=np.corrcoef(g.mean_old, g.mean_new)[0, 1] if len(g) > 2 else np.nan,
                     max_diff_in_sem=g.z.max()))
comp = pd.DataFrame(rows)
comp.to_csv(OUT / "comparison.csv", index=False)
m.to_csv(OUT / "comparison_rows.csv", index=False)

rawc = m[m.method == "Raw HTFA"][["condition", "mean_old", "mean_new", "err_old", "err_new", "diff"]]
rawc.to_csv(OUT / "raw_htfa_old_vs_new.csv", index=False)

# ordering Intact > Paragraph > Word > Rest, by run, per method and target ratio
order = [c for c, _ in CONDS]
ordrows = []
for run, df in (("old", old), ("new", new)):
    df = df.assign(ratio=df.target_ratio)
    for meth, g in df.groupby("method"):
        for ratio, h in g.groupby("ratio"):
            v = h.set_index("condition")["mean"]
            if set(order) <= set(v.index):
                ordrows.append(dict(run=run, method=meth, ratio=ratio,
                                    holds=bool(all(v[order[i]] > v[order[i + 1]] for i in range(3))),
                                    para_between=bool(v["intact"] > v["paragraph"] > v["word"])))
ordt = pd.DataFrame(ordrows)
ordt.to_csv(OUT / "ordering.csv", index=False)

# overlay figure: old dashed, new solid
fig, axes = plt.subplots(1, 4, figsize=(16, 4.4), sharex=True, sharey=True)
for ax, (c, title) in zip(axes, CONDS):
    for name, col in METHODS:
        for run, df, ls in (("old", old, "--"), ("new", new, "-")):
            x = df[(df.condition == c) & (df.method == name)].sort_values("actual_ratio")
            ax.plot(x.actual_ratio, x["mean"], ls, color=col, lw=1.8 if run == "new" else 1.3,
                    alpha=1 if run == "new" else 0.8)
    for run, df, ls in (("old", old, (0, (1, 3))), ("new", new, ":")):
        r = df[(df.condition == c) & (df.method == "Raw HTFA")]["mean"]
        if len(r):
            ax.axhline(float(r.iloc[0]), color="#222222", ls=ls, lw=1.2)
    ax.set_title(title, fontsize=11, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", lw=0.6); ax.set_axisbelow(True)
    ax.set_xlabel("Normalized component ratio", fontsize=9)
axes[0].set_ylabel("Timepoint-identification accuracy", fontsize=9)
handles = [Line2D([], [], color=c, lw=2.5, label=n) for n, c in METHODS] + [
    Line2D([], [], color="0.3", ls="--", lw=1.3, label="August 2026 run"),
    Line2D([], [], color="0.3", ls="-", lw=1.8, label="Current pipeline")]
fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.06))
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"S11_old_vs_new_overlay.{ext}", dpi=200, bbox_inches="tight")

pd.set_option("display.width", 200)
print(comp.round(4).to_string(index=False))
print(rawc.round(4).to_string(index=False))
print(ordt.groupby(["run", "method"])[["holds", "para_between"]].agg(["sum", "count"]).to_string())
