#!/usr/bin/env python3
"""Cross-seed summary of the held-out top-m.  python hsummary.py RUN_DIR MODE FIRST LAST
Writes RUN_DIR/heldout-seeds[-smoke]-results/:
  heldout_topm_runs_all_seeds.csv.gz  every split of every seed, with fit_seed (held-out bookkeeping check)
  heldout_per_K.csv      condition x K x top_m: held-out accuracy (outer-split mean per seed) over seeds
                         (mean, sd, min, max, seed0); also mean full-model accuracy and retention
  selection_by_seed.csv  held-out selection frequency per seed, archetypes renamed to seed-0 IDs
                         (Hungarian match on group-mean spatial maps, map_r = match quality)
  K46_top5.csv           K=46 top-5: per seed-0 archetype and condition, frequency in every seed
  summary.txt            first line SUMMARY COMPLETE / INCOMPLETE
"""
import os
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

RUN, MODE, FIRST, LAST = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
SUB = "msaa_heldout_topm_outputs_spatial_across"
PAPER = [28, 35, 37, 31, 39]
OUT = RUN / ("heldout-seeds-results" if MODE == "full" else "heldout-seeds-smoke-results")
OUT.mkdir(exist_ok=True)
SEEDS = [0] + list(range(FIRST, LAST + 1))
KS = sorted(set(pd.read_csv(RUN / "configs" / "heldout_outer_seeds.csv").K.astype(int)))
if MODE == "smoke":
    KS = [5, 46]
if os.environ.get("HS_KS"):            # set by multiseed/run_local.py --ks
    KS = [k for k in KS if k in {int(x) for x in os.environ["HS_KS"].split(",")}]


def seed_root(s):
    return RUN if s == 0 else RUN / "seeds" / "full" / f"seed{s:02d}"


def slim(s, K):
    p = (RUN / "seeds" / "slim" / "seed00" if s == 0 else seed_root(s) / "slim") / f"spatialAA_across_acrossCond_K{K}.npz"
    return np.load(p)["Smean"].astype(float)


def corr(A, B):
    A = (A - A.mean(1, keepdims=True)) / A.std(1, keepdims=True)
    B = (B - B.mean(1, keepdims=True)) / B.std(1, keepdims=True)
    return A @ B.T / A.shape[1]


lines, ok = [], True
summ, stab, allruns = [], [], []
for s in SEEDS:
    d = seed_root(s) / SUB
    runs = pd.read_csv(d / "heldout_topm_crossfit_runs.csv", usecols=["condition", "K", "outer_rep", "direction"])
    runs = runs[runs.K.isin(KS)]
    full = pd.read_csv(d / "heldout_topm_crossfit_runs.csv")
    allruns.append(full[full.K.isin(KS)].assign(fit_seed=s))
    n = runs.groupby(["condition", "K", "outer_rep"]).direction.nunique().eq(2).sum()
    exp = len(KS) * 3 * 20
    lines.append(f"seed {s}: {n}/{exp} complete outer splits")
    ok &= n == exp
    a = pd.read_csv(d / "heldout_topm_summary.csv"); a = a[a.K.isin(KS)].assign(seed=s); summ.append(a)
    b = pd.read_csv(d / "heldout_archetype_selection_stability.csv"); b = b[b.K.isin(KS)].assign(seed=s); stab.append(b)
summ, stab = pd.concat(summ), pd.concat(stab)
# every held-out split of every seed in one table (fit_seed, outer_rep, direction, condition, K, top_m, ...)
pd.concat(allruns, ignore_index=True).to_csv(OUT / "heldout_topm_runs_all_seeds.csv.gz", index=False)

g = summ.groupby(["condition", "K", "top_m"])
per_K = g["mean"].agg(n_seeds="count", mean="mean", sd="std", min="min", max="max").reset_index()
per_K = per_K.merge(summ[summ.seed == 0][["condition", "K", "top_m", "mean"]].rename(columns={"mean": "seed0"}))
per_K = per_K.merge(g[["mean_full_same_eval", "mean_retention"]].mean().reset_index())
per_K.to_csv(OUT / "heldout_per_K.csv", index=False)

rows = []
for K in KS:
    S0 = slim(0, K)
    for s in SEEDS:
        R = corr(S0, slim(s, K)); i, j = linear_sum_assignment(-R)
        to0 = dict(zip(j, i)); r = dict(zip(j, R[i, j]))
        t = stab[(stab.seed == s) & (stab.K == K)]
        for _, x in t.iterrows():
            rows.append(dict(K=K, condition=x.condition, top_m=int(x.top_m), seed=s, archetype=int(x.archetype),
                             archetype0=int(to0[int(x.archetype)]), map_r=float(r[int(x.archetype)]),
                             selection_frequency=x.selection_frequency))
sel = pd.DataFrame(rows)
sel.to_csv(OUT / "selection_by_seed.csv", index=False)

k46 = sel[(sel.K == 46) & (sel.top_m == 5)]
if len(k46):
    t = k46.pivot_table(index=["archetype0", "condition"], columns="seed", values="selection_frequency", fill_value=0)
    t["mean"] = t.mean(axis=1); t["seeds_ge_0.5"] = (t.drop(columns="mean") >= 0.5).sum(axis=1)
    t = t.sort_values("mean", ascending=False)
    t.to_csv(OUT / "K46_top5.csv")
    lines += ["", "K=46 top-5 held-out selection frequency (seed-0 IDs; * = Fig 5 archetype), top 12 rows:",
              t.head(12).round(2).rename(index=lambda a: f"{a}*" if a in PAPER else a, level=0).to_string()]

for c in ["intact", "word", "rest"]:
    x = per_K[(per_K.condition == c) & (per_K.top_m.isin([5, 10]))]
    lines.append(f"\n{c}: held-out top-5/10 accuracy over seeds, median across K of seed SD = "
                 f"{x.sd.median():.4f}; at K=46: " + ", ".join(
                     f"top-{m} {r['mean']:.4f} +- {r['sd']:.4f} (seed0 {r.seed0:.4f})"
                     for m in (5, 10) for _, r in x[(x.K == 46) & (x.top_m == m)].iterrows()))
head = "SUMMARY COMPLETE" if ok else "SUMMARY INCOMPLETE"
(OUT / "summary.txt").write_text("\n".join([f"{head} ({MODE}, seeds {SEEDS[0]}-{SEEDS[-1]}, {len(KS)} K)"] + lines) + "\n")
print((OUT / "summary.txt").read_text(), flush=True)
with tarfile.open(OUT.with_suffix(".tar.gz"), "w:gz") as tf:
    tf.add(OUT, arcname=OUT.name)
sys.exit(0 if ok else 1)
