#!/usr/bin/env python3
"""Merge one phase and write its results bundle.   python hk_merge.py RUN_DIR smoke|full PHASE WORKERS
Per seed (cached in BASE/merged/<analysis>/seedNN once complete):
  heldout_topm_crossfit_runs.csv.gz      every outer split x direction x top-m (new parts + reused 30-K rows)
  heldout_topm_summary.csv               msaa_pipeline write_heldout_summaries (directions averaged per split)
  heldout_archetype_selection_stability.csv   same logic as write_heldout_summaries (from each ranking's top 10)
  heldout_per_archetype_mean.csv.gz      condition x K x archetype: single-archetype accuracy on ranking halves,
                                         mean over the 40 rankings (the full per-split ranks stay in BASE/parts)
Phase bundle RUN_DIR/heldout-everyk-<analysis>-<seed0|all>[-smoke]-results(.tar.gz):
  heldout_per_K.csv (condition x K x top_m over seeds: n_seeds mean sd min max seed0 + full + retention),
  per-seed summaries/stability, heldout_topm_runs_all_seeds.csv.gz (with fit_seed), and for spatial
  selection_by_seed.csv (seed-0 IDs, Hungarian on slim group-mean maps) + K46_top5.csv; summary.txt
  (first line SUMMARY COMPLETE / INCOMPLETE). sA/tA bundles include seed 0. Exit 1 if incomplete."""
import shutil
import sys
import tarfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

KIT = Path(__file__).resolve().parent
sys.path.insert(0, str(KIT))
from hk_common import CONDS, NREP, OLD_SUB, TOP_M, base_dir, old_complete_ks, phase, seed_root, split_dir  # noqa

RUN, MODE, PH, WORKERS = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
sys.path.insert(0, str(RUN / "src"))
from msaa_pipeline.heldout import write_heldout_summaries  # noqa: E402

BASE = base_dir(RUN, MODE)
AT, SEEDS, KS = phase(MODE, PH)
if PH in ("sA", "tA"):
    SEEDS = [0] + SEEDS
PAPER = [28, 35, 37, 31, 39]
KEY = ["analysis_type", "fit_scope", "condition", "K", "outer_rep", "direction"]
SORT = ["mean_accuracy", "rank_mean", "archetype"]


def stability_long(ranks):
    """Top-10 of each ranking, ordered exactly as write_heldout_summaries orders them."""
    out = []
    for keys, sub in ranks.groupby(KEY):
        top = sub.sort_values(SORT, ascending=[False, False, True])["archetype"].astype(int).tolist()[:max(TOP_M)]
        out += [dict(zip(KEY, keys), rank_position=i, archetype=a) for i, a in enumerate(top, 1)]
    return pd.DataFrame(out)


def per_arch(ranks):
    return ranks.groupby(["analysis_type", "condition", "K", "archetype"], as_index=False).agg(
        mean_accuracy=("mean_accuracy", "mean"), sd_accuracy=("mean_accuracy", "std"), n_rankings=("mean_accuracy", "size"))


def one_k(args):
    at, s, K = args
    runs, ranks, n = [], [], 0
    for c in CONDS:
        for r in range(NREP):
            d = split_dir(BASE, at, s, K, c, r)
            if (d / "runs.csv").exists():
                runs.append(pd.read_csv(d / "runs.csv")); ranks.append(pd.read_csv(d / "ranks.csv")); n += 1
    if not n:
        return K, 0, None, None, None
    runs, ranks = pd.concat(runs, ignore_index=True), pd.concat(ranks, ignore_index=True)
    return K, n, runs, stability_long(ranks), per_arch(ranks)


def merge_seed(s):
    out = BASE / "merged" / AT / f"seed{s:02d}"
    if (out / "COMPLETE").exists():
        print(f"seed {s}: merged earlier, reused", flush=True)
        return out, True, (out / "COMPLETE").read_text().strip()
    out.mkdir(parents=True, exist_ok=True)
    reuse = (old_complete_ks(RUN, s) if AT == "spatial" else set()) & set(KS)
    runs, stab, pa, nsplit, kdone = [], [], [], 0, 0
    if reuse:
        d = seed_root(RUN, s) / OLD_SUB
        r0 = pd.read_csv(d / "heldout_topm_crossfit_runs.csv"); r0 = r0[r0.K.isin(reuse)]
        k0 = pd.read_csv(d / "heldout_per_archetype_ranking_runs.csv"); k0 = k0[k0.K.isin(reuse)]
        runs.append(r0); stab.append(stability_long(k0)); pa.append(per_arch(k0))
        nsplit += len(reuse) * 3 * NREP; kdone += len(reuse)
    todo = [(AT, s, K) for K in KS if K not in reuse]
    with ProcessPoolExecutor(WORKERS) as ex:
        for K, n, r, st, p in ex.map(one_k, todo, chunksize=4):
            nsplit += n; kdone += n == 3 * NREP
            if n:
                runs.append(r); stab.append(st); pa.append(p)
    expect = len(KS) * 3 * NREP
    runs = pd.concat(runs, ignore_index=True)
    # keep only splits with both directions (as hseed_merge / run_parallel merge do)
    runs = runs[runs.groupby(KEY[:5])["direction"].transform("nunique") == 2]
    runs = runs.drop_duplicates(KEY + ["top_m"], keep="last").sort_values(["condition", "K", "outer_rep", "direction", "top_m"])
    runs.to_csv(out / "heldout_topm_crossfit_runs.csv.gz", index=False)
    write_heldout_summaries(runs.reset_index(drop=True), pd.DataFrame(), out, NREP)
    st = pd.concat(stab, ignore_index=True)
    rows = []
    for m in sorted(set(runs.top_m.astype(int))):
        x = st[st.rank_position <= m].copy()
        x["top_m"] = np.minimum(m, x.K)
        rows.append(x)
    st = pd.concat(rows).drop_duplicates(KEY + ["top_m", "archetype"])  # m > K collapses onto m = K
    stab = st.groupby(["analysis_type", "fit_scope", "condition", "K", "top_m", "archetype"], as_index=False).agg(
        selected_count=("archetype", "size"), mean_rank_position=("rank_position", "mean"))
    stab["selection_frequency"] = stab.selected_count / (NREP * 2)
    stab.sort_values(["condition", "K", "top_m", "selection_frequency", "archetype"],
                     ascending=[True, True, True, False, True]).to_csv(out / "heldout_archetype_selection_stability.csv", index=False)
    pd.concat(pa, ignore_index=True).to_csv(out / "heldout_per_archetype_mean.csv.gz", index=False)
    ok = nsplit == expect
    msg = f"seed {s}: {nsplit}/{expect} outer splits, {kdone}/{len(KS)} K complete ({len(reuse)} K reused from the 30-K run)"
    if ok:
        (out / "COMPLETE").write_text(msg + "\n")
    print(msg, flush=True)
    return out, ok, msg


def slim(s, K):
    p = (RUN / "seeds" / "slim" / "seed00" if s == 0 else seed_root(RUN, s) / "slim") / f"spatialAA_across_acrossCond_K{K}.npz"
    return np.load(p)["Smean"].astype(float)


def corr(A, B):
    A = (A - A.mean(1, keepdims=True)) / A.std(1, keepdims=True)
    B = (B - B.mean(1, keepdims=True)) / B.std(1, keepdims=True)
    return A @ B.T / A.shape[1]


def main():
    tag = f"heldout-everyk-{AT}-{'seed0' if PH in ('s0', 't0') else 'all'}" + ("-smoke" if MODE == "smoke" else "")
    OUT = RUN / f"{tag}-results"
    shutil.rmtree(OUT, ignore_errors=True); OUT.mkdir()
    ok_all, lines, summ, stab, allruns = True, [], [], [], []
    for s in SEEDS:
        d, ok, msg = merge_seed(s)
        ok_all &= ok; lines.append(msg)
        a = pd.read_csv(d / "heldout_topm_summary.csv").assign(seed=s); summ.append(a)
        b = pd.read_csv(d / "heldout_archetype_selection_stability.csv").assign(seed=s); stab.append(b)
        allruns.append(pd.read_csv(d / "heldout_topm_crossfit_runs.csv.gz").assign(fit_seed=s))
        for f in ("heldout_topm_summary.csv", "heldout_archetype_selection_stability.csv", "heldout_per_archetype_mean.csv.gz"):
            shutil.copy(d / f, OUT / f"seed{s:02d}_{f}")
    summ, stab = pd.concat(summ), pd.concat(stab)
    pd.concat(allruns, ignore_index=True).to_csv(OUT / "heldout_topm_runs_all_seeds.csv.gz", index=False)
    g = summ.groupby(["condition", "K", "top_m"])
    per_K = g["mean"].agg(n_seeds="count", mean="mean", sd="std", min="min", max="max").reset_index()
    per_K = per_K.merge(summ[summ.seed == 0][["condition", "K", "top_m", "mean"]].rename(columns={"mean": "seed0"}), how="left")
    per_K = per_K.merge(g[["mean_full_same_eval", "mean_retention"]].mean().reset_index())
    per_K.insert(0, "analysis_type", AT)
    per_K.to_csv(OUT / "heldout_per_K.csv", index=False)

    if AT == "spatial":
        rows, nomatch = [], 0
        for K in KS:
            try:
                S0 = slim(0, K)
            except FileNotFoundError:
                nomatch += 1; continue
            from scipy.optimize import linear_sum_assignment
            for s in SEEDS:
                try:
                    R = corr(S0, slim(s, K))
                except FileNotFoundError:
                    nomatch += 1; continue
                i, j = linear_sum_assignment(-R); to0 = dict(zip(j, i)); r = dict(zip(j, R[i, j]))
                t = stab[(stab.seed == s) & (stab.K == K)]
                rows += [dict(K=K, condition=x.condition, top_m=int(x.top_m), seed=s, archetype=int(x.archetype),
                              archetype0=int(to0[int(x.archetype)]), map_r=float(r[int(x.archetype)]),
                              selection_frequency=x.selection_frequency) for x in t.itertuples()]
        sel = pd.DataFrame(rows); sel.to_csv(OUT / "selection_by_seed.csv", index=False)
        if nomatch:
            lines.append(f"selection_by_seed: {nomatch} (seed, K) without slim maps, left out")
        k46 = sel[(sel.K == 46) & (sel.top_m == 5)] if len(sel) else sel
        if len(k46):
            t = k46.pivot_table(index=["archetype0", "condition"], columns="seed", values="selection_frequency", fill_value=0)
            t["mean"] = t.mean(axis=1); t["seeds_ge_0.5"] = (t.drop(columns="mean") >= 0.5).sum(axis=1)
            t = t.sort_values("mean", ascending=False); t.to_csv(OUT / "K46_top5.csv")
            lines += ["", "K=46 top-5 held-out selection frequency (seed-0 IDs; * = Fig 5 archetype), top 12 rows:",
                      t.head(12).round(2).rename(index=lambda a: f"{a}*" if a in PAPER else a, level=0).to_string()]
    for c in CONDS:
        for m in (5, 10):
            x = per_K[(per_K.condition == c) & (per_K.top_m == m) & (per_K.K > m)]
            if len(x):
                b = x.loc[x["mean"].idxmax()]
                lines.append(f"{c} top-{m}: peak of seed-mean at K={int(b.K)} ({b['mean']:.4f}); "
                             f"median seed SD across K {x.sd.median():.4f}")
    head = "SUMMARY COMPLETE" if ok_all else "SUMMARY INCOMPLETE"
    (OUT / "summary.txt").write_text("\n".join(
        [f"{head} ({MODE}, {AT}, seeds {SEEDS[0]}-{SEEDS[-1]}, K {KS[0]}-{KS[-1]}, {len(KS)} K)"] + lines) + "\n")
    print((OUT / "summary.txt").read_text(), flush=True)
    with tarfile.open(OUT.with_suffix(".tar.gz"), "w:gz") as tf:
        tf.add(OUT, arcname=OUT.name)
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
