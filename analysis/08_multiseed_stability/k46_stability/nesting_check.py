"""Raw nesting sanity check BEFORE consensus aggregation and alignment.

For every orientation (spatial/temporal, where results exist) x K x condition x fit seed x outer split x direction:
  assert set(top5_raw_ids) <= set(top10_raw_ids)   (raw, seed-own component IDs; no Hungarian, no consensus)
and additionally that top-5 is the ordered prefix of top-10 (both are ranked[:m] of one ranking).
Requested m is clamped to K (m = min(m, K)), as in heldout.evaluate_topm_on_subjects.
Also checks: (a) every (K, condition, seed) has all 40 splits (20 outer reps x 2 directions);
(b) the per-seed selection_stability CSVs (input to the consensus sets) reproduce the selected counts recomputed from
the raw runs; (c) reports, without asserting, how often a seed's consensus top-5 is not inside its consensus top-10.
Usage: python nesting_check.py [--root RUN_ROOT] [--out OUTDIR]
Reads RUN_ROOT/heldout-everyk-{spatial,temporal}-all-results/ (default: the run root in paths.py);
default OUTDIR outputs/multiseed_stability/nesting_check.
"""
import argparse, re, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis/08_multiseed_stability, for paths.py
import paths  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--root", default=str(paths.ROOT))
ap.add_argument("--out", default=str(paths.OUT / "nesting_check")); a = ap.parse_args()
HE, OUT = Path(a.root), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
log = open(OUT / "nesting_summary.txt", "w")
def say(*a):
    s = " ".join(map(str, a)); print(s, flush=True); log.write(s + "\n"); log.flush()
parse = lambda t: [int(x) for x in re.findall(r"\d+", t)]
allrows = []
for orient in ["spatial", "temporal"]:
    d = HE / f"heldout-everyk-{orient}-all-results"
    f = d / "heldout_topm_runs_all_seeds.csv.gz"
    if not f.exists(): say(orient, "no results, skipped"); continue
    cols = ["top_m", "selected_archetypes", "condition", "K", "outer_rep", "direction", "fit_seed"]
    cR = OUT / f".cache_runs_{orient}.pkl"
    if cR.exists(): R = pd.read_pickle(cR)
    else:
        R = pd.concat([ch[ch.top_m <= 10] for ch in pd.read_csv(f, usecols=cols, chunksize=500000)]); R.to_pickle(cR)
    R = R[(R.top_m == np.minimum(5, R.K)) | (R.top_m == np.minimum(10, R.K))].copy()
    R["want"] = np.where(R.top_m == np.minimum(10, R.K), 10, 5)
    # K <= 5: one row serves both
    small = R[R.K <= 5].copy(); small["want"] = 5; R = pd.concat([R, small]).drop_duplicates(["K", "condition", "fit_seed", "outer_rep", "direction", "want"])
    R["ids"] = R.selected_archetypes.map(parse)
    key = ["K", "condition", "fit_seed", "outer_rep", "direction"]
    W = R.pivot_table(index=key, columns="want", values="ids", aggfunc="first")
    W = W.dropna()
    n = len(W); nK = W.index.get_level_values("K").nunique()
    sub = W.apply(lambda r: set(r[5]) <= set(r[10]), axis=1)
    pre = W.apply(lambda r: r[10][:len(r[5])] == r[5], axis=1)
    sizes = W.apply(lambda r: len(r[5]) == min(5, r.name[0]) and len(r[10]) == min(10, r.name[0]) and len(set(r[10])) == len(r[10]), axis=1)
    cnt = W.groupby(["K", "condition", "fit_seed"]).size()
    say(f"[{orient}] rankings checked: {n}  (K values {nK}: {W.index.get_level_values('K').min()}-{W.index.get_level_values('K').max()}, "
        f"conditions {sorted(W.index.get_level_values('condition').unique())}, seeds {sorted(W.index.get_level_values('fit_seed').unique())})")
    say(f"[{orient}] top5 subset of top10: {int(sub.sum())}/{n} violations {int((~sub).sum())}")
    say(f"[{orient}] top5 == ordered prefix of top10: {int(pre.sum())}/{n}")
    say(f"[{orient}] set sizes correct & no duplicate ids: {int(sizes.sum())}/{n}")
    say(f"[{orient}] splits per (K,cond,seed): min {cnt.min()} max {cnt.max()} (expected 40); groups {len(cnt)}")
    assert sub.all() and pre.all() and sizes.all() and (cnt == 40).all(), orient
    if (~sub).any(): W[~sub].to_csv(OUT / f"violations_{orient}.csv")
    # (b) stability CSV reproduces raw counts
    raw = (R.explode("ids").groupby(["K", "condition", "fit_seed", "want", "ids"]).size().rename("raw_count").reset_index()
           .rename(columns={"ids": "archetype", "fit_seed": "seed"}))
    raw["archetype"] = raw.archetype.astype(int)
    mism, tot = 0, 0; cons_rows = []
    for s in sorted(R.fit_seed.unique()):
        p = d / f"seed{s:02d}_heldout_archetype_selection_stability.csv"
        cF = OUT / f".cache_stab_{orient}_{s}.pkl"
        if cF.exists(): F = pd.read_pickle(cF)
        else:
            F = pd.concat([ch[ch.top_m == np.minimum(np.where(ch.top_m <= 5, 5, 10), ch.K)] for ch in pd.read_csv(p, chunksize=500000)]); F.to_pickle(cF)
        F["want"] = np.where((F.top_m == np.minimum(10, F.K)) & (F.K > 5), 10, 5)
        sm = F[F.K <= 5].copy(); sm["want"] = 10; F = pd.concat([F, sm]).drop_duplicates(["K", "condition", "want", "archetype"])
        m = F.merge(raw[raw.seed == s], on=["K", "condition", "want", "archetype"], how="outer")
        bad = (m.selected_count.fillna(0) != m.raw_count.fillna(0)); mism += int(bad.sum()); tot += len(m)
        # (c) consensus nesting: head(m) by selection_frequency desc, mean_rank_position asc (as in topm_stability_neighbourhood.py)
        for (K, c), g in F.groupby(["K", "condition"]):
            top = {w: set(gg.sort_values(["selection_frequency", "mean_rank_position"], ascending=[False, True]).head(min(w, K)).archetype)
                   for w, gg in g.groupby("want")}
            if 5 in top and 10 in top:
                cons_rows.append(dict(orientation=orient, K=K, condition=c, seed=s, consensus_top5_in_top10=top[5] <= top[10],
                                      n_top5_outside_top10=len(top[5] - top[10]), outside=sorted(top[5] - top[10])))
    say(f"[{orient}] selection_stability CSV vs raw recount: {tot - mism}/{tot} (archetype x K x cond x seed x m) match, {mism} mismatch")
    C = pd.DataFrame(cons_rows); allrows.append(C)
    say(f"[{orient}] CONSENSUS (after aggregation over 40 splits, not asserted): seed top-5 inside seed top-10 in "
        f"{int(C.consensus_top5_in_top10.sum())}/{len(C)} (K x cond x seed); max outside {C.n_top5_outside_top10.max()}")
    for K in [46, 51, 14]:
        cc = C[C.K == K]
        if len(cc): say(f"   K={K}: {int(cc.consensus_top5_in_top10.sum())}/{len(cc)}  " + "; ".join(
            f"{r.condition} s{r.seed} {r.outside}" for r in cc[~cc.consensus_top5_in_top10].itertuples()))
pd.concat(allrows).to_csv(OUT / "consensus_nesting_by_seed.csv", index=False)
say("ALL RAW NESTING ASSERTIONS PASSED")
