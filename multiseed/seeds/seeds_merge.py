#!/usr/bin/env python3
"""Summary step of the multi-seed run. Reads every seed directory that exists (all stages), plus
seed 0 (the existing every-K fits in $RUN_DIR), and writes $RUN_DIR/seeds/results/:

  fits.csv              one row per fit x seed x decoding condition: final loss, varexpl,
                        iterations, full-K decoding, best single archetype, top-m (across only)
  per_K.csv             per analysis/scope/condition/K and metric: mean, SD, min, max over seeds,
                        the canonical fit's value (seed 0), the lowest-loss fit's value and seed,
                        and seed 0's loss rank (1 = lowest)
  matching_fits.csv     per fit and seed >= 1: Hungarian match of its archetypes to seed 0's
                        (spatial: subject-mean maps; temporal: C over time, best corr within
                        +-LAG TRs), with summaries of the matched correlations, how far above the
                        next-best archetype they are, and whether per-archetype decoding reproduces
  matching_archetypes.csv.gz  across fits only: one row per seed-0 archetype x seed
  paper_K46.csv         the K=46 spatial-across paper archetypes (28, 35, 37, 31, 39) across seeds
  jumpiness.csv         mean |change| from one K to the next for seed 0, the seed mean and the
                        lowest-loss fit, on the K values every seed has
  summary.txt           the headline numbers
and packs them into $RUN_DIR/seeds-results.tar.gz.

  python seeds_merge.py [--workers 16] [--lag 2]
Exit code 1 (and "INCOMPLETE" in the summary) if any expected seed fit, decode table or top-m part
is missing; the tables are still written for what exists.
"""
import argparse
import json
import os
import sys
import tarfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from seeds_common import (CONDITIONS, ROOT, SEEDS_DIR, TOP_M_VALUES, all_seed_dirs, fit_name,  # noqa: E402
                          stage_Ks, write_slim)

HERE = Path(__file__).resolve().parent
OUT = SEEDS_DIR / "results"
SLIM0 = SEEDS_DIR / "slim" / "seed00"
PAPER_K46 = [28, 35, 37, 31, 39]
LAG = 2


# ---------- inputs ----------

def parse_name(name):
    """spatialAA_within_rest_K5.npz -> (spatial, within, rest, 5); across -> cond None."""
    stem = name[:-4]
    at, rest = stem.split("AA_", 1)
    sc, cond, k = rest.split("_")
    return at, sc, (None if sc == "across" else cond), int(k[1:])


def included_stages():
    """Stage this summary covers: smoke or full (the STAGE being merged; without STAGE, full if it
    has a task setup, else smoke)."""
    st = os.environ.get("STAGE")
    if st in ("smoke", "full"):
        return [st]
    return ["full"] if (HERE / "tasks" / "full_expected.json").exists() else ["smoke"]


def seed_dirs():
    inc = included_stages()
    return [(st, s, d) for st, s, d in all_seed_dirs() if st in inc]


def seed_fit_index():
    """{(seed, at, sc, cond, K): (stage, slim path, seed dir)} for every finished seed fit."""
    idx = {}
    for stage, s, d in seed_dirs():
        for p in (d / "slim").glob("*.npz"):
            if ".tmp" in p.name:
                continue
            at, sc, cond, K = parse_name(p.name)
            idx[(s, at, sc, cond, K)] = (stage, p, d)
    return idx


def expected_fits():
    """{(stage, seed, at, sc, K, cond)} that the included stages' task setups expect."""
    want = set()
    for stage in included_stages():
        exp = json.loads((HERE / "tasks" / f"{stage}_expected.json").read_text())
        for s in exp["seeds"]:
            for job in exp["jobs"]:
                at, sc = job.split("_")
                for K in stage_Ks(stage, at):
                    for cond in ([None] if sc == "across" else CONDITIONS):
                        want.add((stage, s, at, sc, cond, K))
    return want


def _slim0(name):
    out = SLIM0 / name
    if not out.exists():
        write_slim(ROOT / "msaa_flexible_outputs_npz" / name, out)
    return name


def load_slim(p):
    z = np.load(p)
    return {k: z[k] for k in z.files}


def decode_tables(root, at, sc):
    """(full-K decoding mean per (condition, K), {(condition, K): per-archetype accuracy by archetype})."""
    d = Path(root) / f"msaa_condrank_decoding_outputs_{at}_{sc}"
    full = pa = None
    if (d / "full_summary.csv").exists():
        f = pd.read_csv(d / "full_summary.csv")
        full = f.groupby(["condition", "K"])["mean"].mean()
    if (d / "per_archetype_mean_accuracy.csv").exists():
        p = pd.read_csv(d / "per_archetype_mean_accuracy.csv")
        p = p.groupby(["condition", "K", "archetype"])["mean_accuracy"].mean()
        pa = {(c, int(K)): v.droplevel([0, 1]) for (c, K), v in p.groupby(level=[0, 1])}
    return full, pa


def topm_table_seed(d):
    parts = sorted((d / "topm_parts").glob("*.csv")) if (d / "topm_parts").exists() else []
    if not parts:
        return None
    t = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    return t.groupby(["analysis_type", "condition", "K", "top_m"])["mean"].mean()


def topm_table_seed0():
    out = []
    for at in ["spatial", "temporal"]:
        p = ROOT / "topm-everyk" / f"msaa_condrank_decoding_outputs_{at}_across" / "topm_summary.csv"
        if p.exists():
            out.append(pd.read_csv(p))
    if not out:
        return None
    t = pd.concat(out, ignore_index=True)
    return t.groupby(["analysis_type", "condition", "K", "top_m"])["mean"].mean()


# ---------- matching ----------

def corr_matrix(A, B):
    """Column-by-column Pearson corr, A (n x Ka), B (n x Kb) -> (Ka x Kb)."""
    A = A - A.mean(0); B = B - B.mean(0)
    na, nb = np.linalg.norm(A, axis=0), np.linalg.norm(B, axis=0)
    R = (A / np.where(na > 0, na, 1)).T @ (B / np.where(nb > 0, nb, 1))
    R[np.outer(na == 0, nb == 0)] = 1.0
    return R


def shift(C, lag):
    """Shift columns of C (time x K) down by lag rows (lag < 0: up), zero-filled."""
    out = np.zeros_like(C)
    if lag > 0:
        out[lag:] = C[:-lag]
    elif lag < 0:
        out[:lag] = C[-lag:]
    else:
        out[:] = C
    return out


def hungarian(R):
    from scipy.optimize import linear_sum_assignment
    r, c = linear_sum_assignment(-R)
    m = np.empty(R.shape[0], int); m[r] = c
    return m


def match(at, a, b, lag):
    """Match seed b's archetypes to seed 0's (a). Returns per-archetype arrays for a's archetypes."""
    Ca, Cb = a["C"].astype(float), b["C"].astype(float)
    if at == "spatial":
        R = corr_matrix(a["Smean"].astype(float).T, b["Smean"].astype(float).T)   # maps over nodes
        m = hungarian(R)
        j = np.arange(len(m))
        out = dict(match=m, corr=R[j, m], C_corr=corr_matrix(Ca, Cb)[j, m])
    else:
        lags = list(range(-lag, lag + 1))
        Rl = np.stack([corr_matrix(Ca, shift(Cb, L)) for L in lags])                  # C over time
        best = Rl.max(0)
        m = hungarian(best)
        j = np.arange(len(m))
        out = dict(match=m, corr=best[j, m], corr_lag0=Rl[lags.index(0)][j, m],
                   lag=np.array(lags)[Rl.argmax(0)[j, m]],
                   Smean_corr=corr_matrix(a["Smean"].astype(float).T, b["Smean"].astype(float).T)[j, m])
        R = best
    Rn = R.copy(); Rn[np.arange(len(m)), m] = -np.inf
    out["next_best"] = Rn.max(1) if R.shape[1] > 1 else np.full(len(m), np.nan)
    return out


def match_job(args):
    key, p0, p1, lag = args
    s, at, sc, cond, K = key
    a, b = load_slim(p0), load_slim(p1)
    r = match(at, a, b, lag)
    return key, r, float(a["cost_last"]), float(b["cost_last"])


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--lag", type=int, default=LAG)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); SLIM0.mkdir(parents=True, exist_ok=True)
    notes = []

    idx = seed_fit_index()
    want = expected_fits()
    have = {(st,) + k for k, (st, _, _) in idx.items()}
    missing_fits = sorted(want - have, key=str)
    print(f"seed fits found: {len(idx)}; expected {len(want)}; missing {len(missing_fits)}", flush=True)

    # seed 0 slim copies for every (at, sc, cond, K) that has a seed fit
    keys0 = sorted({k[1:] for k in idx})
    names0 = [fit_name(at, sc, K, cond) for at, sc, cond, K in keys0]
    todo = [n for n in names0 if not (SLIM0 / n).exists()]
    print(f"seed 0 slim copies to make: {len(todo)} of {len(names0)}", flush=True)
    with ProcessPoolExecutor(args.workers) as ex:
        for i, _ in enumerate(ex.map(_slim0, todo, chunksize=4), 1):
            if i % 200 == 0:
                print(f"  seed 0 slim {i}/{len(todo)}", flush=True)

    # loss, varexpl, iterations for every fit (seed 0 included)
    info = []
    for (s, at, sc, cond, K), (st, p, d) in idx.items():
        z = np.load(p)
        info.append(dict(seed=s, stage=st, analysis_type=at, fit_scope=sc, fit_cond=cond or "acrossCond", K=K,
                         cost_last=float(z["cost_last"]), varexpl=float(z["varexpl"]), n_iter=int(z["n_iter"])))
    for at, sc, cond, K in keys0:
        z = np.load(SLIM0 / fit_name(at, sc, K, cond))
        info.append(dict(seed=0, stage="seed0", analysis_type=at, fit_scope=sc, fit_cond=cond or "acrossCond", K=K,
                         cost_last=float(z["cost_last"]), varexpl=float(z["varexpl"]), n_iter=int(z["n_iter"])))
    info = pd.DataFrame(info)

    # decoding tables per seed root (seed 0 = $RUN_DIR) and top-m
    roots = {(0, None): ROOT}
    for st, s, d in seed_dirs():
        roots[(s, st)] = d
    dec, missing_dec = {}, []
    for (s, st), r in roots.items():
        for at in ["spatial", "temporal"]:
            for sc in ["across", "within"]:
                full, pa = decode_tables(r, at, sc)
                if s and any(k[0] == s and k[1] == at and k[2] == sc and idx[k][0] == st for k in idx) and full is None:
                    missing_dec.append(f"seed{s:02d}/{st} {at}_{sc}")
                dec[(s, st, at, sc)] = (full, pa)
    topm = {(0, None): topm_table_seed0()}
    for st, s, d in seed_dirs():
        topm[(s, st)] = topm_table_seed(d)

    def stage_of(s, at, sc, cond, K):
        return None if s == 0 else idx[(s, at, sc, cond, K)][0]

    # fits.csv: one row per fit x seed x decoding condition
    rows = []
    for r in info.itertuples(index=False):
        cond = None if r.fit_cond == "acrossCond" else r.fit_cond
        st = stage_of(r.seed, r.analysis_type, r.fit_scope, cond, r.K)
        full, pa = dec.get((r.seed, st, r.analysis_type, r.fit_scope), (None, None))
        tm = topm.get((r.seed, st))
        for c in (CONDITIONS if r.fit_scope == "across" else [cond]):
            row = r._asdict(); row["condition"] = c
            row["full_mean"] = full.get((c, r.K), np.nan) if full is not None else np.nan
            if pa is not None and (c, r.K) in pa:
                v = pa[(c, r.K)]
                row["best_single"], row["mean_single"], row["best_archetype"] = v.max(), v.mean(), int(v.idxmax())
            else:
                row["best_single"] = row["mean_single"] = np.nan; row["best_archetype"] = -1
            for m in TOP_M_VALUES:
                row[f"top{m}"] = (tm.get((r.analysis_type, c, r.K, m), np.nan)
                                  if tm is not None and r.fit_scope == "across" else np.nan)
            rows.append(row)
    fits = pd.DataFrame(rows)
    fits.to_csv(OUT / "fits.csv", index=False)
    missing_topm = int(((fits.fit_scope == "across") & (fits.seed > 0) & fits.top1.isna()).sum())

    # per_K.csv
    metrics = ["full_mean", "best_single", "mean_single"] + [f"top{m}" for m in TOP_M_VALUES] + ["varexpl", "cost_last"]
    g = ["analysis_type", "fit_scope", "condition", "K"]
    per = []
    for key, sub in fits.groupby(g):
        sub = sub.sort_values("seed")
        lo = sub.loc[sub.cost_last.idxmin()]
        s0 = sub[sub.seed == 0]
        base = dict(zip(g, key), n_seeds=len(sub), lowest_loss_seed=int(lo.seed),
                    seed0_loss_rank=int(sub.cost_last.rank(method="min")[s0.index].iloc[0]) if len(s0) else -1,
                    loss_spread_rel=(sub.cost_last.max() - sub.cost_last.min()) / abs(sub.cost_last.min()))
        for mt in metrics:
            v = sub[mt]
            if v.notna().sum() == 0:
                continue
            per.append(dict(base, metric=mt, mean=v.mean(), sd=v.std(ddof=1) if v.notna().sum() > 1 else np.nan,
                            min=v.min(), max=v.max(), seed0=s0[mt].iloc[0] if len(s0) else np.nan,
                            lowest_loss=lo[mt], n=int(v.notna().sum())))
    per = pd.DataFrame(per)
    per.to_csv(OUT / "per_K.csv", index=False)

    # matching
    jobs = [(k, SLIM0 / fit_name(k[1], k[2], k[4], k[3]), p, args.lag) for k, (st, p, d) in idx.items()]
    mrows, arows = [], []
    with ProcessPoolExecutor(args.workers) as ex:
        for i, (key, r, c0, c1) in enumerate(ex.map(match_job, jobs, chunksize=8), 1):
            s, at, sc, cond, K = key
            st = idx[key][0]
            corr = r["corr"]
            row = dict(seed=s, stage=st, analysis_type=at, fit_scope=sc, fit_cond=cond or "acrossCond", K=K,
                       corr_median=np.median(corr), corr_mean=corr.mean(), corr_min=corr.min(),
                       frac_ge_090=(corr >= 0.9).mean(), frac_ge_099=(corr >= 0.99).mean(),
                       frac_above_next=np.nanmean(corr > r["next_best"]) if K > 1 else np.nan,
                       n_same_index=int((r["match"] == np.arange(K)).sum()),
                       loss_rel_to_seed0=(c1 - c0) / abs(c0))
            if at == "spatial":
                row["C_corr_median"] = np.median(r["C_corr"])
            else:
                row.update(corr_lag0_median=np.median(r["corr_lag0"]), frac_lag_nonzero=(r["lag"] != 0).mean(),
                           Smean_corr_median=np.median(r["Smean_corr"]))
            # does per-archetype decoding reproduce after matching?
            pa0 = dec[(0, None, at, sc)][1]; pas = dec.get((s, st, at, sc), (None, None))[1]
            accs, same = [], []
            for c in (CONDITIONS if sc == "across" else [cond]):
                if pa0 is None or pas is None or (c, K) not in pa0 or (c, K) not in pas:
                    continue
                v0 = pa0[(c, K)].reindex(range(K)).to_numpy(); v1 = pas[(c, K)].reindex(range(K)).to_numpy()
                v1m = v1[r["match"]]
                if np.isfinite(v0).all() and np.isfinite(v1m).all() and K > 2:
                    accs.append(np.corrcoef(v0, v1m)[0, 1])
                    same.append(float(r["match"][np.argmax(v0)] == np.argmax(v1)))
            row["acc_corr_matched"] = np.mean(accs) if accs else np.nan
            row["best_archetype_same"] = np.mean(same) if same else np.nan
            mrows.append(row)
            if sc == "across":
                a = dict(seed=s, analysis_type=at, K=K, archetype=np.arange(K), matched=r["match"],
                         corr=corr, next_best=r["next_best"])
                a.update({k: r[k] for k in ["C_corr", "corr_lag0", "lag", "Smean_corr"] if k in r})
                arows.append(pd.DataFrame(a))
            if i % 1000 == 0:
                print(f"  matched {i}/{len(jobs)}", flush=True)
    M = pd.DataFrame(mrows)
    M.to_csv(OUT / "matching_fits.csv", index=False)
    A = pd.concat(arows, ignore_index=True) if arows else pd.DataFrame()
    if len(A):
        A.to_csv(OUT / "matching_archetypes.csv.gz", index=False)
        P = A[(A.analysis_type == "spatial") & (A.K == 46) & A.archetype.isin(PAPER_K46)]
        P.to_csv(OUT / "paper_K46.csv", index=False)
    else:
        P = pd.DataFrame()

    # jumpiness: K-to-K change on the K values every seed has
    jump = []
    for key, sub in per[per.metric.isin(["full_mean", "best_single", "top1", "top3", "top5", "top10"])].groupby(
            ["analysis_type", "fit_scope", "condition", "metric"]):
        sub = sub[sub.n == sub.n.max()].sort_values("K")
        if len(sub) < 3:
            continue
        d = lambda col: np.nanmean(np.abs(np.diff(sub[col].to_numpy())))
        jump.append(dict(zip(["analysis_type", "fit_scope", "condition", "metric"], key), n_K=len(sub),
                         n_seeds=int(sub.n.max()), jump_seed0=d("seed0"), jump_seed_mean=d("mean"),
                         jump_lowest_loss=d("lowest_loss"), median_seed_sd=sub.sd.median()))
    J = pd.DataFrame(jump)
    J.to_csv(OUT / "jumpiness.csv", index=False)

    # summary.txt
    complete = not missing_fits and not missing_dec and missing_topm == 0
    L = [f"Multi-seed MS-AA run: {'COMPLETE' if complete else 'INCOMPLETE'}",
         f"seed fits: {len(idx)} found, {len(missing_fits)} missing; seeds {sorted({k[0] for k in idx})}; "
         f"stages {sorted({v[0] for v in idx.values()})}",
         f"decode tables missing: {len(missing_dec)} {missing_dec[:6]}", f"top-m rows missing (across, seeds >= 1): {missing_topm}", ""]
    if missing_fits:
        L += ["first missing fits: " + "; ".join(" ".join(map(str, f)) for f in missing_fits[:8]), ""]
    if len(per):
        lo = per[per.metric == "cost_last"]
        L += ["Seed 0 (current paper fit) loss rank among the seeds at each K (1 = lowest loss):",
              lo.groupby(["analysis_type", "fit_scope"]).seed0_loss_rank.describe()[["count", "mean", "50%", "max"]].to_string(),
              f"share of fits where seed 0 is the lowest-loss fit: {(lo.seed0_loss_rank == 1).mean():.2f}",
              "relative loss spread across seeds (max-min)/min, median by analysis: " +
              lo.groupby("analysis_type").loss_spread_rel.median().round(5).to_dict().__repr__(), ""]
    if len(M):
        M["K_band"] = pd.cut(M.K, [0, 20, 50, 100, 200, 400, 700])
        L += ["Matched archetype corr to seed 0 (median of per-fit medians; spatial = maps, temporal = C within "
              f"+-{args.lag} TRs):",
              M.groupby(["analysis_type", "fit_scope", "K_band"], observed=True)[["corr_median", "frac_ge_090",
              "frac_above_next", "acc_corr_matched", "best_archetype_same"]].median().round(3).to_string(), ""]
        if "frac_lag_nonzero" in M:
            L += [f"temporal: median share of matches at a nonzero lag {M.frac_lag_nonzero.median():.2f}", ""]
    if len(P):
        L += ["K=46 spatial-across paper archetypes vs their matches in each seed:",
              P.groupby("archetype")[["corr", "C_corr", "next_best"]].agg(["median", "min"]).round(3).to_string(), ""]
    if len(J):
        L += ["K-to-K jumpiness (mean |change| between neighbouring K on the common grid):",
              J.round(4).to_string(index=False), ""]
    L += [f"results: {OUT}"]
    text = "\n".join(L)
    (OUT / "summary.txt").write_text(text + "\n")
    print(text, flush=True)

    with tarfile.open(ROOT / "seeds-results.tar.gz", "w:gz") as tf:
        tf.add(OUT, arcname="seeds-results")
    print("wrote", ROOT / "seeds-results.tar.gz", flush=True)
    sys.exit(0 if complete else 1)


if __name__ == "__main__":
    main()
