#!/usr/bin/env python3
"""Parallel driver for the from-scratch rerun (not part of the repo).

Low-rank workspace (2026-10-27): reuses the validated rerun's fits (linked or
renumbered by rerun/reorder_fits.py), so only decode, heldout and merge run.
decode_subjects now uses the exact low-rank path for <= 10 archetypes.

Calls the unchanged msaa_pipeline functions, one task per (job, K[, condition]),
so results are identical to scripts/run_pipeline.py --profile full: every seed
depends only on K / outer_rep / direction, never on execution order.

Decode and held-out tasks each write into their own shard directory (to avoid
concurrent appends to one CSV); `merge` concatenates shards into the canonical
paths under the rerun root.

Usage (from the rerun root):
  python rerun/run_parallel.py fit     --jobs spatial_across [--workers 7 --threads 2]
  python rerun/run_parallel.py decode  --jobs spatial_across,temporal_across
  python rerun/run_parallel.py heldout
  python rerun/run_parallel.py merge
"""
from __future__ import annotations
import argparse, os, sys, time, traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHARDS = ROOT / "rerun" / "shards"
FIT_DIR = ROOT / "msaa_flexible_outputs_npz"


def _set_threads(n):
    for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        os.environ[v] = str(n)


def _init(threads):
    _set_threads(threads)
    os.chdir(ROOT)  # heldout reads configs/heldout_outer_seeds.csv relative to cwd
    sys.path.insert(0, str(ROOT / "src"))


def _log(msg):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


# ---------------------------------------------------------------- tasks
def fit_task(analysis_type, fit_scope, K, condition, maxiter):
    from msaa_pipeline.data import load_inputs
    from msaa_pipeline.fit import fit_one, fit_filename, save_fit, load_frozen_init
    out = FIT_DIR / fit_filename(analysis_type, fit_scope, K, condition)
    if out.exists():
        return f"skip {out.name}"
    condition_data, centers, widths = load_inputs(ROOT, trim_rest_by=100, max_subjects_per_condition=None)
    if fit_scope == "across":
        subjects, labels = [], []
        for c in ["intact", "word", "rest"]:
            subjects.extend(condition_data[c]); labels.extend([c] * len(condition_data[c]))
    else:
        subjects = condition_data[condition]; labels = [condition] * len(subjects)
    t = time.time()
    # Sequential archetype numbering: frozen FurthestSum picks and no post-fit sort,
    # so archetype j is the one initialised at init_index[j] on every machine.
    init_idx = load_frozen_init(ROOT, analysis_type, fit_scope, K, condition)
    res = fit_one(subjects, K, analysis_type, ROOT, maxiter=maxiter, seed=42, verbose=False,
                  init_idx=init_idx, sort_crit="none")
    tmp = out.with_name(out.stem + ".partial.npz")
    save_fit(tmp, res, labels, analysis_type, fit_scope, centers, widths)
    tmp.rename(out)
    return f"fit {out.name} varexpl={res['varexpl']:.6f} secs={res['secs']:.0f} wall={time.time()-t:.0f}"


def _shard(kind, analysis_type, fit_scope, K):
    d = SHARDS / kind / f"{analysis_type}_{fit_scope}_K{K}"
    d.mkdir(parents=True, exist_ok=True)
    link = d / "msaa_flexible_outputs_npz"
    if not link.exists():
        link.symlink_to(FIT_DIR)
    return d


def decode_task(analysis_type, fit_scope, K, nfolds, nreps):
    from msaa_pipeline.decode import run_decode_job
    d = _shard("decode", analysis_type, fit_scope, K)
    run_decode_job(d, analysis_type, fit_scope, [K], nfolds=nfolds, nreps=nreps, seed=42)
    return f"decode {analysis_type}_{fit_scope} K={K}"


def heldout_task(analysis_type, K, top_m, n_outer, nfolds, nreps):
    from msaa_pipeline.heldout import run_heldout_topm_job
    d = _shard("heldout", analysis_type, "across", K)
    run_heldout_topm_job(d, analysis_type, [K], top_m=top_m, n_outer_reps=n_outer,
                         nfolds=nfolds, nreps=nreps, seed=42)
    return f"heldout {analysis_type} K={K}"


def _split_done(K, condition, outer_rep):
    """True if this split is complete in the per-K shard (old layout) or its own split file."""
    import pandas as pd
    if (SHARDS / "heldout_split" / f"K{K}_{condition}_r{outer_rep}" / "runs.csv").exists():
        return True
    p = SHARDS / "heldout" / f"spatial_across_K{K}" / "msaa_heldout_topm_outputs_spatial_across" / "heldout_topm_crossfit_runs.csv"
    q = p.with_name("heldout_per_archetype_ranking_runs.csv")
    if not (p.exists() and q.exists()):
        return False
    for f in (p, q):
        d = pd.read_csv(f, usecols=["condition", "outer_rep", "direction"])
        d = d[(d.condition == condition) & (d.outer_rep == outer_rep)]
        if not {"A_to_B", "B_to_A"}.issubset(set(d.direction)):
            return False
    return True


def heldout_split_task(analysis_type, K, condition, outer_rep, top_m, nfolds, nreps, seed=42):
    """One outer split of run_heldout_topm_job, with the identical seeds and steps."""
    import pandas as pd
    import numpy as np
    from msaa_pipeline.decode import get_subjects
    from msaa_pipeline.heldout import (_load_frozen_outer_seeds, _split_outer_subject_indices,
                                       _subset_subjects, rank_archetypes_on_subjects,
                                       evaluate_topm_on_subjects)
    d = SHARDS / "heldout_split" / f"K{K}_{condition}_r{outer_rep}"
    if (d / "runs.csv").exists():
        return f"skip heldout K={K} {condition} rep={outer_rep}"
    subjects, source = get_subjects(ROOT, analysis_type, "across", K, condition)
    outer_seed = int(_load_frozen_outer_seeds()[(str(condition), int(K), int(outer_rep))])
    A_idx, B_idx = _split_outer_subject_indices(len(subjects), np.random.default_rng(outer_seed))
    new_runs, new_ranks = [], []
    for direction, ranking_idx, evaluation_idx in [("A_to_B", A_idx, B_idx), ("B_to_A", B_idx, A_idx)]:
        ranking_subjects = _subset_subjects(subjects, ranking_idx)
        evaluation_subjects = _subset_subjects(subjects, evaluation_idx)
        direction_offset = 0 if direction == "A_to_B" else 1
        ranking_seed = seed + 100000 * outer_rep + 1000 * K + direction_offset
        evaluation_seed = seed + 200000 * outer_rep + 1000 * K + direction_offset
        ranked, per_arch = rank_archetypes_on_subjects(ranking_subjects, analysis_type,
                                                       nfolds=nfolds, nreps=nreps, seed=ranking_seed)
        per_arch = per_arch.copy()
        per_arch["analysis_type"] = analysis_type
        per_arch["fit_scope"] = "across"
        per_arch["condition"] = condition
        per_arch["K"] = K
        per_arch["outer_rep"] = outer_rep
        per_arch["direction"] = direction
        per_arch["outer_seed"] = outer_seed
        per_arch["n_ranking_subjects"] = len(ranking_subjects)
        new_ranks.append(per_arch)
        heldout = evaluate_topm_on_subjects(evaluation_subjects, analysis_type, ranked,
                                            top_m_values=top_m, nfolds=nfolds, nreps=nreps,
                                            seed=evaluation_seed)
        heldout["analysis_type"] = analysis_type
        heldout["fit_scope"] = "across"
        heldout["condition"] = condition
        heldout["K"] = K
        heldout["outer_rep"] = outer_rep
        heldout["direction"] = direction
        heldout["outer_seed"] = outer_seed
        heldout["n_ranking_subjects"] = len(ranking_subjects)
        heldout["n_evaluation_subjects"] = len(evaluation_subjects)
        heldout["fit_source_file"] = str(source)
        heldout["fit_source_basename"] = source.name
        new_runs.append(heldout)
    d.mkdir(parents=True, exist_ok=True)
    pd.concat(new_ranks, ignore_index=True).to_csv(d / "ranks.csv", index=False)
    pd.concat(new_runs, ignore_index=True).to_csv(d / "runs.partial.csv", index=False)
    (d / "runs.partial.csv").rename(d / "runs.csv")  # runs.csv marks the split complete
    return f"heldout K={K} {condition} rep={outer_rep}"


def _run(tasks, workers, threads):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    ctx = mp.get_context("spawn")
    _log(f"{len(tasks)} tasks, {workers} workers x {threads} threads")
    failed = 0
    with ProcessPoolExecutor(workers, mp_context=ctx, initializer=_init, initargs=(threads,)) as ex:
        futs = {ex.submit(fn, *args): (fn.__name__, args) for fn, args in tasks}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                _log(f"[{i}/{len(tasks)}] {f.result()}")
            except Exception:
                failed += 1
                _log(f"[{i}/{len(tasks)}] FAILED {futs[f]}\n{traceback.format_exc()}")
    _log(f"done, {failed} failed")
    return failed


# ---------------------------------------------------------------- merge
def merge():
    import pandas as pd
    _init(1)
    from msaa_pipeline.heldout import write_heldout_summaries
    from msaa_pipeline.profiles import load_profile
    prof = load_profile(ROOT, "full")

    def cat(paths):
        dfs = [pd.read_csv(p) for p in paths if p.exists()]
        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    def fix_src(df):
        if "fit_source_file" in df:
            df["fit_source_file"] = df["fit_source_basename"].map(lambda b: str(FIT_DIR / b))
        return df

    for job in prof["jobs"]:
        at, sc = job.split("_")
        shards = sorted((SHARDS / "decode").glob(f"{at}_{sc}_K*"))
        if not shards:
            continue
        out = ROOT / f"msaa_condrank_decoding_outputs_{at}_{sc}"
        out.mkdir(exist_ok=True)
        full = fix_src(cat([s / f"msaa_condrank_decoding_outputs_{at}_{sc}" / "full_summary.csv" for s in shards]))
        arch = fix_src(cat([s / f"msaa_condrank_decoding_outputs_{at}_{sc}" / "per_archetype_mean_accuracy.csv" for s in shards]))
        full.sort_values(["analysis_type", "fit_scope", "condition", "K"]).to_csv(out / "full_summary.csv", index=False)
        arch.sort_values(["analysis_type", "fit_scope", "condition", "K", "archetype"]).to_csv(out / "per_archetype_mean_accuracy.csv", index=False)
        _log(f"merged decode {job}: {len(shards)} K, {len(full)} full rows, {len(arch)} archetype rows")

    shards = sorted((SHARDS / "heldout").glob("spatial_across_K*"))
    splits = sorted((SHARDS / "heldout_split").glob("K*_r*"))
    if shards or splits:
        out = ROOT / "msaa_heldout_topm_outputs_spatial_across"
        out.mkdir(exist_ok=True)
        sub = "msaa_heldout_topm_outputs_spatial_across"
        runs = fix_src(cat([s / sub / "heldout_topm_crossfit_runs.csv" for s in shards] + [s / "runs.csv" for s in splits]))
        ranks = cat([s / sub / "heldout_per_archetype_ranking_runs.csv" for s in shards] + [s / "ranks.csv" for s in splits if (s / "runs.csv").exists()])
        key = ["analysis_type", "fit_scope", "condition", "K", "outer_rep", "direction"]
        # Old per-K shards may hold a split with only one direction (killed mid-split); keep complete splits only.
        full = ranks.groupby(key[:5])["direction"].transform("nunique") == 2
        ranks = ranks[full]
        full = runs.groupby(key[:5])["direction"].transform("nunique") == 2
        runs = runs[full]
        runs = runs.drop_duplicates(key + ["top_m"], keep="last")
        ranks = ranks.drop_duplicates(key + ["archetype"], keep="last")
        runs = runs.sort_values(["condition", "K", "outer_rep", "direction", "top_m"]).reset_index(drop=True)
        ranks = ranks.sort_values(["condition", "K", "outer_rep", "direction", "archetype"]).reset_index(drop=True)
        runs.to_csv(out / "heldout_topm_crossfit_runs.csv", index=False)
        ranks.to_csv(out / "heldout_per_archetype_ranking_runs.csv", index=False)
        write_heldout_summaries(runs, ranks, out, prof["heldout_outer_reps"])
        _log(f"merged heldout: {len(runs.groupby(['condition', 'K', 'outer_rep']))} splits, {len(runs)} run rows")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["fit", "decode", "heldout", "merge"])
    ap.add_argument("--jobs", default=None, help="comma list of profile job keys (default: all)")
    ap.add_argument("--ks", default=None, help="comma list restricting K values")
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--threads", type=int, default=1)
    # ben 2026-10-27: use a tie tolerance so results agree across machines/BLAS
    # (temporal S is mostly exact zeros, giving exact correlation ties).
    ap.add_argument("--tie-tol", type=float, default=1e-8,
                    help="correlation tie tolerance for decode scoring (msaa_pipeline.decode.TIE_TOL); "
                         "0 = exact comparison as in the validated rerun")
    args = ap.parse_args()
    os.environ["MSAA_TIE_TOL"] = repr(args.tie_tol)  # inherited by the spawned workers
    _log(f"MSAA_TIE_TOL={os.environ.get('MSAA_TIE_TOL', '0')}")

    if args.stage == "merge":
        return merge()

    _init(args.threads)
    from msaa_pipeline.profiles import load_profile
    prof = load_profile(ROOT, "full")
    keep_k = set(map(int, args.ks.split(","))) if args.ks else None
    tasks = []
    if args.stage in ("fit", "decode"):
        jobs = args.jobs.split(",") if args.jobs else list(prof["jobs"])
        for job in jobs:
            at, sc = job.split("_")
            for K in prof["jobs"][job]:
                if keep_k and K not in keep_k:
                    continue
                if args.stage == "fit":
                    conds = [None] if sc == "across" else ["intact", "word", "rest"]
                    for c in conds:
                        tasks.append((K, (fit_task, (at, sc, K, c, prof["fit_maxiter"]))))
                else:
                    tasks.append((K, (decode_task, (at, sc, K, prof["decode_nfolds"], prof["decode_nreps"]))))
    else:
        # profiles.json lists 40 held-out K but configs/heldout_outer_seeds.csv
        # (and the legacy outputs) only cover 30; the others would raise.
        import pandas as pd
        seeded = set(pd.read_csv(ROOT / "configs" / "heldout_outer_seeds.csv")["K"].astype(int))
        for K in prof["heldout_jobs"]["spatial_across"]:
            if keep_k and K not in keep_k:
                continue
            if K not in seeded:
                _log(f"skipping held-out K={K}: no frozen outer seeds")
                continue
            for condition in ["intact", "word", "rest"]:
                for rep in range(prof["heldout_outer_reps"]):
                    if not _split_done(K, condition, rep):
                        tasks.append((K, (heldout_split_task, ("spatial", K, condition, rep, prof["top_m"],
                                                               prof["decode_nfolds"], prof["decode_nreps"]))))
    # Longest (largest K) first so the big tasks don't straggle at the end.
    tasks = [t for _, t in sorted(tasks, key=lambda x: -x[0])]
    sys.exit(1 if _run(tasks, args.workers, args.threads) else 0)


if __name__ == "__main__":
    main()
