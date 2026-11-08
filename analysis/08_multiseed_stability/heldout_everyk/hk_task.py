#!/usr/bin/env python3
"""One task: every unfinished outer split of its (analysis, seed, K list).
   python hk_task.py RUN_DIR smoke|full PHASE INDEX(1-based) WORKERS
The split itself is rerun/run_parallel.py heldout_split_task line for line (same msaa_pipeline functions,
ranking/evaluation seeds, nfolds/nreps, top-m), with the outer seed read from the every-K table and
analysis = spatial or temporal. Output per split: BASE/parts/<analysis>/seedNN/K<K>_<cond>_r<rep>/{ranks,runs}.csv
(runs.csv is written last and marks the split complete; reruns skip it)."""
import os
# one BLAS thread per worker: set before numpy is imported here or in the spawned workers
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[_v] = "1"
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from pathlib import Path

KIT = Path(__file__).resolve().parent
sys.path.insert(0, str(KIT))
from hk_common import BASE_SEED, CONDS, NFOLDS, NREP, NREPS, TOP_M, base_dir, load_outer_seeds, seed_root, split_dir  # noqa


def log(msg):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def init(run):
    for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        os.environ[v] = "1"
    sys.path.insert(0, str(Path(run) / "src"))
    sys.path.insert(0, str(KIT))


def unit(run, base, at, s, K, condition, reps):
    """Outer splits `reps` of one (seed, K, condition); loads the fit once."""
    import numpy as np
    import pandas as pd
    from msaa_pipeline.decode import get_subjects
    from msaa_pipeline.heldout import (_split_outer_subject_indices, _subset_subjects,
                                       rank_archetypes_on_subjects, evaluate_topm_on_subjects)
    todo = [r for r in reps if not (split_dir(base, at, s, K, condition, r) / "runs.csv").exists()]
    if not todo:
        return f"skip {at} seed{s} K={K} {condition}"
    seeds = load_outer_seeds(base)
    subjects, source = get_subjects(seed_root(run, s), at, "across", K, condition)
    for outer_rep in todo:
        d = split_dir(base, at, s, K, condition, outer_rep)
        outer_seed = seeds[(condition, K, outer_rep)]
        A_idx, B_idx = _split_outer_subject_indices(len(subjects), np.random.default_rng(outer_seed))
        new_runs, new_ranks = [], []
        for direction, ranking_idx, evaluation_idx in [("A_to_B", A_idx, B_idx), ("B_to_A", B_idx, A_idx)]:
            ranking_subjects = _subset_subjects(subjects, ranking_idx)
            evaluation_subjects = _subset_subjects(subjects, evaluation_idx)
            direction_offset = 0 if direction == "A_to_B" else 1
            ranking_seed = BASE_SEED + 100000 * outer_rep + 1000 * K + direction_offset
            evaluation_seed = BASE_SEED + 200000 * outer_rep + 1000 * K + direction_offset
            ranked, per_arch = rank_archetypes_on_subjects(ranking_subjects, at, nfolds=NFOLDS, nreps=NREPS,
                                                           seed=ranking_seed)
            per_arch = per_arch.copy()
            per_arch["analysis_type"] = at
            per_arch["fit_scope"] = "across"
            per_arch["condition"] = condition
            per_arch["K"] = K
            per_arch["outer_rep"] = outer_rep
            per_arch["direction"] = direction
            per_arch["outer_seed"] = outer_seed
            per_arch["n_ranking_subjects"] = len(ranking_subjects)
            new_ranks.append(per_arch)
            heldout = evaluate_topm_on_subjects(evaluation_subjects, at, ranked, top_m_values=TOP_M,
                                                nfolds=NFOLDS, nreps=NREPS, seed=evaluation_seed)
            heldout["analysis_type"] = at
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
        (d / "runs.partial.csv").rename(d / "runs.csv")
    return f"{at} seed{s} K={K} {condition} reps {todo[0]}-{todo[-1]} ({len(todo)})"


def main():
    run, mode, ph, idx, workers = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    base = base_dir(run, mode)
    at, s, ks = (base / f"tasks_{ph}.txt").read_text().splitlines()[idx - 1].split()
    s, ks = int(s), [int(k) for k in ks.split(":")]
    units = [(K, c, list(range(r, min(r + 5, NREP)))) for K in sorted(ks, reverse=True) for c in CONDS
             for r in range(0, NREP, 5)]
    log(f"{ph} task {idx}: {at} seed {s}, K {ks[0]}-{ks[-1]} ({len(ks)} K), {len(units)} units, {workers} workers")
    failed = 0
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("spawn"), initializer=init, initargs=(str(run),)) as ex:
        futs = {ex.submit(unit, str(run), str(base), at, s, K, c, reps): (K, c, reps) for K, c, reps in units}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                log(f"[{i}/{len(units)}] {f.result()}")
            except Exception:
                failed += 1
                log(f"[{i}/{len(units)}] FAILED {futs[f]}\n{traceback.format_exc()}")
    log(f"done, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
