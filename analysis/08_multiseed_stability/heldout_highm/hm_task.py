#!/usr/bin/env python3
"""One task: every unfinished outer split of one (analysis, seed, K).
   python hm_task.py RUN_DIR smoke|full ANALYSIS INDEX(1-based) WORKERS
For each split x direction: rebuild the ranking from the every-K run's saved per-archetype rows
(order: mean_accuracy desc, rank_mean desc, archetype asc = msaa_pipeline.heldout.rank_archetypes_on_subjects),
redo the same participant split from its outer_seed, and score the top-m on the evaluation half with
msaa_pipeline.heldout.evaluate_topm_on_subjects and the every-K run's evaluation seed.
Output: BASE/parts/<analysis>/seedNN/K<K>_<cond>_r<rep>/runs.csv (written last; reruns skip finished splits)."""
import os
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[_v] = "1"
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import multiprocessing as mp
from pathlib import Path

KIT = Path(__file__).resolve().parent
sys.path.insert(0, str(KIT))
from hm_common import (BASE_SEED, CONDS, NFOLDS, NREPS, OLD_SUB, base_dir, everyk_split_dir, m_values,  # noqa
                       out_dir, seed_root, tie_tol_ok)


def log(msg):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def init(run):
    for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        os.environ[v] = "1"
    sys.path.insert(0, str(Path(run) / "src"))
    sys.path.insert(0, str(KIT))


@lru_cache(maxsize=2)
def old_ranks(path, K):
    import pandas as pd
    x = pd.read_csv(path)
    return x[x.K == K]


def ranking_rows(run, at, s, K, c, r):
    import pandas as pd
    f = everyk_split_dir(run, at, s, K, c, r) / "ranks.csv"
    if f.exists():
        return pd.read_csv(f), str(f)
    old = seed_root(run, s) / OLD_SUB / "heldout_per_archetype_ranking_runs.csv"
    if at == "spatial" and old.exists():
        x = old_ranks(str(old), K)
        x = x[(x.condition == c) & (x.outer_rep == r)]
        if len(x):
            return x, str(old)
    raise FileNotFoundError(f"no saved ranking for {at} seed{s} K={K} {c} r{r}")


def unit(run, base, at, s, K, condition, reps):
    import numpy as np
    import pandas as pd
    from msaa_pipeline.decode import get_subjects
    from msaa_pipeline.heldout import _split_outer_subject_indices, _subset_subjects, evaluate_topm_on_subjects
    todo = [r for r in reps if not (out_dir(base, at, s, K, condition, r) / "runs.csv").exists()]
    if not todo:
        return f"skip {at} seed{s} K={K} {condition}"
    ms = m_values(at, K)
    subjects, source = get_subjects(seed_root(Path(run), s), at, "across", K, condition)
    for outer_rep in todo:
        ranks, rsrc = ranking_rows(Path(run), at, s, K, condition, outer_rep)
        seeds = ranks.outer_seed.unique()
        assert len(seeds) == 1, (at, s, K, condition, outer_rep, seeds)
        outer_seed = int(seeds[0])
        A_idx, B_idx = _split_outer_subject_indices(len(subjects), np.random.default_rng(outer_seed))
        rows = []
        for direction, ranking_idx, evaluation_idx in [("A_to_B", A_idx, B_idx), ("B_to_A", B_idx, A_idx)]:
            pa = ranks[ranks.direction == direction]
            assert len(pa) == K and pa.archetype.nunique() == K, (at, s, K, condition, outer_rep, direction, len(pa))
            assert int(pa.n_ranking_subjects.iloc[0]) == len(ranking_idx)
            ranked = (pa.sort_values(["mean_accuracy", "rank_mean", "archetype"], ascending=[False, False, True])
                      ["archetype"].astype(int).tolist())
            direction_offset = 0 if direction == "A_to_B" else 1
            evaluation_seed = BASE_SEED + 200000 * outer_rep + 1000 * K + direction_offset
            ev = _subset_subjects(subjects, evaluation_idx)
            h = evaluate_topm_on_subjects(ev, at, ranked, top_m_values=ms, nfolds=NFOLDS, nreps=NREPS,
                                          seed=evaluation_seed)
            h = h.assign(analysis_type=at, fit_scope="across", condition=condition, K=K, outer_rep=outer_rep,
                         direction=direction, outer_seed=outer_seed, fit_seed=s,
                         n_ranking_subjects=len(ranking_idx), n_evaluation_subjects=len(ev),
                         fit_source_basename=source.name, ranking_source=Path(rsrc).name)
            rows.append(h)
        d = out_dir(base, at, s, K, condition, outer_rep); d.mkdir(parents=True, exist_ok=True)
        pd.concat(rows, ignore_index=True).to_csv(d / "runs.partial.csv", index=False)
        (d / "runs.partial.csv").rename(d / "runs.csv")
    return f"{at} seed{s} K={K} {condition} reps {todo[0]}-{todo[-1]} ({len(todo)}), m {ms}"


def main():
    run, mode, at, idx, workers = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    if not tie_tol_ok():
        sys.exit("MSAA_TIE_TOL is not set (export MSAA_TIE_TOL=1e-08): refusing")
    base = base_dir(run, mode)
    _, s, K = (base / f"tasks_{at}.txt").read_text().split("\n")[idx - 1].split()
    s, K = int(s), int(K)
    nrep = int((base / "nrep.txt").read_text())
    step = 5 if nrep >= 5 else nrep
    units = [(c, list(range(r, min(r + step, nrep)))) for c in CONDS for r in range(0, nrep, step)]
    log(f"task {idx}: {at} seed {s} K={K}, m {m_values(at, K)}, {len(units)} units, {workers} workers, "
        f"MSAA_TIE_TOL={os.environ['MSAA_TIE_TOL']}")
    failed = 0
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("spawn"), initializer=init, initargs=(str(run),)) as ex:
        futs = {ex.submit(unit, str(run), str(base), at, s, K, c, reps): (c, reps) for c, reps in units}
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
