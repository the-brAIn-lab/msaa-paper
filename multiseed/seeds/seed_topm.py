#!/usr/bin/env python3
"""Same-data top-m for the seed fits (across scope, as topm-everyk did for seed 0), one array task.
Exactly topm-everyk/topm_task.py's procedure, run in each seed directory: rank archetypes per
condition by that seed's own single-archetype accuracy (its per_archetype_mean_accuracy.csv), then
decode the top m (1, 2, 3, 5, 10) with msaa_pipeline.decode (2 folds x 20 reps, seed 42).
One CSV per fit in <seed dir>/topm_parts/; existing parts are skipped (resumable).

  python seed_topm.py STAGE TASK_ID N_TASKS
"""
import os, sys, time
for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ.setdefault(v, "1")
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from seeds_common import CONDITIONS, TOP_M_VALUES, fit_name, fits_of, seed_dir  # noqa: E402


def ranking(pa, cond, K):
    sub = pa[(pa.condition.astype(str) == cond) & (pa.K == K)]
    m = sub.groupby("archetype")["mean_accuracy"].mean().reset_index()
    m = m.sort_values(["mean_accuracy", "archetype"], ascending=[False, True], kind="stable")
    return m.archetype.astype(int).tolist()


def main():
    stage, task, ntasks = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    jobs = sorted([f for f in fits_of(stage) if f[2] == "across"], key=lambda f: (-f[3], f[1], f[0]))
    mine = jobs[task - 1::ntasks]
    if not mine:
        return
    sys.path.insert(0, str(seed_dir(stage, mine[0][0]) / "src"))
    import numpy as np
    import pandas as pd
    from msaa_pipeline.decode import decode_subjects, summarize
    pa_cache = {}
    for s, at, sc, K, _ in mine:
        d = seed_dir(stage, s)
        part = d / "topm_parts" / f"{at}_{sc}_acrossCond_K{K}.csv"
        if part.exists():
            print("skip", s, part.name, flush=True); continue
        part.parent.mkdir(exist_ok=True)
        if (s, at) not in pa_cache:
            pa_cache[s, at] = pd.read_csv(d / f"msaa_condrank_decoding_outputs_{at}_{sc}" / "per_archetype_mean_accuracy.csv")
        pa = pa_cache[s, at]
        src = d / "msaa_flexible_outputs_npz" / fit_name(at, sc, K, None)
        z = np.load(src, allow_pickle=True)
        subs_all = z["results_subj"].tolist()
        labels = np.asarray(z["condition_labels_str"].tolist(), dtype=object)
        t, rows = time.time(), []
        for c in CONDITIONS:
            subjects = [subs_all[i] for i in np.where(labels == c)[0]]
            ranked = ranking(pa, c, K)
            if len(ranked) != K:
                raise RuntimeError(f"seed {s} {at} K{K} {c}: {len(ranked)} ranked archetypes, expected {K}")
            for m in TOP_M_VALUES:
                sel = ranked[:m]
                rows.append(dict(seed=s, analysis_type=at, fit_scope=sc, condition=c, K=K, top_m=m,
                                 selected_archetypes=str(sel), **summarize(decode_subjects(subjects, at, sel, 2, 20, 42))))
        tmp = part.with_suffix(".tmp")
        pd.DataFrame(rows).to_csv(tmp, index=False); tmp.rename(part)
        print(f"done seed {s} {at} K{K} {time.time() - t:.1f}s", flush=True)


if __name__ == "__main__":
    main()
