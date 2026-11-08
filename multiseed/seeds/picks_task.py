#!/usr/bin/env python3
"""FurthestSum picks for one (seed, analysis type, scope, condition) at every nch-th K of a stage: line N
of tasks/<STAGE>_picks.txt. Same call as rerun/freeze_init.py / everyk/freeze_everyk.py, with
seed=<seed> instead of 42. Writes <seed dir>/picks_parts/<at>_<scope>_<cond>_<c>of<nch>.csv.

  python picks_task.py STAGE N
"""
import os, sys, time
for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ.setdefault(v, "1")
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from seeds_common import CONDITIONS, seed_dir, stage_Ks  # noqa: E402

HERE = Path(__file__).resolve().parent


def main():
    stage, n = sys.argv[1], int(sys.argv[2])
    s, at, sc, cond, c, nch = (HERE / "tasks" / f"{stage}_picks.txt").read_text().splitlines()[n - 1].split()
    s, cond, c, nch = int(s), (None if cond == "-" else cond), int(c), int(nch)
    d = seed_dir(stage, s)
    part = d / "picks_parts" / f"{at}_{sc}_{cond or 'acrossCond'}_{c}of{nch}.csv"
    if part.exists():
        print("skip", part.name); return
    sys.path.insert(0, str(d / "src"))
    import pandas as pd
    from msaa_pipeline.data import load_inputs
    from msaa_pipeline.fit import furthest_sum_picks
    D, _, _ = load_inputs(d, trim_rest_by=100, max_subjects_per_condition=None)
    subjects = [x for c in CONDITIONS for x in D[c]] if sc == "across" else D[cond]
    rows, t = [], time.time()
    for K in sorted(stage_Ks(stage, at), reverse=True)[c::nch]:
        idx = furthest_sum_picks(subjects, K, at, d, seed=s)
        rows.append(pd.DataFrame(dict(analysis_type=at, fit_scope=sc, condition=cond or "acrossCond", K=K,
                                      archetype=range(K), init_index=idx)))
    tmp = part.with_suffix(".tmp")
    pd.concat(rows, ignore_index=True).to_csv(tmp, index=False); tmp.rename(part)
    print(f"done seed {s} {at} {sc} {cond or 'acrossCond'} chunk {c}/{nch}: {len(rows)} K in {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    main()
