#!/usr/bin/env python3
"""Compute and freeze the FurthestSum initialisation picks for every fit.

Writes configs/furthestsum_init.csv.gz with one row per (analysis_type,
fit_scope, condition, K, archetype): init_index is the data row/column
(node for spatial, time point for temporal, 0-based) archetype j starts from.
Fits read these instead of recomputing FurthestSum, so the archetype
numbering is identical on every machine. Existing rows are kept unchanged;
only missing fits are added, so extending the K grid never renumbers old fits.

Usage (from the repo root):
  python scripts/freeze_init.py [--profile full] [--workers 8]
"""
from __future__ import annotations
import argparse, os, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "configs" / "furthestsum_init.csv.gz"
KEY = ["analysis_type", "fit_scope", "condition", "K"]


def _picks(at, sc, cond, K):
    sys.path.insert(0, str(ROOT / "src"))
    from msaa_pipeline.data import load_inputs
    from msaa_pipeline.fit import furthest_sum_picks
    data, _, _ = load_inputs(ROOT, trim_rest_by=100, max_subjects_per_condition=None)
    subjects = [x for c in ["intact", "word", "rest"] for x in data[c]] if sc == "across" else data[cond]
    idx = furthest_sum_picks(subjects, K, at, ROOT, seed=42)
    return pd.DataFrame(dict(analysis_type=at, fit_scope=sc, condition=cond or "acrossCond", K=K,
                             archetype=range(K), init_index=idx))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="full")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        os.environ[v] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from msaa_pipeline.profiles import load_profile
    prof = load_profile(ROOT, args.profile)
    old = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=KEY + ["archetype", "init_index"])
    have = set(map(tuple, old[KEY].astype(str).drop_duplicates().itertuples(index=False, name=None)))
    todo = []
    for job, ks in prof["jobs"].items():
        at, sc = job.split("_")
        for K in ks:
            for cond in ([None] if sc == "across" else ["intact", "word", "rest"]):
                if (at, sc, cond or "acrossCond", str(K)) not in have:
                    todo.append((at, sc, cond, int(K)))
    print(f"{len(todo)} fits need picks ({len(have)} already frozen)", flush=True)
    with ProcessPoolExecutor(args.workers) as ex:
        new = list(ex.map(_picks, *zip(*todo))) if todo else []
    df = pd.concat([old, *new], ignore_index=True).sort_values(KEY + ["archetype"]).reset_index(drop=True)
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {df.groupby(KEY).ngroups} fits, {len(df)} rows")


if __name__ == "__main__":
    main()
