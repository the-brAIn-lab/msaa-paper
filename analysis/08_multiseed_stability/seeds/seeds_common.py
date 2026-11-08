"""Shared pieces of the multi-seed run: where things live, which K each stage covers, fit names.

Layout ($RUN_DIR = the run root, by default the repo root):
  $RUN_DIR/msaa_flexible_outputs_npz/      seed 0 = the existing every-K seqnum fits (canonical)
  $RUN_DIR/seeds/<stage>/seed01 .. seed09/ one self-contained run root per seed and stage, a copy of
                                           the seqnum code with its own configs/furthestsum_init.csv.gz
                                           (FurthestSum picks drawn with that seed) and its own outputs
  $RUN_DIR/seeds/slim/seed00/              compact copies of the canonical fits, made by the merge
  $RUN_DIR/seeds/results/                  merged tables (seeds-results)
Each seed directory looks exactly like $RUN_DIR to the seqnum code, so fitting and decoding run
the unmodified run_parallel.py there.
"""
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]          # run root = repo root (this file is in analysis/08_multiseed_stability/seeds)
SEEDS_DIR = ROOT / "seeds"
CONDITIONS = ["intact", "word", "rest"]
GRID = {"spatial": range(2, 701), "temporal": range(2, 301)}
JOBS = ["spatial_across", "spatial_within", "temporal_across", "temporal_within"]
TOP_M_VALUES = [1, 2, 3, 5, 10]


def stage_Ks(stage, at):
    """K values a stage covers. full = every K (spatial 2-700, temporal 2-300); smoke = 5 and 10.
    SEED_KS (comma list, set by analysis/08_multiseed_stability/run_local.py --ks) restricts the full grid."""
    if stage == "full":
        keep = os.environ.get("SEED_KS")
        ks = list(GRID[at])
        return [k for k in ks if k in {int(x) for x in keep.split(",")}] if keep else ks
    if stage == "smoke":
        return [5, 10]
    raise SystemExit(f"unknown stage {stage!r} (full, smoke)")


def stage_seeds(stage):
    if stage == "smoke":
        return [1, 2]
    return [int(s) for s in os.environ.get("SEEDS", "1 2 3 4 5 6 7 8 9").split()]


def stage_jobs():
    return os.environ.get("SEED_JOBS", ",".join(JOBS)).split(",")


def seed_dir(stage, seed):
    return SEEDS_DIR / stage / f"seed{seed:02d}"


def fit_name(at, sc, K, cond):
    return f"{at}AA_across_acrossCond_K{K}.npz" if sc == "across" else f"{at}AA_within_{cond}_K{K}.npz"


def fits_of(stage):
    """[(seed, at, sc, K, cond-or-None)] for a stage."""
    out = []
    for seed in stage_seeds(stage):
        for job in stage_jobs():
            at, sc = job.split("_")
            for K in stage_Ks(stage, at):
                for cond in ([None] if sc == "across" else CONDITIONS):
                    out.append((seed, at, sc, K, cond))
    return out


def profile(root):
    return json.loads((Path(root) / "configs" / "profiles.json").read_text())["full"]


def all_seed_dirs():
    """Every seed directory of every stage that exists, as [(stage, seed, path)]."""
    out = []
    for st in sorted(SEEDS_DIR.glob("*")):
        if st.is_dir() and st.name in ("full", "smoke"):
            for d in sorted(st.glob("seed[0-9][0-9]")):
                out.append((st.name, int(d.name[4:]), d))
    return out


def write_slim(fit_path, out_path):
    """Compact copy of one fit, enough for matching and loss comparisons (a few MB at most):
    C (shared generator), subject-mean S, final cost, iterations, varexpl, 0.5*SSE, init_index."""
    import numpy as np
    z = np.load(fit_path, allow_pickle=True)
    subj = z["results_subj"].tolist()
    cost = np.asarray(z["cost_fun"].tolist(), float).ravel()
    d = dict(C=np.asarray(z["C"].tolist(), float).astype(np.float32),
             Smean=np.mean([np.asarray(r["S"], float) for r in subj], axis=0).astype(np.float32),
             cost_last=np.float64(cost[-1] if cost.size else np.nan), n_iter=np.int64(cost.size),
             varexpl=np.float64(z["varexpl"]) if "varexpl" in z.files else np.float64(np.nan),
             sse_half=np.float64(0.5 * sum(float(r["SSE"]) for r in subj)) if "SSE" in subj[0] else np.float64(np.nan),
             init_index=np.asarray(z["init_index"] if "init_index" in z.files else [-1]))
    out_path = Path(out_path)
    tmp = out_path.with_name(out_path.stem + ".tmp.npz")
    np.savez(tmp, **d)
    tmp.rename(out_path)
