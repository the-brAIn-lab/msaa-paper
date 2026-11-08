#!/usr/bin/env python3
"""Run the multi-seed stability analysis on one machine, without a cluster.

This is the same chain that produced the 10-seed results in the paper (it was run on an HPC
cluster as job arrays); here every step runs in order on local cores. Each step is resumable:
finished fits, decodes and held-out splits are skipped when you rerun.

  seed0      reference fits (seed 0) with the frozen FurthestSum picks in
             configs/furthestsum_init.csv.gz, then decode, held-out top-m (30 paper K) and merge
             (scripts/run_parallel.py fit / decode / heldout / merge in the repo root)
  picks      FurthestSum starting picks drawn with seeds 1..9           (seeds/picks_task.py, picks_merge.py)
  fits       one fit per seed x analysis x scope x K [x condition]      (seeds/seed_fit.py)
  decode     decode each seed's fits                                     (run_parallel.py decode + merge per seed dir)
  topm       same-data top-m for each seed's across fits                 (seeds/seed_topm.py)
  merge      cross-seed tables: matching to seed 0, per-K mean/SD, ...   (seeds/seeds_merge.py)
  heldout30  held-out top-m at the 30 paper K for every seed              (run_parallel.py heldout, heldout_seeds/)
  heldoutK   held-out top-m at every K for every seed                     (heldout_everyk/)

Outputs land in the repo root, in the same layout as the original run:
  seeds/full/seedNN/   per-seed run roots      seeds/results/, seeds-results.tar.gz
  heldout-seeds-results/   heldout-everyk/ and heldout-everyk-<analysis>-<seed0|all>-results/

Usage (from anywhere):
  python analysis/08_multiseed_stability/run_local.py --workers 8                       # everything (weeks of CPU time)
  python analysis/08_multiseed_stability/run_local.py --ks 5,10 --seeds 1,2 --jobs temporal_across --workers 4   # small test
  python analysis/08_multiseed_stability/run_local.py --steps decode,topm,merge          # just some steps

The full run is very large (about 36,000 fits plus held-out splits at every K; it took several
days on ~100 cluster cores). --ks / --seeds / --jobs restrict it for testing.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STEPS = ["seed0", "picks", "fits", "decode", "topm", "merge", "heldout30", "heldoutK"]
JOBS = ["spatial_across", "spatial_within", "temporal_across", "temporal_within"]
HK_PHASES = ["s0", "sA", "t0", "tA"]


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def run(cmd, cwd=ROOT, env=None, check=True):
    log("$ " + " ".join(map(str, cmd)))
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env)
    if check and r.returncode != 0:
        raise SystemExit(f"step failed (exit {r.returncode}): {' '.join(map(str, cmd))}")
    return r.returncode


def run_many(cmds, workers, env, what):
    """Run independent commands in parallel; fail at the end if any failed."""
    if not cmds:
        log(f"{what}: nothing to do")
        return
    log(f"{what}: {len(cmds)} tasks on {workers} workers")
    failed, done = [], 0
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(subprocess.run, [str(c) for c in cmd], cwd=HERE / "seeds", env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True): cmd for cmd in cmds}
        for f in as_completed(futs):
            done += 1
            r = f.result()
            if r.returncode != 0:
                failed.append(futs[f])
                log(f"[{done}/{len(cmds)}] FAILED {' '.join(map(str, futs[f]))}\n{r.stdout[-3000:]}")
            elif done % max(1, len(cmds) // 20) == 0 or done == len(cmds):
                log(f"[{done}/{len(cmds)}] {what}")
    if failed:
        raise SystemExit(f"{what}: {len(failed)} of {len(cmds)} tasks failed; rerun to retry them")


def n_lines(p):
    p = Path(p)
    return len(p.read_text().splitlines()) if p.exists() else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", default=",".join(STEPS), help=f"comma list of steps (default: all): {','.join(STEPS)}")
    ap.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9", help="extra fit seeds (seed 0 is the reference)")
    ap.add_argument("--ks", default=None, help="comma list restricting K (default: every K)")
    ap.add_argument("--jobs", default=",".join(JOBS), help="comma list of analysis_scope jobs")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--threads", type=int, default=1, help="BLAS threads per fit")
    a = ap.parse_args()

    steps = a.steps.split(",")
    bad = [s for s in steps if s not in STEPS]
    if bad:
        raise SystemExit(f"unknown steps {bad}; choose from {STEPS}")
    seeds = [int(s) for s in a.seeds.split(",")]
    py = sys.executable
    env = dict(os.environ, SEEDS=" ".join(map(str, seeds)), SEED_JOBS=a.jobs, MSAA_THREADS=str(a.threads),
               MSAA_TIE_TOL=os.environ.get("MSAA_TIE_TOL", "1e-08"), PYTHONUNBUFFERED="1")
    one = dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if a.ks:
        env.update(SEED_KS=a.ks, HK_KS=a.ks, HS_KS=a.ks)
    env["HK_SEEDS"] = ",".join(map(str, [0] + seeds))
    types = sorted({j.split("_")[0] for j in a.jobs.split(",") if j.endswith("_across")})
    env["HK_TYPES"] = ",".join(types)
    ks = ["--ks", a.ks] if a.ks else []
    rp = ROOT / "scripts" / "run_parallel.py"
    if not (ROOT / "data" / "input" / "pieman_data.mat").exists():
        raise SystemExit("data/input/pieman_data.mat is missing (see data/input/README.md)")

    if "seed0" in steps:
        run([py, rp, "fit", "--jobs", a.jobs, *ks, "--workers", a.workers, "--threads", a.threads], env=env)
        run([py, rp, "decode", "--jobs", a.jobs, *ks, "--workers", a.workers], env=env)
        if "spatial_across" in a.jobs.split(","):
            run([py, rp, "heldout", *ks, "--workers", a.workers], env=env)
        run([py, rp, "merge"], env=env)

    if "picks" in steps:
        run([py, "seeds_setup.py", "full"], cwd=HERE / "seeds", env=env)
        n = n_lines(HERE / "seeds" / "tasks" / "full_picks.txt")
        run_many([[py, "picks_task.py", "full", i] for i in range(1, n + 1)], a.workers, {**env, **one}, "picks")
        run([py, "picks_merge.py", "full"], cwd=HERE / "seeds", env=env)

    if "fits" in steps:
        run([py, "seeds_setup.py", "full"], cwd=HERE / "seeds", env=env)   # refresh the fit list
        n = n_lines(HERE / "seeds" / "tasks" / "full_fit.txt")
        th = dict(OMP_NUM_THREADS=str(a.threads), OPENBLAS_NUM_THREADS=str(a.threads), MKL_NUM_THREADS=str(a.threads))
        run_many([[py, "seed_fit.py", "full", i] for i in range(1, n + 1)], a.workers, {**env, **th}, "fits")

    if "decode" in steps:
        for s in seeds:
            d = ROOT / "seeds" / "full" / f"seed{s:02d}"
            run([py, "rerun/run_parallel.py", "decode", "--workers", a.workers], cwd=d, env=env, check=False)
            run([py, "rerun/run_parallel.py", "merge"], cwd=d, env=env, check=False)

    if "topm" in steps:
        ntask = max(1, a.workers * 4)
        run_many([[py, "seed_topm.py", "full", i, ntask] for i in range(1, ntask + 1)], a.workers, {**env, **one}, "top-m")

    if "merge" in steps:
        run([py, "seeds_merge.py", "--workers", a.workers], cwd=HERE / "seeds", env={**env, **one}, check=False)

    if "heldout30" in steps and "spatial_across" in a.jobs.split(","):
        for s in seeds:
            d = ROOT / "seeds" / "full" / f"seed{s:02d}"
            shutil.copy2(ROOT / "configs" / "heldout_outer_seeds.csv", d / "configs" / "heldout_outer_seeds.csv")
            run([py, "rerun/run_parallel.py", "heldout", *ks, "--workers", a.workers], cwd=d, env=env)
            run([py, HERE / "heldout_seeds" / "hseed_merge.py", d], env=env)
        run([py, HERE / "heldout_seeds" / "hsummary.py", ROOT, "full", min(seeds), max(seeds)], env=env, check=False)

    if "heldoutK" in steps:
        kit = HERE / "heldout_everyk"
        run([py, kit / "hk_setup.py", ROOT, "full", 150000], env=env)
        for ph in [p for p in HK_PHASES if {"s": "spatial", "t": "temporal"}[p[0]] in types]:
            n = n_lines(ROOT / "heldout-everyk" / f"tasks_{ph}.txt")
            for i in range(1, n + 1):
                run([py, kit / "hk_task.py", ROOT, "full", ph, i, a.workers], env=env)
            run([py, kit / "hk_merge.py", ROOT, "full", ph, a.workers], env=env, check=False)
    log("done")


if __name__ == "__main__":
    main()
