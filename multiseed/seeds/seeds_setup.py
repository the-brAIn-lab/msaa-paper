#!/usr/bin/env python3
"""Make the seed directories for a stage and write its task lists. Safe to rerun: existing
directories keep their outputs, and the task lists only hold what is still missing.

  python seeds_setup.py STAGE        (full | smoke)

Each seed directory gets a copy of src/, helpers/, rerun/*.py and configs/profiles.json (with the
"full" profile's K grid cut down to this stage), and a link to data/input. Copies rather than
links, because the seqnum code finds its root from its own file location.
Writes tasks/<STAGE>_picks.txt ("seed at scope cond chunk nchunks"), tasks/<STAGE>_fit.txt ("seed at scope K cond")
and tasks/<STAGE>_expected.json, and prints the counts.
"""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from seeds_common import (CONDITIONS, ROOT, fit_name, fits_of, seed_dir, stage_jobs, stage_Ks,  # noqa: E402
                          stage_seeds)

HERE = Path(__file__).resolve().parent
PICKS_K_PER_TASK = 35          # K values per picks task (interleaved, so each task gets small and large K)


def make_seed_dir(stage, seed):
    d = seed_dir(stage, seed)
    (d / "rerun").mkdir(parents=True, exist_ok=True)
    for sub in ["src", "helpers"]:
        if (ROOT / sub).exists() and not (d / sub).exists():
            shutil.copytree(ROOT / sub, d / sub, ignore=shutil.ignore_patterns("__pycache__"))
    # the driver lives in scripts/ in the repo; each seed dir keeps it under rerun/, as on Hellgate
    shutil.copy2(ROOT / "scripts" / "run_parallel.py", d / "rerun" / "run_parallel.py")
    (d / "data").mkdir(exist_ok=True)
    if not (d / "data" / "input").exists():
        (d / "data" / "input").symlink_to(ROOT / "data" / "input")
    (d / "configs").mkdir(exist_ok=True)
    prof = json.loads((ROOT / "configs" / "profiles.json").read_text())
    prof["full"]["jobs"] = {job: stage_Ks(stage, job.split("_")[0]) for job in stage_jobs()}
    (d / "configs" / "profiles.json").write_text(json.dumps(prof, indent=1))
    (d / "picks_parts").mkdir(exist_ok=True)
    (d / "slim").mkdir(exist_ok=True)
    (d / "msaa_flexible_outputs_npz").mkdir(exist_ok=True)
    (d / "SEED").write_text(f"{seed}\n")
    return d


def main():
    stage = sys.argv[1]
    seeds = stage_seeds(stage)
    for s in seeds:
        make_seed_dir(stage, s)

    picks, fits = [], []
    for s in seeds:
        for job in stage_jobs():
            at, sc = job.split("_")
            nch = -(-len(stage_Ks(stage, at)) // PICKS_K_PER_TASK)
            for cond in ([None] if sc == "across" else CONDITIONS):
                for c in range(nch):
                    part = seed_dir(stage, s) / "picks_parts" / f"{at}_{sc}_{cond or 'acrossCond'}_{c}of{nch}.csv"
                    if not part.exists():
                        picks.append(f"{s} {at} {sc} {cond or '-'} {c} {nch}\n")
    all_fits = fits_of(stage)
    for s, at, sc, K, cond in all_fits:
        if not (seed_dir(stage, s) / "slim" / fit_name(at, sc, K, cond)).exists():
            fits.append((s, at, sc, K, cond))
    # temporal first, then largest K first (long fits start early); seeds interleaved
    fits.sort(key=lambda t: (t[1] != "temporal", -t[3], t[2], t[4] or "", t[0]))
    tasks = HERE / "tasks"
    tasks.mkdir(exist_ok=True)
    (tasks / f"{stage}_picks.txt").write_text("".join(picks))
    (tasks / f"{stage}_fit.txt").write_text("".join(f"{s} {at} {sc} {K} {c or '-'}\n" for s, at, sc, K, c in fits))
    n_topm = sum(1 for f in all_fits if f[2] == "across")
    exp = dict(stage=stage, seeds=seeds, jobs=stage_jobs(), fits=len(all_fits), topm=n_topm,
               picks=len(seeds) * sum((1 if j.endswith("across") else 3) * -(-len(stage_Ks(stage, j.split("_")[0])) // PICKS_K_PER_TASK)
                                      for j in stage_jobs()))
    (tasks / f"{stage}_expected.json").write_text(json.dumps(exp))
    print(f"{stage}: seeds {seeds}; picks to do {len(picks)}/{exp['picks']}, fits to do {len(fits)}/{exp['fits']}, "
          f"top-m fits {n_topm}")


if __name__ == "__main__":
    main()
