#!/usr/bin/env python3
"""Merge and validate.   python hm_merge.py RUN_DIR smoke|full
Writes RUN_DIR/heldout-highm[-landscape][-smoke]-results/:
  heldout_highm_runs_all_seeds.csv.gz  every seed x analysis x condition x K x split x direction x top_m
                                       (the new m plus m=10), long form, with fit_seed
  heldout_highm_per_K.csv              directions averaged within split -> splits within seed -> seeds:
                                       n_seeds, mean, sd, min, max, seed0 per analysis x condition x K x top_m
  validation.csv / validation.txt      m=10 and the full model vs the every-K run on the same cells
                                       (must be identical: same rankings, splits, decoder seeds, code)
  summary.txt                          first line COMPLETE / INCOMPLETE and VALIDATION PASS / FAIL
Exit 1 if incomplete or validation fails."""
import ast
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hm_common import CONDS, OLD_SUB, base_dir, seed_root  # noqa: E402

RUN, MODE = Path(sys.argv[1]), sys.argv[2]
BASE = base_dir(RUN, MODE); OUT = RUN / (BASE.name + "-results"); OUT.mkdir(parents=True, exist_ok=True)
NREP = int((BASE / "nrep.txt").read_text())
KEY = ["condition", "K", "outer_rep", "direction", "top_m"]

expected, missing, frames = 0, [], []
for f in sorted(BASE.glob("tasks_*.txt")):
    for line in f.read_text().strip().splitlines():
        at, s, K = line.split(); s, K = int(s), int(K)
        for c in CONDS:
            for r in range(NREP):
                expected += 1
                p = BASE / "parts" / at / f"seed{s:02d}" / f"K{K}_{c}_r{r}" / "runs.csv"
                if p.exists():
                    frames.append(pd.read_csv(p))
                else:
                    missing.append(str(p))
runs = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
complete = not missing and len(runs) > 0
if len(runs):
    runs.to_csv(OUT / "heldout_highm_runs_all_seeds.csv.gz", index=False)

# ---- validation against the every-K run (m=10 and full model on the same cells)
val = []
for (at, s), sub in (runs.groupby(["analysis_type", "fit_seed"]) if len(runs) else []):
    ref = RUN / "heldout-everyk" / "merged" / at / f"seed{int(s):02d}" / "heldout_topm_crossfit_runs.csv.gz"
    if not ref.exists():
        ref = seed_root(RUN, int(s)) / OLD_SUB / "heldout_topm_crossfit_runs.csv"      # spatial 30-K (smoke)
    if not ref.exists():
        print(f"{at} seed {s}: no every-K reference found"); continue
    old = pd.read_csv(ref)
    if "analysis_type" in old:
        old = old[old.analysis_type == at]
    m = sub.merge(old[KEY + ["mean", "full_mean_same_eval", "selected_archetypes"]], on=KEY, suffixes=("", "_ref"))
    m["fit_seed"] = s
    val.append(m)
ok = False
if val:
    v = pd.concat(val, ignore_index=True)
    v["d_mean"] = (v["mean"] - v["mean_ref"]).abs()
    v["d_full"] = (v["full_mean_same_eval"] - v["full_mean_same_eval_ref"]).abs()
    v["same_set"] = [ast.literal_eval(str(a)) == ast.literal_eval(str(b))
                     for a, b in zip(v.selected_archetypes, v.selected_archetypes_ref)]
    v.to_csv(OUT / "validation.csv", index=False)
    ok = bool(v.d_mean.max() <= 1e-9 and v.d_full.max() <= 1e-9 and v.same_set.all())
    vtxt = (f"VALIDATION {'PASS' if ok else 'FAIL'}: {len(v)} overlapping cells (m=10 + full model), "
            f"max |diff| top-m {v.d_mean.max():.2e}, full model {v.d_full.max():.2e}, "
            f"selected sets identical {v.same_set.mean():.0%}")
else:
    vtxt = "VALIDATION FAIL: no overlapping cells found"
(OUT / "validation.txt").write_text(vtxt + "\n")

# ---- per-K summary: directions within split -> splits within seed -> seeds
if len(runs):
    g = ["analysis_type", "condition", "K", "top_m"]
    split = runs.groupby(g + ["fit_seed", "outer_rep"], as_index=False)[["mean", "full_mean_same_eval"]].mean()
    seed = split.groupby(g + ["fit_seed"], as_index=False).agg(mean=("mean", "mean"), full=("full_mean_same_eval", "mean"),
                                                               n_splits=("outer_rep", "nunique"))
    per = seed.groupby(g).agg(n_seeds=("fit_seed", "nunique"), mean=("mean", "mean"), sd=("mean", "std"),
                              min=("mean", "min"), max=("mean", "max"), full_mean=("full", "mean")).reset_index()
    s0 = seed[seed.fit_seed == 0].set_index(g)["mean"].rename("seed0")
    per = per.join(s0, on=g)
    per.to_csv(OUT / "heldout_highm_per_K.csv", index=False)

head = f"{'COMPLETE' if complete else 'INCOMPLETE'} | {expected - len(missing)}/{expected} splits | {vtxt}"
lines = [head]
if missing:
    lines += [f"missing {len(missing)}, first few:"] + ["  " + m for m in missing[:8]]
(OUT / "summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(0 if (complete and ok) else 1)
