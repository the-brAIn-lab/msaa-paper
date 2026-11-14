#!/usr/bin/env python3
"""Write the task lists and check every input exists.   python hm_setup.py RUN_DIR smoke|full
One task = one (analysis, seed, K): all 3 conditions x 20 outer splits x 2 directions.
smoke = seed 0, spatial K=35 + temporal K=15 (m 7, 10, 15), 2 outer splits.
HM_TYPES=spatial (or temporal) restricts to one analysis. Exit 1 if anything is missing."""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hm_common import (CONDS, GRID, NREP, OLD_SUB, base_dir, everyk_split_dir, m_values, pair_ks, seed_root,  # noqa
                       tie_tol_ok)

run, mode = Path(sys.argv[1]), sys.argv[2]
types = os.environ.get("HM_TYPES", "spatial,temporal").split(",")
if mode == "smoke":
    plan = ({"spatial": ([0], [50]), "temporal": ([0], [])} if GRID == "landscape"
            else {"spatial": ([0], [35]), "temporal": ([0], [15])}); nrep = 2
else:
    plan = {at: (list(range(10)), [K for K in pair_ks(at) if m_values(at, K)]) for at in ("spatial", "temporal")}
    nrep = NREP
base = base_dir(run, mode); base.mkdir(parents=True, exist_ok=True)
if not tie_tol_ok():
    print("MSAA_TIE_TOL is not set (export MSAA_TIE_TOL=1e-08): refusing, results would not match the every-K run")
    sys.exit(1)
missing, total = [], 0
for at in types:
    seeds, ks = plan[at]
    if not ks:
        (base / f"tasks_{at}.txt").write_text("")
        continue
    lines = []
    for s in seeds:
        old = seed_root(run, s) / OLD_SUB / "heldout_per_archetype_ranking_runs.csv"
        for K in ks:
            fit = seed_root(run, s) / "msaa_flexible_outputs_npz" / f"{at}AA_across_acrossCond_K{K}.npz"
            if not fit.exists():
                missing.append(str(fit))
            for c in CONDS:
                for r in range(nrep):
                    if not (everyk_split_dir(run, at, s, K, c, r) / "ranks.csv").exists() and not (at == "spatial" and old.exists()):
                        missing.append(str(everyk_split_dir(run, at, s, K, c, r) / "ranks.csv"))
            lines.append(f"{at} {s} {K}")
            total += sum(len(m_values(at, K)) for _ in CONDS) * nrep * 2
    (base / f"tasks_{at}.txt").write_text("\n".join(lines) + "\n")
    print(f"{at}: {len(lines)} tasks (seeds {seeds[0]}-{seeds[-1]}, K {ks})")
(base / "nrep.txt").write_text(f"{nrep}\n")
print(f"{total} top-m decodes (+1 full-model decode per split x direction)")
if missing:
    print(f"MISSING {len(missing)} inputs, first few:"); [print("  ", m) for m in missing[:10]]
    sys.exit(1)
print("preflight OK")
