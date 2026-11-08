#!/usr/bin/env python3
"""Outer-seed table, preflight and task lists.   python hk_setup.py RUN_DIR smoke|full GROUP_UNITS
Writes BASE/configs/heldout_outer_seeds_everyk.csv (verified identical if it already exists) and
BASE/tasks_<phase>.txt, one line per task: "<analysis> <seed> <K>:<K>:...".
Skips splits that are already finished (BASE/parts/...) and, for spatial, the K already in each
seed's existing held-out results (the 30 paper K), which the merge reuses. Exit 1 = preflight failed."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hk_common import CONDS, NREP, PHASES, base_dir, fit_file, old_complete_ks, outer_seed_table, phase, split_dir  # noqa

RUN, MODE, UNITS = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
BASE = base_dir(RUN, MODE)
(BASE / "configs").mkdir(parents=True, exist_ok=True)
t = outer_seed_table(RUN)
f = BASE / "configs" / "heldout_outer_seeds_everyk.csv"
if f.exists():
    old = pd.read_csv(f)
    if not old[["condition", "K", "outer_rep", "outer_seed"]].equals(t[["condition", "K", "outer_rep", "outer_seed"]]):
        sys.exit(f"{f} differs from the regenerated table; not touching it")
else:
    t.to_csv(f, index=False)
print(f"outer seeds: {len(t)} rows ({(t.source == 'paper_frozen').sum()} frozen paper rows kept)")

missing, summary = [], []
for ph in PHASES:
    at, seeds, ks = phase(MODE, ph)
    lines, nsplit = [], 0
    for s in seeds:
        reuse = old_complete_ks(RUN, s) if at == "spatial" else set()
        todo = []
        for K in ks:
            if not fit_file(RUN, at, s, K).exists():
                missing.append(str(fit_file(RUN, at, s, K)))
                continue
            if K in reuse:
                continue
            n = sum(not (split_dir(BASE, at, s, K, c, r) / "runs.csv").exists() for c in CONDS for r in range(NREP))
            if n:
                todo.append((K, n))
        # group consecutive K into array elements of ~UNITS cost; cost of a split ~ 2K single-archetype decodes + 2 full decodes (~20 single ones)
        group, cost = [], 0
        for K, n in todo:
            group.append(K); cost += n * (K + 20); nsplit += n
            if cost >= UNITS:
                lines.append(f"{at} {s} " + ":".join(map(str, group))); group, cost = [], 0
        if group:
            lines.append(f"{at} {s} " + ":".join(map(str, group)))
        summary.append(f"{ph} seed {s}: {len(todo)} K to run, {len(reuse & set(ks))} K reused")
    (BASE / f"tasks_{ph}.txt").write_text("".join(x + "\n" for x in lines))
    (BASE / f"tasks_{ph}.env").write_text(f"NTASK={len(lines)}\nNSPLIT={nsplit}\n")
    print(f"{ph} ({at}, seeds {seeds[0]}-{seeds[-1]}): {len(lines)} array tasks, {nsplit} splits to run")
for x in summary:
    print("  " + x)
if missing:
    print(f"PREFLIGHT FAILED: {len(missing)} fits missing, e.g.\n  " + "\n  ".join(missing[:10]))
    sys.exit(1)
