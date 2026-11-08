#!/usr/bin/env python3
"""Merge one seed's held-out split shards into SEED_DIR/msaa_heldout_topm_outputs_spatial_across,
exactly as the held-out part of rerun/run_parallel.py merge() (decode tables are not touched)."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ROOT / "src"))
from msaa_pipeline.heldout import write_heldout_summaries  # noqa: E402
from msaa_pipeline.profiles import load_profile  # noqa: E402

prof = load_profile(ROOT, "full")
splits = sorted((ROOT / "rerun" / "shards" / "heldout_split").glob("K*_r*"))
splits = [s for s in splits if (s / "runs.csv").exists()]
runs = pd.concat([pd.read_csv(s / "runs.csv") for s in splits], ignore_index=True)
ranks = pd.concat([pd.read_csv(s / "ranks.csv") for s in splits], ignore_index=True)
key = ["analysis_type", "fit_scope", "condition", "K", "outer_rep", "direction"]
runs = runs[runs.groupby(key[:5])["direction"].transform("nunique") == 2]
ranks = ranks[ranks.groupby(key[:5])["direction"].transform("nunique") == 2]
runs = runs.drop_duplicates(key + ["top_m"], keep="last").sort_values(["condition", "K", "outer_rep", "direction", "top_m"])
ranks = ranks.drop_duplicates(key + ["archetype"], keep="last").sort_values(["condition", "K", "outer_rep", "direction", "archetype"])
out = ROOT / "msaa_heldout_topm_outputs_spatial_across"
out.mkdir(exist_ok=True)
runs.to_csv(out / "heldout_topm_crossfit_runs.csv", index=False)
ranks.to_csv(out / "heldout_per_archetype_ranking_runs.csv", index=False)
write_heldout_summaries(runs.reset_index(drop=True), ranks.reset_index(drop=True), out, prof["heldout_outer_reps"])
print(f"merged {len(splits)} splits, {len(runs)} run rows -> {out}", flush=True)
