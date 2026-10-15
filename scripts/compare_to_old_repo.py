#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from msaa_pipeline.compare import compare_outputs

ap = argparse.ArgumentParser()
ap.add_argument("old_repo", help="Path to the original MS_AA repository")
ap.add_argument("--profile", default="mini-paper", choices=["smoke", "mini-paper", "full"])
ap.add_argument("--new-root", default=None)
args = ap.parse_args()
new_root = Path(args.new_root).resolve() if args.new_root else REPO / "outputs" / args.profile.replace("-", "_")
res = compare_outputs(new_root, Path(args.old_repo).resolve())
out = new_root / "comparison_to_original.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
print("Saved:", out)
