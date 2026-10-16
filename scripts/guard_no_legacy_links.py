#!/usr/bin/env python3
"""Refuse write-producing full jobs while portable legacy links are active."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

manifest = Path("configs/legacy_outputs_manifest.csv")
if not manifest.is_file():
    sys.exit(0)

with manifest.open(newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
roots = sorted({row.get("root", "") for row in rows if row.get("root")})
active = [Path(root) for root in roots if Path(root).is_symlink()]
if active:
    print("Refusing to run a full write-producing analysis while validated legacy links are active:")
    for p in active:
        print(f"  {p} -> {p.readlink()}")
    print("Run `make deactivate-legacy` first. The archived bytes remain safe in reference_outputs/legacy_full/.")
    sys.exit(2)
