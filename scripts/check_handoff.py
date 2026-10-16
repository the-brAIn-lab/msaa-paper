#!/usr/bin/env python3
"""Check that a copied project folder contains the data needed for student handoff."""
from __future__ import annotations

import csv
from pathlib import Path

required_inputs = [
    Path("data/input/pieman_data.mat"),
    Path("data/input/pieman_posterior_K700.mat"),
]
manifest = Path("configs/legacy_outputs_manifest.csv")
store = Path("reference_outputs/legacy_full")

problems = []
for p in required_inputs:
    if not p.is_file():
        problems.append(f"missing input: {p}")
if not manifest.is_file():
    problems.append(f"missing manifest: {manifest}")
if not store.is_dir():
    problems.append(f"missing legacy store: {store}")

if not problems and manifest.is_file():
    with manifest.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    missing = [row["relative_path"] for row in rows if not (store / row["relative_path"]).is_file()]
    if missing:
        problems.append(f"{len(missing)} manifest file(s) missing from legacy store")

if problems:
    print("Handoff check FAILED:")
    for p in problems:
        print(f"  - {p}")
    raise SystemExit(2)

print("Handoff check passed.")
print("This folder contains both HTFA MAT inputs and the portable validated legacy-output store.")
print("On the receiving machine, run `make activate-legacy` before historical downstream notebooks.")
