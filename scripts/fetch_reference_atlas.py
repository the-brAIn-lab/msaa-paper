#!/usr/bin/env python3
"""Fetch the standard Schaefer-1000 / Yeo-7 atlas used for node labeling.

This is a standard reference resource, not an empirical study input. The script
copies the fetched atlas into data/reference/ and writes the simple tabular label
file expected by the historical notebooks.
"""
from __future__ import annotations
from pathlib import Path
import shutil
from nilearn.datasets import fetch_atlas_schaefer_2018

REPO = Path(__file__).resolve().parents[1]
outdir = REPO / "data" / "reference"
outdir.mkdir(parents=True, exist_ok=True)

atlas = fetch_atlas_schaefer_2018(n_rois=1000, yeo_networks=7, resolution_mm=2)
src = Path(atlas.maps)
dst = outdir / "Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
if src.resolve() != dst.resolve():
    shutil.copy2(src, dst)

labels = list(atlas.labels)
# Nilearn versions may include a background label; the NIfTI parcel IDs are 1..1000.
if len(labels) == 1001:
    labels = labels[1:]
if len(labels) != 1000:
    raise RuntimeError(f"Expected 1000 parcel labels, got {len(labels)}")

label_path = outdir / "Schaefer2018_1000Parcels_7Networks_order.txt"
with label_path.open("w") as f:
    for i, label in enumerate(labels, start=1):
        if isinstance(label, bytes):
            label = label.decode("utf-8")
        # Historical notebooks read six tab-separated fields but use only id/name.
        f.write(f"{i}\t{label}\t0\t0\t0\t0\n")

print(dst)
print(label_path)
