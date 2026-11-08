#!/usr/bin/env python3
"""Build configs/node_yeo7.csv: the Yeo-7 network label of each of the 700 HTFA nodes.

Each HTFA node is a Gaussian (RBF) source with a centre and a width
(data/input/pieman_posterior_K700.mat). Its label is the Yeo-7 network with the most RBF
weight over the voxels of the Schaefer-2018 1000-parcel / 7-network atlas (fetched by
scripts/fetch_reference_atlas.py). This is the same rule the interpretation notebooks
use (build_node_codes in analysis/04_interpretation/revision_06_*.ipynb), so the table
reproduces the labels behind the paper figures.

Usage (from the repo root):
  python scripts/fetch_reference_atlas.py          # once; needs network access
  python scripts/build_node_yeo7.py [--check LEGACY_allnodes.csv]

Output columns: node (0-based), x, y, z (MNI mm, node centre), width, network, code (1-7,
in the order Visual, Somatomotor, Dorsal attention, Ventral attention, Limbic,
Frontoparietal, Default mode).
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

REPO = Path(__file__).resolve().parents[1]
POSTERIOR_MAT = REPO / "data/input/pieman_posterior_K700.mat"
SCHAEFER_TXT = REPO / "data/reference/Schaefer2018_1000Parcels_7Networks_order.txt"
SCHAEFER_NII = REPO / "data/reference/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
OUT = REPO / "configs/node_yeo7.csv"

NETWORK_ORDER = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic",
                 "Frontoparietal", "Default mode"]
LOOKUP = {"Vis": "Visual", "SomMot": "Somatomotor", "DorsAttn": "Dorsal attention",
          "SalVentAttn": "Ventral attention", "Limbic": "Limbic", "Cont": "Frontoparietal",
          "Default": "Default mode"}
CODES = {name: i + 1 for i, name in enumerate(NETWORK_ORDER)}


def node_centres_widths(mat=POSTERIOR_MAT):
    post = loadmat(mat)["posterior"]
    centers = np.asarray(post["centers"][0][0][0][0][0], dtype=float)
    widths = np.asarray(list(post["widths"][0][0][0][0][0][:, 0].T), dtype=float).ravel()
    return centers, widths


def atlas_voxels(nii=SCHAEFER_NII, txt=SCHAEFER_TXT):
    """Voxel coordinates (mm) and Yeo-7 code of every in-mask atlas voxel, as in the notebooks."""
    import nibabel as nib
    try:
        from nilearn.maskers import NiftiMasker
    except ImportError:
        from nilearn.input_data import NiftiMasker

    key = pd.read_csv(txt, sep="\t", header=None, names=["id", "name", "x", "y", "z", "t"])
    key["code"] = key["name"].map(lambda n: CODES[LOOKUP[n.split("_")[2]]])
    code_of = dict(zip(key["id"], key["code"])); code_of[0] = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        img = nib.load(str(nii))
        masker = NiftiMasker(mask_strategy="background").fit(str(nii))
        labels = np.asarray(masker.transform(str(nii)), dtype=np.float32).ravel()
    aff = img.get_sform()
    mask_flat = np.asarray(masker.mask_img_.dataobj).reshape(-1, order="C")
    vmask = np.nonzero(mask_flat)[0]
    # index grid with the first axis varying fastest (MATLAB-style fullfact), columns reversed,
    # 0-based -- kept exactly as in the notebooks so the labels match the published figures
    dims = img.shape[:3]
    grid = np.stack(np.unravel_index(np.arange(np.prod(dims)), dims, order="F"), axis=1) + 1
    vox = grid[vmask, ::-1] - 1
    r = vox @ aff[:3, :3] + aff[:3, 3]
    codes = np.array([code_of[int(v)] for v in labels])
    return r, codes


def label_nodes(centers, widths, r, codes):
    out = []
    for c, w in zip(centers, widths):
        wt = np.exp(-np.sum((r - c) ** 2, axis=1) / w)
        out.append(int(np.argmax([wt[codes == k].sum() for k in range(1, 8)]) + 1))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--check", default=None, help="legacy table with node,x,y,z,network columns to compare against")
    a = ap.parse_args()
    for p in (POSTERIOR_MAT, SCHAEFER_NII, SCHAEFER_TXT):
        if not p.exists():
            raise SystemExit(f"missing {p} (run scripts/fetch_reference_atlas.py for the atlas files)")
    centers, widths = node_centres_widths()
    r, codes = atlas_voxels()
    lab = label_nodes(centers, widths, r, codes)
    T = pd.DataFrame({"node": np.arange(len(centers)), "x": centers[:, 0], "y": centers[:, 1], "z": centers[:, 2],
                      "width": widths, "network": [NETWORK_ORDER[k - 1] for k in lab], "code": lab})
    T.to_csv(a.out, index=False)
    print(f"wrote {a.out}: {len(T)} nodes"); print(T.network.value_counts().reindex(NETWORK_ORDER).to_string())
    if a.check:
        L = pd.read_csv(a.check).drop_duplicates("node").set_index("node").sort_index()
        same_net = (L["network"].to_numpy() == T["network"].to_numpy())
        dxyz = np.abs(L[["x", "y", "z"]].to_numpy() - T[["x", "y", "z"]].to_numpy()).max()
        print(f"check vs {a.check}: same network {same_net.sum()}/{len(T)}, max |xyz diff| {dxyz:.2e}")
        if not same_net.all() or dxyz > 1e-6:
            raise SystemExit("node table does not match the legacy labels")


if __name__ == "__main__":
    main()
