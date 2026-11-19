"""Convert pieman2_htfa.pkl into the two HTFA input files the pipeline reads.

Usage:
    python scripts/pkl_to_mat.py <path-to-pieman2_htfa.pkl> [--out data/input]

Writes, in the layout that src/msaa_pipeline/data.py and the analysis scripts expect:
  - pieman_data.mat: 1 x N cell arrays intact, paragraph, rest and word, one T x 700 array
    of HTFA weights per scan.
  - pieman_posterior_K700.mat: struct `posterior` with centers.mean (700 x 3) and
    widths.mean (700 x 1).

The pkl's paragraph-scrambled condition holds 17 scans because a duplicated entry was
removed. The original MAT file kept that duplicate at column 3 (identical to column 9),
and the four-condition analysis (analysis/06_four_condition) drops column 3 by position.
To keep every script working unchanged, the converter restores the original 18-column
layout by placing a copy of the duplicated scan back at column 3.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import savemat

CONDITIONS = ["intact", "paragraph", "rest", "word"]
PARAGRAPH_DUPLICATE_COLUMN = 3  # in the original MAT file
PARAGRAPH_DUPLICATE_SOURCE = 8  # position of the same scan in the pkl's 17-scan list


def to_cell(scans):
    cell = np.empty((1, len(scans)), dtype=object)
    for i, s in enumerate(scans):
        cell[0, i] = s
    return cell


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pkl", type=Path, help="path to pieman2_htfa.pkl")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "data" / "input",
                    help="output directory (default: data/input)")
    a = ap.parse_args()

    htfa = pd.read_pickle(a.pkl)
    weights = htfa["weights"]

    data = {}
    for c in CONDITIONS:
        scans = [np.asarray(x, dtype=np.float64) for x in weights[c]]
        if c == "paragraph":
            if len(scans) != 17:
                raise ValueError(f"expected 17 paragraph scans in the pkl, found {len(scans)}")
            scans.insert(PARAGRAPH_DUPLICATE_COLUMN, scans[PARAGRAPH_DUPLICATE_SOURCE].copy())
        data[c] = to_cell(scans)

    centers = np.asarray(htfa["centers"], dtype=np.float64)
    widths = np.asarray(htfa["widths"], dtype=np.float64).reshape(-1, 1)
    posterior = {"posterior": {"centers": {"mean": centers}, "widths": {"mean": widths}}}

    a.out.mkdir(parents=True, exist_ok=True)
    savemat(a.out / "pieman_data.mat", data)
    savemat(a.out / "pieman_posterior_K700.mat", posterior)
    for c in CONDITIONS:
        print(f"{c}: {data[c].shape[1]} scans, {data[c][0, 0].shape}")
    print(f"centers {centers.shape}, widths {widths.shape} -> {a.out}")


if __name__ == "__main__":
    main()
