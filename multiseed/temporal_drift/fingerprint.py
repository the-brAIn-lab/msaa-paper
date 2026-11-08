#!/usr/bin/env python3
"""Compact fingerprint of a fit folder, small enough to copy between machines (the full fits are
~4.8 GB). Per fit it keeps everything compare_hellgate.py needs:
  C (shared generator, float64), the subject-mean S (float32), subject 0's S and sXC (float32),
  init_index, varexpl, final cost, number of subjects.

  python fingerprint.py --root RUN_ROOT --out fp.npz [--only name1,name2]
"""
import argparse
from pathlib import Path
import numpy as np


def fit_fingerprint(path):
    z = np.load(path, allow_pickle=True)
    subj = z["results_subj"].tolist()
    S = [np.asarray(s["S"], float) for s in subj]
    cost = np.asarray(z["cost_fun"].tolist(), float).ravel()
    return {
        "C": np.asarray(z["C"].tolist(), float),
        "Smean": np.mean(S, axis=0).astype(np.float32),
        "S0": S[0].astype(np.float32),
        "sXC0": np.asarray(subj[0]["sXC"], float).astype(np.float32),
        "init_index": np.asarray(z["init_index"] if "init_index" in z.files else [-1]),
        "varexpl": np.asarray(float(z["varexpl"])),
        "cost_last": np.asarray(cost[-1] if cost.size else np.nan),
        "n_subj": np.asarray(len(subj)),
    }


def build(root, only=None):
    """{fit_stem: fingerprint dict} for every fit under root/msaa_flexible_outputs_npz."""
    fit_dir = Path(root) / "msaa_flexible_outputs_npz"
    out = {}
    for p in sorted(fit_dir.glob("*.npz")):
        if ".partial" in p.name or (only and p.stem not in only):
            continue
        out[p.stem] = fit_fingerprint(p)
    return out


def save(fp, path):
    np.savez_compressed(path, **{f"{fit}|{k}": v for fit, d in fp.items() for k, v in d.items()})


def load(path):
    z = np.load(path)
    fp = {}
    for key in z.files:
        fit, k = key.split("|")
        fp.setdefault(fit, {})[k] = z[key]
    return fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="run root containing msaa_flexible_outputs_npz/")
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default=None, help="comma list of fit names (without .npz)")
    a = ap.parse_args()
    only = set(a.only.split(",")) if a.only else None
    fp = build(a.root, only)
    if only and set(fp) != only:
        raise SystemExit(f"missing fits: {sorted(only - set(fp))}")
    save(fp, a.out)
    print(f"{len(fp)} fits -> {a.out} ({Path(a.out).stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
