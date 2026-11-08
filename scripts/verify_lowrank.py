#!/usr/bin/env python3
"""Check the low-rank decoder against the validated rerun's dense decoder.

1. The vectorised _decode_corr gives bit-identical results to the original
   loop (dense path, new module vs the validated rerun's module).
2. Low-rank vs dense correlation matrices agree up to rounding, and the
   per-fold accuracy / rank / error match, for single archetypes and top-m
   subsets, spatial and temporal, across and within fits.
3. Single-archetype summaries reproduce the validated rerun's
   per_archetype_mean_accuracy.csv.

Usage (from the workspace root):
  python scripts/verify_lowrank.py [--quick]
  MSAA_REF_RUN=/path/to/reference_run (a run root with the dense-decode code and outputs)
"""
from __future__ import annotations
import argparse, importlib.util, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = Path(os.environ.get("MSAA_REF_RUN", ROOT))
sys.path.insert(0, str(ROOT / "src"))
from msaa_pipeline import decode as new  # noqa: E402

spec = importlib.util.spec_from_file_location("ref_decode", REF / "src/msaa_pipeline/decode.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)

CASES = [  # analysis_type, fit_scope, K, condition
    ("spatial", "across", 46, "intact"),
    ("spatial", "across", 300, "rest"),
    ("spatial", "within", 14, "intact"),
    ("spatial", "within", 252, "rest"),
    ("temporal", "across", 45, "word"),
    ("temporal", "across", 300, "intact"),
    ("temporal", "within", 35, "word"),
    ("temporal", "within", 300, "rest"),
]


def dense_corrs(subjects, analysis_type, idx, folds_list, nfolds):
    stack = np.stack([ref.reconstruct(s, analysis_type, idx) for s in subjects])
    out = []
    for folds in folds_list:
        for fold in range(nfolds):
            c = ref._row_corr(stack[folds == fold].mean(axis=0), stack[folds != fold].mean(axis=0))
            out.append(np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fewer archetypes per case")
    args = ap.parse_args()
    fits = REF / "msaa_flexible_outputs_npz"
    report = []
    for at, sc, K, cond in CASES:
        subjects, _ = new.get_subjects(REF, at, sc, K, cond)
        n_arch = K if not args.quick else min(K, 8)
        arch_ids = list(range(n_arch)) if n_arch == K else list(np.linspace(0, K - 1, n_arch).astype(int))
        refcsv = pd.read_csv(REF / f"msaa_condrank_decoding_outputs_{at}_{sc}" / "per_archetype_mean_accuracy.csv")
        refcsv = refcsv[(refcsv.K == K) & (refcsv.condition == cond)].set_index("archetype")

        # 1. dense path: new module (vectorised scoring) vs reference module, bit for bit
        t = time.time()
        d_ref = ref.decode_subjects(subjects, at, None, 2, 20, 42)
        t_full = time.time() - t
        d_new = new.decode_subjects(subjects, at, None, 2, 20, 42)
        full_identical = d_ref.equals(d_new)

        # 2./3. single archetypes and top-m subsets: low-rank vs dense
        max_corr_diff = 0.0
        fold_mismatch = 0
        fold_mismatch_tol = 0
        max_dacc = max_drank = 0.0
        csv_diff = 0.0
        t_dense = t_low = 0.0
        subsets = [[a] for a in arch_ids]
        rng = np.random.default_rng(0)
        for m in (2, 3, 5, 10):
            if m <= K:
                subsets.append(sorted(rng.choice(K, m, replace=False).tolist()))
        for idx in subsets:
            t = time.time(); dd = ref.decode_subjects(subjects, at, idx, 2, 20, 42); t_dense += time.time() - t
            new.TIE_TOL = 0.0
            t = time.time(); dl = new.decode_subjects(subjects, at, idx, 2, 20, 42); t_low += time.time() - t
            cols = ["accuracy", "rank", "error"]
            fold_mismatch += int((np.abs(dd[cols].to_numpy() - dl[cols].to_numpy()) > 1e-12).any(axis=1).sum())
            max_dacc = max(max_dacc, abs(dd.accuracy.mean() - dl.accuracy.mean()))
            max_drank = max(max_drank, abs(dd["rank"].mean() - dl["rank"].mean()))
            new.TIE_TOL = 1e-8
            dt = new.decode_subjects(subjects, at, idx, 2, 20, 42, method="dense")
            dlt = new.decode_subjects(subjects, at, idx, 2, 20, 42)
            new.TIE_TOL = 0.0
            fold_mismatch_tol += int((np.abs(dt[cols].to_numpy() - dlt[cols].to_numpy()) > 1e-12).any(axis=1).sum())
            if len(idx) == 1:
                csv_diff = max(csv_diff, abs(dl.accuracy.mean() - refcsv.loc[idx[0], "mean_accuracy"]),
                               abs(dl["rank"].mean() - refcsv.loc[idx[0], "rank_mean"]))
        # correlation matrices themselves, for a few subsets
        frng = np.random.default_rng(42)
        folds_list = [new._assignments(len(subjects), 2, frng) for _ in range(3)]
        for idx in subsets[:3] + subsets[-2:]:
            P = [new._lowrank_factors(s, at, np.asarray(idx)) for s in subjects]
            Pm = np.concatenate([p[0] for p in P], axis=1); Wm = np.concatenate([p[1] for p in P], axis=1)
            Wc = Wm - Wm.mean(axis=0); G = Wc.T @ Wc; m = len(idx)
            dcs = dense_corrs(subjects, at, idx, folds_list, 2)
            k = 0
            for folds in folds_list:
                cols = np.repeat(folds, m)
                for fold in range(2):
                    ia, ib = np.flatnonzero(cols == fold), np.flatnonzero(cols != fold)
                    num = Pm[:, ia] @ G[np.ix_(ia, ib)] @ Pm[:, ib].T
                    na = np.sqrt(np.einsum("tr,rt->t", Pm[:, ia], G[np.ix_(ia, ia)] @ Pm[:, ia].T))
                    nb = np.sqrt(np.einsum("tr,rt->t", Pm[:, ib], G[np.ix_(ib, ib)] @ Pm[:, ib].T))
                    c = np.nan_to_num(num / np.outer(na, nb))
                    max_corr_diff = max(max_corr_diff, float(np.abs(c - dcs[k]).max())); k += 1
        r = dict(case=f"{at}_{sc}_K{K}_{cond}", n_subj=len(subjects), n_single=len(arch_ids),
                 n_topm=len(subsets) - len(arch_ids), full_dense_bit_identical=full_identical,
                 max_abs_corr_diff=max_corr_diff, folds_with_metric_diff=fold_mismatch,
                 max_dmean_acc=max_dacc, max_dmean_rank=max_drank, folds_diff_tol1e8=fold_mismatch_tol,
                 max_diff_vs_validated_csv=csv_diff, dense_s=round(t_dense, 2), lowrank_s=round(t_low, 2),
                 speedup=round(t_dense / max(t_low, 1e-9), 1), full_decode_s=round(t_full, 2))
        print(r, flush=True)
        report.append(r)
    out = pd.DataFrame(report)
    out.to_csv(ROOT / "rerun" / "verify_lowrank.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
