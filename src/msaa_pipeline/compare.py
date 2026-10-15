from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from .decode import load_fit, reconstruct


def _compare_csv(new_path: Path, old_path: Path, key_cols, value_cols, atol=1e-10):
    if not new_path.exists() or not old_path.exists():
        return {"status": "missing", "new": str(new_path), "old": str(old_path)}
    new = pd.read_csv(new_path); old = pd.read_csv(old_path)
    keys = [c for c in key_cols if c in new.columns and c in old.columns]
    vals = [c for c in value_cols if c in new.columns and c in old.columns]
    m = new.merge(old, on=keys, suffixes=("_new", "_old")) if keys else pd.DataFrame()
    diffs = {}
    for c in vals:
        if len(m): diffs[c] = float(np.nanmax(np.abs(m[f"{c}_new"] - m[f"{c}_old"])))
    return {"status": "ok", "matched_rows": int(len(m)), "max_abs_diff": diffs,
            "within_tolerance": bool(diffs and all(v <= atol for v in diffs.values()))}


def compare_outputs(new_root: Path, old_root: Path):
    results = {"decoding": {}, "fits": {}}
    for d in new_root.glob("msaa_condrank_decoding_outputs_*"):
        oldd = old_root / d.name
        results["decoding"][d.name] = {
            "full": _compare_csv(d/"full_summary.csv", oldd/"full_summary.csv",
                                 ["analysis_type","fit_scope","condition","K"], ["mean","err","rank_mean","error_mean"], atol=1e-8),
            "per_archetype": _compare_csv(d/"per_archetype_mean_accuracy.csv", oldd/"per_archetype_mean_accuracy.csv",
                                 ["analysis_type","fit_scope","condition","K","archetype"], ["mean_accuracy","sem_accuracy","rank_mean","error_mean"], atol=1e-8),
            "topm": _compare_csv(d/"topm_summary.csv", oldd/"topm_summary.csv",
                                 ["analysis_type","fit_scope","condition","K","top_m"], ["mean","err","rank_mean","error_mean"], atol=1e-8),
        }

    newfits = new_root / "msaa_flexible_outputs_npz"
    oldfits = old_root / "msaa_flexible_outputs_npz"
    if newfits.exists() and oldfits.exists():
        for p in sorted(newfits.glob("*.npz")):
            q = oldfits / p.name
            if not q.exists():
                results["fits"][p.name] = {"status":"missing_old"}; continue
            nf, of = load_fit(p), load_fit(q)
            nsub, osub = nf["results_subj"][0], of["results_subj"][0]
            # infer analysis from filename; reconstruction is invariant to archetype permutation
            typ = "spatial" if p.name.startswith("spatial") else "temporal"
            nr, orr = reconstruct(nsub, typ), reconstruct(osub, typ)
            nd = np.load(p, allow_pickle=True); od = np.load(q, allow_pickle=True)
            nv = float(np.asarray(nd["varexpl"]).squeeze()) if "varexpl" in nd else None
            ov = float(np.asarray(od["varexpl"]).squeeze()) if "varexpl" in od else None
            results["fits"][p.name] = {
                "status":"ok",
                "shape_equal": nr.shape == orr.shape,
                "variance_explained_new": nv,
                "variance_explained_old": ov,
                "variance_explained_abs_diff": abs(nv-ov) if nv is not None and ov is not None else None,
                "reconstruction_max_abs_diff": float(np.max(np.abs(nr-orr))) if nr.shape == orr.shape else None,
            }
    return results
