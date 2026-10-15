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
        if len(m):
            diffs[c] = float(np.nanmax(np.abs(m[f"{c}_new"] - m[f"{c}_old"])))
    return {
        "status": "ok",
        "matched_rows": int(len(m)),
        "max_abs_diff": diffs,
        "within_tolerance": bool(diffs and all(v <= atol for v in diffs.values())),
        "atol": float(atol),
    }


def compare_outputs(new_root: Path, old_root: Path):
    results = {"decoding": {}, "heldout_topm": {}, "fits": {}}

    # Core comparison: full reconstruction and single-archetype diagnostics only.
    # Historical same-data top-m outputs are intentionally ignored.
    for d in new_root.glob("msaa_condrank_decoding_outputs_*"):
        oldd = old_root / d.name
        results["decoding"][d.name] = {
            "full": _compare_csv(
                d/"full_summary.csv", oldd/"full_summary.csv",
                ["analysis_type","fit_scope","condition","K"],
                ["mean","err","rank_mean","error_mean"], atol=5e-4,
            ),
            "per_archetype": _compare_csv(
                d/"per_archetype_mean_accuracy.csv", oldd/"per_archetype_mean_accuracy.csv",
                ["analysis_type","fit_scope","condition","K","archetype"],
                ["mean_accuracy","sem_accuracy","rank_mean","error_mean"], atol=5e-3,
            ),
        }

    # Canonical sparse comparison: held-out A->B/B->A cross-fitting only.
    for d in new_root.glob("msaa_heldout_topm_outputs_*_across"):
        oldd = old_root / d.name
        results["heldout_topm"][d.name] = {
            "crossfit_runs": _compare_csv(
                d/"heldout_topm_crossfit_runs.csv", oldd/"heldout_topm_crossfit_runs.csv",
                ["analysis_type","fit_scope","condition","K","outer_rep","direction","top_m"],
                ["mean","full_mean_same_eval","retention_above_chance"], atol=5e-3,
            )
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
