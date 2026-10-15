from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from .data import load_inputs


def _load_msaa(repo_root: Path):
    helpers = str(repo_root / "helpers")
    if helpers not in sys.path:
        sys.path.insert(0, helpers)
    from MultiSubject_AA import Subject, multi_subject_aa
    return Subject, multi_subject_aa


def orient_subject_matrix(X_raw, analysis_type):
    Xtv = np.asarray(X_raw, dtype=float)
    return Xtv if analysis_type == "spatial" else Xtv.T


def build_subject_objects(subject_list, analysis_type, Subject):
    return [Subject(X=orient_subject_matrix(x, analysis_type), sX=orient_subject_matrix(x, analysis_type)) for x in subject_list]


def fit_one(subject_list, K, analysis_type, repo_root: Path, maxiter=200, seed=42, verbose=True):
    Subject, multi_subject_aa = _load_msaa(repo_root)
    opts = dict(maxiter=int(maxiter), conv_crit=1e-6, fix_var_iter=5,
                heteroscedastic=False, use_gpu=False, rngSEED=int(seed))
    subj_objs = build_subject_objects(subject_list, analysis_type, Subject)
    results_subj, C, cost_fun, varexpl, secs = multi_subject_aa(subj_objs, noc=int(K), opts=opts)
    if verbose:
        print(f"[{analysis_type}] K={K}: variance explained={varexpl:.6f}, seconds={secs:.1f}")
    return dict(results_subj=results_subj, C=C, cost_fun=cost_fun, varexpl=varexpl, secs=secs, K=int(K))


def fit_filename(analysis_type, fit_scope, K, condition=None):
    if fit_scope == "across":
        return f"{analysis_type}AA_across_acrossCond_K{K}.npz"
    if condition is None:
        raise ValueError("condition is required for within fits")
    return f"{analysis_type}AA_within_{condition}_K{K}.npz"


def save_fit(path: Path, fit_res: dict, labels, analysis_type, fit_scope, centers, widths):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        K=np.array(fit_res["K"]),
        results_subj=np.array(fit_res["results_subj"], dtype=object),
        C=np.array(fit_res["C"], dtype=object),
        cost_fun=np.array(fit_res["cost_fun"], dtype=object),
        varexpl=np.array(fit_res["varexpl"]),
        secs=np.array(fit_res["secs"]),
        condition_labels_str=np.array(labels, dtype=object),
        condition_codes=np.array(pd.Categorical(labels).codes),
        condition_names=np.array(sorted(np.unique(labels)), dtype=object),
        analysis_type=np.array(analysis_type, dtype=object),
        fit_scope=np.array(fit_scope, dtype=object),
        centers=np.array(centers),
        widths=np.array(widths),
    )


def run_fit_job(repo_root: Path, output_root: Path, analysis_type: str, fit_scope: str, K_values,
                *, max_subjects_per_condition=None, maxiter=200, seed=42, skip_existing=True):
    condition_data, centers, widths = load_inputs(
        repo_root, trim_rest_by=100, max_subjects_per_condition=max_subjects_per_condition
    )
    fit_dir = output_root / "msaa_flexible_outputs_npz"
    manifest_rows = []

    if fit_scope == "across":
        all_subjects, labels = [], []
        for condition in ["intact", "word", "rest"]:
            all_subjects.extend(condition_data[condition])
            labels.extend([condition] * len(condition_data[condition]))
        jobs = [(None, all_subjects, labels)]
    else:
        jobs = [(c, condition_data[c], [c] * len(condition_data[c])) for c in ["intact", "word", "rest"]]

    for K in K_values:
        for condition, subjects, labels in jobs:
            out = fit_dir / fit_filename(analysis_type, fit_scope, int(K), condition)
            if skip_existing and out.exists():
                print("Skipping existing", out)
                status = "existing"
            else:
                res = fit_one(subjects, int(K), analysis_type, repo_root, maxiter=maxiter, seed=seed)
                save_fit(out, res, labels, analysis_type, fit_scope, centers, widths)
                status = "fit"
            manifest_rows.append({
                "analysis_type": analysis_type,
                "fit_scope": fit_scope,
                "condition": condition or "acrossCond",
                "K": int(K),
                "n_subjects": len(subjects),
                "path": str(out.relative_to(output_root)),
                "status": status,
            })
    return pd.DataFrame(manifest_rows)
