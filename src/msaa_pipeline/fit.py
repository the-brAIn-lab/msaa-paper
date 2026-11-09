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


INIT_PATH = Path("configs/furthestsum_init.csv.gz")


def furthest_sum_picks(subject_list, K, analysis_type, repo_root: Path, seed=42):
    """The FurthestSum picks multi_subject_aa would compute for this fit.

    Same data stacking, rng draw and call as the initialisation inside
    helpers/MultiSubject_AA.py, so passing the result back as init_idx gives a
    bit-identical fit.
    """
    helpers = str(repo_root / "helpers")
    if helpers not in sys.path:
        sys.path.insert(0, helpers)
    from furthest_sum import furthest_sum
    X = np.vstack([orient_subject_matrix(x, analysis_type) for x in subject_list])
    first = np.random.default_rng(int(seed)).integers(low=0, high=X.shape[1])
    return [int(i) for i in furthest_sum(X, noc=int(K), i=int(first), exclude=None, treat_as_kernel=False)]


def load_frozen_init(repo_root: Path, analysis_type, fit_scope, K, condition=None):
    """Frozen FurthestSum picks for one fit from configs/furthestsum_init.csv.gz.

    Freezing the picks makes archetype numbering identical on every machine:
    archetype j is always the one initialised at row/column init_index[j], and
    with the post-fit sort off, nothing reorders the columns afterwards.
    """
    path = repo_root / INIT_PATH
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; create it with rerun/freeze_init.py")
    df = pd.read_csv(path)
    cond = "acrossCond" if fit_scope == "across" else condition
    sub = df[(df.analysis_type == analysis_type) & (df.fit_scope == fit_scope)
             & (df.condition == cond) & (df.K == int(K))].sort_values("archetype")
    if len(sub) != int(K):
        raise ValueError(f"No frozen FurthestSum picks for {analysis_type}/{fit_scope}/{cond}/K={K} in {path}")
    return sub["init_index"].astype(int).tolist()


def fit_one(subject_list, K, analysis_type, repo_root: Path, maxiter=200, seed=42, verbose=True,
            init_idx=None, sort_crit="corr"):
    Subject, multi_subject_aa = _load_msaa(repo_root)
    opts = dict(maxiter=int(maxiter), conv_crit=1e-6, fix_var_iter=5,
                heteroscedastic=False, use_gpu=False, rngSEED=int(seed), sort_crit=sort_crit)
    if init_idx is not None:
        opts["init_idx"] = list(init_idx)
    subj_objs = build_subject_objects(subject_list, analysis_type, Subject)
    results_subj, C, cost_fun, varexpl, secs = multi_subject_aa(subj_objs, noc=int(K), opts=opts)
    if verbose:
        print(f"[{analysis_type}] K={K}: variance explained={varexpl:.6f}, seconds={secs:.1f}")
    return dict(results_subj=results_subj, C=C, cost_fun=cost_fun, varexpl=varexpl, secs=secs, K=int(K),
                init_index=None if init_idx is None else np.asarray(init_idx, dtype=int), sort_crit=sort_crit)


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
        **({"init_index": fit_res["init_index"]} if fit_res.get("init_index") is not None else {}),
        **({"sort_crit": np.array(fit_res["sort_crit"], dtype=object)} if "sort_crit" in fit_res else {}),
    )


def run_fit_job(repo_root: Path, output_root: Path, analysis_type: str, fit_scope: str, K_values,
                *, max_subjects_per_condition=None, maxiter=200, seed=42, skip_existing=True,
                frozen_init=False):
    """Fit every K (and condition, for within fits) of one job.

    frozen_init=True uses the frozen FurthestSum picks in configs/furthestsum_init.csv.gz and no
    post-fit sort, exactly as scripts/run_parallel.py does, so archetype numbering matches the
    paper's (canonical) fits. frozen_init=False recomputes FurthestSum and sorts archetypes after
    fitting (the legacy numbering of the original analysis outputs).
    """
    if frozen_init and max_subjects_per_condition is not None:
        raise ValueError("frozen FurthestSum picks were drawn on the full data; "
                         "use frozen_init=False with max_subjects_per_condition")
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
                if frozen_init:
                    init_idx = load_frozen_init(repo_root, analysis_type, fit_scope, int(K), condition)
                    res = fit_one(subjects, int(K), analysis_type, repo_root, maxiter=maxiter, seed=seed,
                                  init_idx=init_idx, sort_crit="none")
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
