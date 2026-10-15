from __future__ import annotations
from pathlib import Path
import importlib.util
import numpy as np
from scipy.io import loadmat

CONDITIONS = ["intact", "word", "rest"]


def _to_float_array(x):
    return np.asarray(x, dtype=float)


def coerce_TxV(x):
    x = _to_float_array(x)
    if x.ndim != 2:
        raise ValueError(f"Expected a 2D subject matrix; got {x.shape}")
    T, V = x.shape
    if V < T:
        x = x.T
    x = x - x.mean(axis=0, keepdims=True)
    return x.astype(np.float32, copy=False)


def _normalize_across_fallback(subjects, eps=1e-8):
    subj = [coerce_TxV(x) for x in subjects]
    stacked = np.vstack(subj)
    mu = stacked.mean(axis=0, keepdims=True)
    sd = stacked.std(axis=0, keepdims=True) + eps
    return [((x - mu) / sd).astype(np.float32, copy=False) for x in subj]


def _load_legacy_normalizer(repo_root: Path):
    path = repo_root / "hypertools" / "normalize.py"
    if not path.exists():
        return None
    try:
        spec = importlib.util.spec_from_file_location("msaa_legacy_normalize", path)
        mod = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(mod)
        return mod
    except Exception:
        return None


def normalize_subject_list(subjects, repo_root: Path):
    subj = [coerce_TxV(x) for x in subjects]
    ht = _load_legacy_normalizer(repo_root)
    if ht is not None:
        try:
            return ht.normalize(subj, normalize="across")
        except Exception:
            pass
    return _normalize_across_fallback(subj)


def load_inputs(repo_root: Path, *, trim_rest_by: int = 100, max_subjects_per_condition=None):
    data_path = repo_root / "data" / "input" / "pieman_data.mat"
    posterior_path = repo_root / "data" / "input" / "pieman_posterior_K700.mat"
    if not data_path.exists() or not posterior_path.exists():
        missing = [str(p) for p in (data_path, posterior_path) if not p.exists()]
        raise FileNotFoundError(
            "Missing required starting inputs:\n  " + "\n  ".join(missing) +
            "\nSee data/README.md."
        )

    pieman = loadmat(data_path)
    condition_data = {}
    for condition in CONDITIONS:
        if condition not in pieman:
            raise KeyError(f"{data_path} does not contain condition {condition!r}")
        arr = pieman[condition]
        subjects = [arr[:, i][0] for i in range(arr.shape[1])]
        if condition == "rest" and trim_rest_by > 0:
            subjects = [np.asarray(x)[:-trim_rest_by] for x in subjects]
        if max_subjects_per_condition is not None:
            subjects = subjects[: int(max_subjects_per_condition)]
        condition_data[condition] = normalize_subject_list(subjects, repo_root)

    posterior = loadmat(posterior_path)

    def unwrap_singleton_object(x):
        # MATLAB structs saved by different scipy/MATLAB versions can have
        # different numbers of singleton object-array wrappers.
        while isinstance(x, np.ndarray) and x.dtype == object and x.size == 1:
            x = x.flat[0]
        return x

    try:
        post = posterior["posterior"]
        if getattr(post, "dtype", None) is not None and post.dtype.names:
            centers_raw = post["centers"][0, 0]
            widths_raw = post["widths"][0, 0]
        else:
            raise TypeError("posterior is not a MATLAB struct")
        centers = _to_float_array(unwrap_singleton_object(centers_raw))
        widths = _to_float_array(unwrap_singleton_object(widths_raw)).ravel()
    except Exception as e:
        raise ValueError(
            "pieman_posterior_K700.mat did not expose posterior.centers/widths in a supported MATLAB-struct layout."
        ) from e
    return condition_data, centers, widths


def describe_inputs(repo_root: Path) -> dict:
    condition_data, centers, widths = load_inputs(repo_root, trim_rest_by=100, max_subjects_per_condition=1)
    data_path = repo_root / "data" / "input" / "pieman_data.mat"
    raw = loadmat(data_path)
    variables = sorted(k for k in raw.keys() if not k.startswith("__"))
    paragraph_aliases = {"paragraph", "paragraph_scrambled", "paragraph-scrambled", "para", "paragraphscrambled"}
    has_paragraph = any(str(k).lower() in paragraph_aliases for k in variables)
    return {
        "mat_variables": variables,
        "four_condition_paragraph_available": has_paragraph,
        "conditions": {k: {"n_loaded_for_check": len(v), "shape": list(np.asarray(v[0]).shape)} for k, v in condition_data.items()},
        "centers_shape": list(np.asarray(centers).shape),
        "widths_shape": list(np.asarray(widths).shape),
    }
