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

    def unwrap_singletons(x):
        """Remove MATLAB/scipy singleton wrappers regardless of ndarray dtype.

        The original analysis file has several nested 1x1 arrays around the
        numeric centers/widths arrays. Restricting this to dtype=object is too
        narrow because scipy can expose some wrappers as structured/numeric
        ndarrays depending on scipy/MATLAB versions.
        """
        while isinstance(x, np.ndarray) and x.size == 1:
            x = x.reshape(-1)[0]
        return x

    def extract_field_recursive(mat_struct, field):
        """Extract one field from a scipy-loaded MATLAB struct and unwrap 1x1 layers."""
        x = mat_struct
        # Unwrap outer singleton arrays until we reach either a structured
        # scalar/array or the actual numeric value.
        x = unwrap_singletons(x)

        # scipy can represent a MATLAB struct either as a structured ndarray
        # or as a numpy.void structured scalar after singleton unwrapping.
        if isinstance(x, np.void) and x.dtype.names and field in x.dtype.names:
            return unwrap_singletons(x[field])
        if isinstance(x, np.ndarray) and x.dtype.names and field in x.dtype.names:
            return unwrap_singletons(x[field])

        raise KeyError(f"MATLAB struct does not contain field {field!r}")

    extraction_errors = []
    centers = widths = None

    # First reproduce the exact indexing used by the original fitting notebook.
    # Keeping this path is important because these two MAT files define the
    # reproducibility boundary for the project.
    try:
        centers = _to_float_array(
            posterior["posterior"]["centers"][0][0][0][0][0]
        )
        widths = _to_float_array(
            posterior["posterior"]["widths"][0][0][0][0][0]
        ).ravel()
    except Exception as e:
        extraction_errors.append(f"legacy extraction: {type(e).__name__}: {e}")

    # Fallback for equivalent MAT files loaded with a scipy/MATLAB version that
    # exposes a slightly different number/type of singleton wrappers.
    if centers is None or widths is None:
        try:
            post = posterior["posterior"]
            centers = _to_float_array(extract_field_recursive(post, "centers"))
            widths = _to_float_array(extract_field_recursive(post, "widths")).ravel()
        except Exception as e:
            extraction_errors.append(f"recursive extraction: {type(e).__name__}: {e}")
            centers = widths = None

    if centers is None or widths is None:
        top_vars = sorted(k for k in posterior.keys() if not k.startswith("__"))
        post = posterior.get("posterior")
        post_desc = (
            f"type={type(post).__name__}, shape={getattr(post, 'shape', None)}, "
            f"dtype={getattr(post, 'dtype', None)}, "
            f"fields={getattr(getattr(post, 'dtype', None), 'names', None)}"
        )
        raise ValueError(
            "Could not extract centers/widths from pieman_posterior_K700.mat.\n"
            f"Top-level variables: {top_vars}\n"
            f"posterior structure: {post_desc}\n"
            + "\n".join(extraction_errors)
        )

    centers = np.asarray(centers, dtype=float).squeeze()
    widths = np.asarray(widths, dtype=float).squeeze().ravel()

    # Canonical orientation is K x xyz (normally 700 x 3).
    if centers.ndim == 1 and widths.size == 1:
        centers = centers.reshape(1, -1)
    if centers.ndim != 2:
        raise ValueError(f"Expected centers to be 2D after extraction; got {centers.shape}")
    if centers.shape[0] != widths.size and centers.shape[1] == widths.size:
        centers = centers.T
    if centers.shape[0] != widths.size:
        raise ValueError(
            f"Posterior centers/widths disagree: centers={centers.shape}, widths={widths.shape}"
        )

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
