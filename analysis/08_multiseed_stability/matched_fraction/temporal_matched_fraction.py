#!/usr/bin/env python3
"""Matched-retained-fraction analysis for temporal MS-AA.

At a fixed FRACTION of the temporal basis (m = round(fraction * K)), how much of the full model's
above-chance held-out accuracy do the top-m archetypes retain, and does that peak at some K?
Three changes from the version first written for this question keep it identical to the held-out
top-m pipeline (validate.py checks this): the 1e-8 decoder tie tolerance, the archetype ranking order,
and keeping one fit per seed even when two seeds converged to the same fit.

The workflow has four subcommands:

1. build-manifests
   Discover temporal fit NPZ files and create fit/task manifests.
2. audit
   Check raw held-out ranking coverage for each task.
3. run-task
   Run or resume one fit-seed x K x condition task. Each task processes all
   outer splits and both ranking/evaluation directions, checkpointing after
   every completed direction.
4. aggregate
   Merge task outputs and summarize direction -> outer split -> fit seed.

The script never refits MS-AA. It expects existing temporal fit NPZ files and
raw participant-held-out per-archetype ranking rows.
"""

from __future__ import annotations

import argparse
import ast
import glob
import hashlib
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


VERSION = "1.0.0"
DEFAULT_K_VALUES = [5, 8, 10, 12, 14, 16, 20, 25, 30, 40, 50]
DEFAULT_FRACTIONS = [0.25, 0.50, 0.70, 1.00]
DEFAULT_CONDITIONS = ["intact", "word"]
ALL_CONDITIONS = ["intact", "word", "rest"]
CONDITION_INDEX = {condition: i for i, condition in enumerate(ALL_CONDITIONS)}
FIT_FILENAME_RE = re.compile(r"temporalAA_across_acrossCond_K(\d+)\.npz$")
SEED_RE = re.compile(r"seed[_-]?(\d+)", flags=re.IGNORECASE)
EPS = 1e-8


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------


def parse_csv_ints(text: str) -> List[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def parse_csv_floats(text: str) -> List[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_csv_strings(text: str) -> List[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def expand_paths(patterns: Optional[Sequence[str]]) -> List[Path]:
    if not patterns:
        return []
    paths: List[Path] = []
    seen = set()
    for pattern in patterns:
        matches = glob.glob(pattern, recursive=True)
        if not matches and Path(pattern).exists():
            matches = [pattern]
        for match in matches:
            path = Path(match).expanduser().resolve()
            if path in seen or not path.exists() or not path.is_file():
                continue
            seen.add(path)
            paths.append(path)
    return paths


def atomic_write_csv(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=path.name + ".",
        suffix=".tmp",
        dir=str(path.parent),
    )
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        df.to_csv(tmp_path, index=False)
        tmp_path.replace(path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def load_csv(path: Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path)


def sem(values: Iterable[float]) -> float:
    series = pd.Series(values, dtype=float).dropna()
    if len(series) <= 1:
        return float("nan")
    return float(series.std(ddof=1) / math.sqrt(len(series)))


def normalize_labels(values: np.ndarray) -> np.ndarray:
    labels = []
    for value in np.asarray(values).reshape(-1):
        if isinstance(value, bytes):
            labels.append(value.decode("utf-8"))
        else:
            labels.append(str(value))
    return np.asarray(labels, dtype=object)


def normalize_seed_value(value: object) -> Optional[int]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value)
    match = re.search(r"(?:seed[_-]?)?(\d+)$", text, flags=re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def safe_literal_list(value: object) -> Optional[List[int]]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (list, tuple, np.ndarray)):
        return [int(x) for x in value]
    try:
        parsed = ast.literal_eval(str(value))
    except Exception:
        return None
    if isinstance(parsed, (list, tuple, np.ndarray)):
        return [int(x) for x in parsed]
    return None


# -----------------------------------------------------------------------------
# Fit discovery
# -----------------------------------------------------------------------------


def load_fit(path: Path) -> Dict[str, object]:
    path = Path(path)
    data = np.load(path, allow_pickle=True)
    if "results_subj" not in data or "K" not in data:
        raise ValueError(f"Fit file lacks required arrays: {path}")
    results = data["results_subj"].tolist()
    if isinstance(results, np.ndarray):
        results = results.tolist()
    labels = (
        normalize_labels(data["condition_labels_str"])
        if "condition_labels_str" in data
        else None
    )
    return {
        "path": path,
        "K": int(np.asarray(data["K"]).item()),
        "results_subj": results,
        "condition_labels_str": labels,
        "C": np.asarray(data["C"], dtype=float) if "C" in data else None,
    }


def fit_fingerprint(path: Path) -> str:
    fit = load_fit(path)
    digest = hashlib.sha1()
    C = fit.get("C")
    if C is not None:
        digest.update(np.round(np.asarray(C, dtype=float), 10).tobytes())
    else:
        first = fit["results_subj"][0]
        digest.update(
            np.round(np.asarray(first["sXC"], dtype=float), 10).tobytes()
        )
    return digest.hexdigest()[:16]


def infer_fit_seed(path: Path) -> Optional[int]:
    match = SEED_RE.search(str(path))
    if match:
        return int(match.group(1))
    # The historical canonical fit directory is the seed-0 paper fit.
    if path.parent.name == "msaa_flexible_outputs_npz":
        return 0
    return None


def infer_fit_label(path: Path, fingerprint: str) -> str:
    seed = infer_fit_seed(path)
    if seed is not None:
        return f"seed{seed}"
    parent = re.sub(r"[^A-Za-z0-9_.-]+", "_", path.parent.name).strip("_")
    return parent or f"fit_{fingerprint}"


def discover_fit_manifest(
    fit_roots: Sequence[Path],
    k_values: Sequence[int],
) -> pd.DataFrame:
    rows = []
    seen_paths = set()
    seen_fingerprints = set()
    allowed_k = set(int(k) for k in k_values)

    for root in fit_roots:
        root = Path(root).expanduser().resolve()
        if not root.exists():
            continue
        candidates: Iterable[Path]
        if root.is_file():
            candidates = [root]
        else:
            candidates = root.glob("**/temporalAA_across_acrossCond_K*.npz")
        for path in candidates:
            path = Path(path).resolve()
            if path in seen_paths:
                continue
            seen_paths.add(path)
            match = FIT_FILENAME_RE.search(path.name)
            if not match:
                continue
            K = int(match.group(1))
            if K not in allowed_k:
                continue
            try:
                fit = load_fit(path)
                fingerprint = fit_fingerprint(path)
            except Exception as exc:
                print(f"WARNING: skipping unreadable fit {path}: {exc}", file=sys.stderr)
                continue
            fingerprint = f"{infer_fit_label(path, fingerprint)}-{fingerprint}"  # keep seeds whose fits are identical (one task per seed)
            if fingerprint in seen_fingerprints:
                continue
            seen_fingerprints.add(fingerprint)
            seed = infer_fit_seed(path)
            label = infer_fit_label(path, fingerprint)
            labels = fit.get("condition_labels_str")
            condition_counts = {}
            if labels is not None:
                label_series = pd.Series(labels.astype(str))
                condition_counts = label_series.value_counts().to_dict()
            rows.append(
                {
                    "fit_label": label,
                    "fit_seed": seed,
                    "K": int(fit["K"]),
                    "fit_path": str(path),
                    "fit_fingerprint": fingerprint,
                    "n_observations": len(fit["results_subj"]),
                    "condition_counts_json": json.dumps(condition_counts, sort_keys=True),
                }
            )

    manifest = pd.DataFrame(rows)
    if len(manifest) == 0:
        return manifest
    manifest["_seed_sort"] = manifest["fit_seed"].fillna(10**9)
    manifest = (
        manifest.sort_values(["K", "_seed_sort", "fit_label", "fit_path"])
        .drop(columns="_seed_sort")
        .reset_index(drop=True)
    )
    return manifest


def build_task_manifest(
    fit_manifest: pd.DataFrame,
    conditions: Sequence[str],
) -> pd.DataFrame:
    rows = []
    task_id = 0
    for fit_row in fit_manifest.itertuples(index=False):
        for condition in conditions:
            rows.append(
                {
                    "task_id": task_id,
                    "fit_label": str(fit_row.fit_label),
                    "fit_seed": fit_row.fit_seed,
                    "K": int(fit_row.K),
                    "condition": str(condition),
                    "fit_path": str(fit_row.fit_path),
                    "fit_fingerprint": str(fit_row.fit_fingerprint),
                }
            )
            task_id += 1
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# CSV normalization and fit-row matching
# -----------------------------------------------------------------------------


def add_source_metadata(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    out = df.copy()
    out["_source_csv"] = str(path)
    source_seed = infer_fit_seed(path)
    out["_source_fit_seed"] = source_seed
    out["_source_fit_label"] = (
        f"seed{source_seed}" if source_seed is not None else path.parent.name
    )
    return out


def load_many_csv(patterns: Optional[Sequence[str]]) -> pd.DataFrame:
    frames = []
    for path in expand_paths(patterns):
        try:
            frame = pd.read_csv(path)
        except Exception as exc:
            print(f"WARNING: could not read {path}: {exc}", file=sys.stderr)
            continue
        frames.append(add_source_metadata(frame, path))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def alias_column(df: pd.DataFrame, target: str, aliases: Sequence[str]) -> pd.DataFrame:
    if target in df.columns:
        return df
    for alias in aliases:
        if alias in df.columns:
            return df.rename(columns={alias: target})
    return df


def normalize_ranking_df(df: pd.DataFrame) -> pd.DataFrame:
    if len(df) == 0:
        return df
    out = df.copy()
    aliases = {
        "K": ["k", "model_order"],
        "condition": ["cond"],
        "outer_rep": ["outer_split", "split", "outer_fold"],
        "direction": ["split_direction", "rank_direction"],
        "archetype": ["component", "component_id", "archetype_id"],
        "mean_accuracy": ["accuracy", "mean", "decode_accuracy"],
        "rank_mean": ["mean_rank", "rank"],
        "error_mean": ["mean_error", "error"],
    }
    for target, options in aliases.items():
        out = alias_column(out, target, options)
    required = {"condition", "K", "outer_rep", "direction", "archetype", "mean_accuracy"}
    missing = sorted(required - set(out.columns))
    if missing:
        raise ValueError(f"Ranking CSV is missing required columns: {missing}")
    return out


def normalize_run_df(df: pd.DataFrame) -> pd.DataFrame:
    if len(df) == 0:
        return df
    out = df.copy()
    aliases = {
        "K": ["k", "model_order"],
        "condition": ["cond"],
        "outer_rep": ["outer_split", "split", "outer_fold"],
        "direction": ["split_direction", "rank_direction"],
        "top_m": ["m", "topm"],
        "mean": ["accuracy", "mean_accuracy", "decode_accuracy"],
        "err": ["sem", "sem_accuracy", "stderr"],
    }
    for target, options in aliases.items():
        out = alias_column(out, target, options)
    required = {"condition", "K", "outer_rep", "direction", "top_m"}
    missing = sorted(required - set(out.columns))
    if missing:
        raise ValueError(f"Existing-run CSV is missing required columns: {missing}")
    return out


def match_seed_series(series: pd.Series, fit_seed: Optional[int]) -> pd.Series:
    if fit_seed is None:
        return pd.Series(False, index=series.index)
    normalized = series.map(normalize_seed_value)
    return normalized == int(fit_seed)


def filter_for_fit(
    df: pd.DataFrame,
    task: pd.Series,
    *,
    context: str,
) -> pd.DataFrame:
    if len(df) == 0:
        return df.copy()
    base = df.copy()
    base = base[base["K"].astype(int) == int(task["K"])].copy()
    if len(base) == 0:
        return base

    fingerprint = str(task["fit_fingerprint"])
    fit_path = Path(str(task["fit_path"])).resolve()
    fit_label = str(task["fit_label"])
    fit_seed_raw = task.get("fit_seed")
    fit_seed = None if pd.isna(fit_seed_raw) else int(fit_seed_raw)

    # Strongest identifier: exact fit fingerprint.
    if "fit_fingerprint" in base.columns and base["fit_fingerprint"].notna().any():
        selected = base[base["fit_fingerprint"].astype(str) == fingerprint].copy()
        if len(selected):
            return selected

    # Exact or basename fit path.
    for column in ["fit_path", "fit_source_file", "source_fit", "npz_path"]:
        if column not in base.columns or not base[column].notna().any():
            continue
        values = base[column].astype(str)
        selected = base[
            (values == str(fit_path))
            | (values.map(lambda x: Path(x).name) == fit_path.name)
        ].copy()
        if len(selected):
            return selected

    # Seed columns in combined ten-seed tables.
    for column in ["fit_seed", "seed", "seed_id", "rng_seed", "fit_label"]:
        if column not in base.columns or not base[column].notna().any():
            continue
        selected = base[match_seed_series(base[column], fit_seed)].copy()
        if len(selected):
            return selected
        if column == "fit_label":
            selected = base[base[column].astype(str) == fit_label].copy()
            if len(selected):
                return selected

    # Separate per-seed CSV paths.
    if "_source_fit_seed" in base.columns:
        selected = base[base["_source_fit_seed"] == fit_seed].copy()
        if len(selected):
            return selected
    if "_source_fit_label" in base.columns:
        selected = base[base["_source_fit_label"].astype(str) == fit_label].copy()
        if len(selected):
            return selected

    # Safe fallback only when there is one apparent fit source.
    source_columns = [
        c
        for c in ["fit_fingerprint", "fit_seed", "seed", "fit_label", "_source_csv"]
        if c in base.columns and base[c].notna().any()
    ]
    unique_sources = {
        tuple(str(row[c]) for c in source_columns)
        for _, row in base[source_columns].drop_duplicates().iterrows()
    }
    if len(unique_sources) <= 1:
        return base

    raise ValueError(
        f"Could not map {context} rows to {fit_label} (seed={fit_seed}, K={task['K']}). "
        "Provide a fit_fingerprint/fit_seed/fit_label column, separate seed-labeled CSV paths, "
        "or a fit_source_file column."
    )


# -----------------------------------------------------------------------------
# Held-out splitting and decoding
# -----------------------------------------------------------------------------


def outer_seed_for(base_seed: int, condition: str, outer_rep: int, K: int) -> int:
    sequence = np.random.SeedSequence(
        [
            int(base_seed),
            1,  # temporal-analysis tag used by the optimized pipeline
            int(K),
            int(CONDITION_INDEX[condition]),
            int(outer_rep),
        ]
    )
    return int(sequence.generate_state(1, dtype=np.uint32)[0])


def evaluation_seed_for(base_seed: int, K: int, outer_rep: int, direction_offset: int) -> int:
    return int(base_seed + 200000 * int(outer_rep) + 1000 * int(K) + int(direction_offset))


def split_outer_subject_indices(n_subjects: int, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(int(seed))
    indices = np.arange(n_subjects)
    rng.shuffle(indices)
    midpoint = n_subjects // 2
    return indices[:midpoint], indices[midpoint:]


def condition_subjects(fit: Dict[str, object], condition: str) -> List[dict]:
    labels = fit.get("condition_labels_str")
    if labels is None:
        raise ValueError(f"Fit lacks condition_labels_str: {fit['path']}")
    indices = np.where(np.asarray(labels).astype(str) == str(condition))[0]
    if len(indices) == 0:
        raise ValueError(f"No {condition!r} observations in fit: {fit['path']}")
    return [fit["results_subj"][int(i)] for i in indices]


def get_xval_assignments(ndata: int, nfolds: int, rng: np.random.Generator) -> np.ndarray:
    assignments = np.zeros(ndata, dtype=int)
    groupsize = int(np.ceil(ndata / nfolds))
    for fold in range(1, nfolds):
        indices = np.arange(fold * groupsize, min((fold + 1) * groupsize, ndata))
        assignments[indices] = fold
    rng.shuffle(assignments)
    return assignments


# Same tie tolerance as msaa_pipeline.decode._decode_corr with MSAA_TIE_TOL=1e-08, which
# every held-out run used. Temporal S is
# mostly exact zeros, so without it ties are broken by last-bit rounding and the numbers
# would not be comparable with the existing held-out top-m.
TIE_TOL = float(os.environ.get("MSAA_TIE_TOL", "1e-08"))


def decoder_metrics(correlations: np.ndarray) -> Dict[str, float]:
    T = correlations.shape[0]
    true_time = np.arange(T)
    if TIE_TOL > 0:
        predicted_time = np.argmax(
            correlations >= correlations.max(axis=1, keepdims=True) - TIE_TOL, axis=1
        )
    else:
        predicted_time = np.argmax(correlations, axis=1)
    diagonal = np.diag(correlations)
    return {
        "accuracy": float(np.mean(predicted_time == true_time)),
        "error": float(np.mean(np.abs(predicted_time - true_time)) / T),
        "rank": float(np.mean(np.mean(correlations <= diagonal[:, None] + TIE_TOL, axis=1))),
    }


def summarize_metric_rows(rows: List[Dict[str, float]]) -> Dict[str, float]:
    df = pd.DataFrame(rows)
    return {
        "mean": float(df["accuracy"].mean()),
        "err": float(df["accuracy"].std(ddof=1) / math.sqrt(max(len(df), 1))),
        "std": float(df["accuracy"].std(ddof=1)),
        "count": int(len(df)),
        "rank_mean": float(df["rank"].mean()),
        "error_mean": float(df["error"].mean()),
    }


def batch_decode_stack(
    representation_stack: np.ndarray,
    *,
    nfolds: int,
    nreps: int,
    seed: int,
) -> Dict[str, float]:
    if int(nfolds) != 2:
        raise ValueError("Optimized decoder currently requires nfolds=2.")
    X = np.asarray(representation_stack, dtype=float)
    N, T, F = X.shape
    if F < 2:
        return {
            "mean": float("nan"),
            "err": float("nan"),
            "std": float("nan"),
            "count": 0,
            "rank_mean": float("nan"),
            "error_mean": float("nan"),
        }

    rng = np.random.default_rng(int(seed))
    masks0 = np.zeros((nreps, N), dtype=bool)
    for repetition in range(nreps):
        masks0[repetition] = get_xval_assignments(N, 2, rng) == 0
    n0 = masks0.sum(axis=1).astype(float)
    n1 = (N - masks0.sum(axis=1)).astype(float)
    A = np.einsum("rn,ntf->rtf", masks0 / n0[:, None], X, optimize=True)
    B = np.einsum("rn,ntf->rtf", (~masks0) / n1[:, None], X, optimize=True)
    A = A - A.mean(axis=2, keepdims=True)
    B = B - B.mean(axis=2, keepdims=True)
    norm_a = np.linalg.norm(A, axis=2)
    norm_b = np.linalg.norm(B, axis=2)
    A = A / np.where(norm_a[:, :, None] > 0, norm_a[:, :, None], 1.0)
    B = B / np.where(norm_b[:, :, None] > 0, norm_b[:, :, None], 1.0)
    correlations = np.matmul(A, np.swapaxes(B, 1, 2))
    correlations = np.nan_to_num(correlations, nan=0.0, posinf=0.0, neginf=0.0)

    rows: List[Dict[str, float]] = []
    for repetition in range(nreps):
        forward = decoder_metrics(correlations[repetition])
        forward.update({"rep": repetition, "fold": 0})
        rows.append(forward)
        reverse = decoder_metrics(correlations[repetition].T)
        reverse.update({"rep": repetition, "fold": 1})
        rows.append(reverse)
    return summarize_metric_rows(rows)


def temporal_contribution_stack(subjects: Sequence[dict], archetype: int) -> np.ndarray:
    temporal = np.stack(
        [np.asarray(subject["S"], dtype=float)[int(archetype), :] for subject in subjects],
        axis=0,
    )
    spatial = np.stack(
        [np.asarray(subject["sXC"], dtype=float)[:, int(archetype)] for subject in subjects],
        axis=0,
    )
    return temporal[:, :, None] * spatial[:, None, :]


def normalize_stack_over_time(raw_stack: np.ndarray) -> np.ndarray:
    mean = raw_stack.mean(axis=1, keepdims=True)
    std = raw_stack.std(axis=1, keepdims=True)
    return (raw_stack - mean) / (std + EPS)


def full_temporal_reconstruction_stack(subjects: Sequence[dict]) -> np.ndarray:
    reconstructed = []
    for subject in subjects:
        spatial_maps = np.asarray(subject["sXC"], dtype=float)
        temporal_weights = np.asarray(subject["S"], dtype=float)
        Xhat = (spatial_maps @ temporal_weights).T
        Xhat = (Xhat - Xhat.mean(axis=0, keepdims=True)) / (
            Xhat.std(axis=0, keepdims=True) + EPS
        )
        reconstructed.append(Xhat)
    return np.stack(reconstructed, axis=0)


def evaluate_missing_m(
    subjects: Sequence[dict],
    ranked_archetypes: Sequence[int],
    missing_m: Sequence[int],
    *,
    nfolds: int,
    nreps: int,
    seed: int,
    full_summary_override: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    desired = sorted({int(m) for m in missing_m if 1 <= int(m) <= len(ranked_archetypes)})
    if not desired:
        return pd.DataFrame()

    if full_summary_override is None:
        full_summary = batch_decode_stack(
            full_temporal_reconstruction_stack(subjects),
            nfolds=nfolds,
            nreps=nreps,
            seed=seed,
        )
    else:
        full_summary = full_summary_override

    T = np.asarray(subjects[0]["S"]).shape[1]
    V = np.asarray(subjects[0]["sXC"]).shape[0]
    chance = 1.0 / float(T)
    denominator = float(full_summary["mean"]) - chance
    raw = np.zeros((len(subjects), T, V), dtype=float)
    rows = []
    max_m = max(desired)
    checkpoints = set(desired)

    for position, archetype in enumerate(ranked_archetypes[:max_m], start=1):
        raw += temporal_contribution_stack(subjects, int(archetype))
        if position not in checkpoints:
            continue
        summary = batch_decode_stack(
            normalize_stack_over_time(raw),
            nfolds=nfolds,
            nreps=nreps,
            seed=seed,
        )
        retention = (
            (float(summary["mean"]) - chance) / denominator
            if abs(denominator) > 1e-12
            else float("nan")
        )
        rows.append(
            {
                "top_m": int(position),
                "selected_archetypes": str([int(x) for x in ranked_archetypes[:position]]),
                **summary,
                "full_mean_same_eval": float(full_summary["mean"]),
                "full_err_same_eval": float(full_summary.get("err", float("nan"))),
                "retention_above_chance": retention,
                "chance": chance,
                "result_source": "computed",
            }
        )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Ranking/run extraction for one task
# -----------------------------------------------------------------------------


def task_ranking_rows(ranking_df: pd.DataFrame, task: pd.Series) -> pd.DataFrame:
    fit_rows = filter_for_fit(ranking_df, task, context="ranking")
    current = fit_rows[
        (fit_rows["condition"].astype(str) == str(task["condition"]))
        & (fit_rows["K"].astype(int) == int(task["K"]))
    ].copy()
    return current


def ranking_for_direction(
    ranking_rows: pd.DataFrame,
    *,
    outer_rep: int,
    direction: str,
    K: int,
) -> Tuple[List[int], Optional[int]]:
    current = ranking_rows[
        (ranking_rows["outer_rep"].astype(int) == int(outer_rep))
        & (ranking_rows["direction"].astype(str) == str(direction))
    ].copy()
    if len(current) == 0:
        raise ValueError(f"No ranking rows for outer_rep={outer_rep}, direction={direction}")

    # Collapse accidental duplicates conservatively.
    agg_spec = {"mean_accuracy": "mean"}
    if "rank_mean" in current.columns:
        agg_spec["rank_mean"] = "mean"
    if "error_mean" in current.columns:
        agg_spec["error_mean"] = "mean"
    if "outer_seed" in current.columns:
        agg_spec["outer_seed"] = "first"
    current = current.groupby("archetype", as_index=False).agg(agg_spec)
    if current["archetype"].nunique() != int(K):
        raise ValueError(
            f"Incomplete ranking for K={K}, outer_rep={outer_rep}, direction={direction}: "
            f"{current['archetype'].nunique()}/{K} archetypes"
        )

    sort_columns = ["mean_accuracy"]
    ascending = [False]
    if "rank_mean" in current.columns:
        sort_columns.append("rank_mean")
        ascending.append(False)
    # No error_mean tie-break: msaa_pipeline.heldout.rank_archetypes_on_subjects
    # sorts by mean_accuracy, rank_mean, archetype, and the held-out top-m used that order.
    sort_columns.append("archetype")
    ascending.append(True)
    ranked = (
        current.sort_values(sort_columns, ascending=ascending)["archetype"]
        .astype(int)
        .tolist()
    )
    outer_seed = None
    if "outer_seed" in current.columns and current["outer_seed"].notna().any():
        unique = current["outer_seed"].dropna().astype(np.uint64).unique()
        if len(unique) == 1:
            outer_seed = int(unique[0])
    return ranked, outer_seed


def task_existing_run_rows(run_df: pd.DataFrame, task: pd.Series) -> pd.DataFrame:
    if len(run_df) == 0:
        return run_df.copy()
    fit_rows = filter_for_fit(run_df, task, context="existing-run")
    return fit_rows[
        (fit_rows["condition"].astype(str) == str(task["condition"]))
        & (fit_rows["K"].astype(int) == int(task["K"]))
    ].copy()


def existing_direction_rows(
    existing_rows: pd.DataFrame,
    outer_rep: int,
    direction: str,
) -> pd.DataFrame:
    if len(existing_rows) == 0:
        return existing_rows.copy()
    return existing_rows[
        (existing_rows["outer_rep"].astype(int) == int(outer_rep))
        & (existing_rows["direction"].astype(str) == str(direction))
    ].copy()


def summary_from_existing_row(row: pd.Series) -> Dict[str, float]:
    def get(name: str, default: float = float("nan")) -> float:
        return float(row[name]) if name in row.index and pd.notna(row[name]) else default

    return {
        "mean": get("mean"),
        "err": get("err"),
        "std": get("std"),
        "count": int(get("count", 0)),
        "rank_mean": get("rank_mean"),
        "error_mean": get("error_mean"),
    }


def full_summary_from_existing(current: pd.DataFrame) -> Optional[Dict[str, float]]:
    if len(current) == 0 or "full_mean_same_eval" not in current.columns:
        return None
    valid = current[current["full_mean_same_eval"].notna()]
    if len(valid) == 0:
        return None
    row = valid.iloc[0]
    return {
        "mean": float(row["full_mean_same_eval"]),
        "err": float(row["full_err_same_eval"]) if "full_err_same_eval" in row.index and pd.notna(row["full_err_same_eval"]) else float("nan"),
        "std": float("nan"),
        "count": 0,
        "rank_mean": float("nan"),
        "error_mean": float("nan"),
    }


# -----------------------------------------------------------------------------
# Subcommands
# -----------------------------------------------------------------------------


def command_build_manifests(args: argparse.Namespace) -> None:
    roots = [Path(x) for x in args.fit_root]
    k_values = parse_csv_ints(args.k_values)
    conditions = parse_csv_strings(args.conditions)
    invalid = sorted(set(conditions) - set(ALL_CONDITIONS))
    if invalid:
        raise ValueError(f"Unknown conditions: {invalid}")
    fit_manifest = discover_fit_manifest(roots, k_values)
    if len(fit_manifest) == 0:
        raise RuntimeError("No temporal fit files were discovered.")
    task_manifest = build_task_manifest(fit_manifest, conditions)
    atomic_write_csv(fit_manifest, Path(args.fit_manifest))
    atomic_write_csv(task_manifest, Path(args.task_manifest))
    print(f"Wrote {len(fit_manifest)} fits to {args.fit_manifest}")
    print(f"Wrote {len(task_manifest)} tasks to {args.task_manifest}")
    coverage = fit_manifest.groupby("K").agg(
        n_fits=("fit_label", "size"),
        labels=("fit_label", lambda x: ",".join(map(str, x))),
    )
    print("\nFit coverage by K:")
    print(coverage.to_string())


def command_audit(args: argparse.Namespace) -> None:
    tasks = pd.read_csv(args.task_manifest)
    ranking_df = normalize_ranking_df(load_many_csv(args.ranking_csv))
    rows = []
    for task in tasks.itertuples(index=False):
        task_series = pd.Series(task._asdict())
        try:
            current = task_ranking_rows(ranking_df, task_series)
            n_outer = current["outer_rep"].nunique() if len(current) else 0
            n_directions = current["direction"].nunique() if len(current) else 0
            group_counts = (
                current.groupby(["outer_rep", "direction"])["archetype"].nunique()
                if len(current)
                else pd.Series(dtype=int)
            )
            complete_groups = int((group_counts == int(task.K)).sum())
            expected_groups = int(args.outer_reps) * 2
            error = ""
        except Exception as exc:
            n_outer = n_directions = complete_groups = 0
            expected_groups = int(args.outer_reps) * 2
            error = repr(exc)
        rows.append(
            {
                "task_id": int(task.task_id),
                "fit_label": task.fit_label,
                "fit_seed": task.fit_seed,
                "K": int(task.K),
                "condition": task.condition,
                "n_outer_reps": n_outer,
                "n_directions": n_directions,
                "complete_direction_groups": complete_groups,
                "expected_direction_groups": expected_groups,
                "complete": complete_groups == expected_groups,
                "error": error,
            }
        )
    audit = pd.DataFrame(rows)
    if args.output:
        atomic_write_csv(audit, Path(args.output))
    print(audit.groupby(["K", "condition"])["complete"].agg(["sum", "count"]).to_string())
    failures = audit[~audit["complete"]]
    if len(failures):
        print("\nIncomplete tasks:")
        print(failures.to_string(index=False))
        raise SystemExit(2)
    print("\nAll task rankings are complete.")


def result_file_for_task(output_dir: Path, task: pd.Series) -> Path:
    safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(task["fit_label"]))
    return output_dir / "tasks" / f"task_{int(task['task_id']):04d}_{safe_label}_K{int(task['K'])}_{task['condition']}.csv"


def command_run_task(args: argparse.Namespace) -> None:
    tasks = pd.read_csv(args.task_manifest)
    matches = tasks[tasks["task_id"].astype(int) == int(args.task_id)]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one task_id={args.task_id}; found {len(matches)}")
    task = matches.iloc[0]
    fractions = parse_csv_floats(args.fractions)
    if not fractions:
        raise ValueError("No fractions were requested.")
    output_dir = Path(args.output_dir)
    result_path = result_file_for_task(output_dir, task)
    checkpoint = load_csv(result_path)

    ranking_df = normalize_ranking_df(load_many_csv(args.ranking_csv))
    ranking_rows = task_ranking_rows(ranking_df, task)
    existing_run_df = normalize_run_df(load_many_csv(args.existing_run_csv)) if args.existing_run_csv else pd.DataFrame()
    existing_rows = task_existing_run_rows(existing_run_df, task) if len(existing_run_df) else pd.DataFrame()

    fit = load_fit(Path(task["fit_path"]))
    if int(fit["K"]) != int(task["K"]):
        raise ValueError(f"Task K={task['K']} but fit contains K={fit['K']}")
    subjects = condition_subjects(fit, str(task["condition"]))
    K = int(task["K"])
    mapping = pd.DataFrame(
        [
            {
                "fraction_requested": float(fraction),
                "top_m": int(np.clip(np.round(float(fraction) * K), 1, K)),
            }
            for fraction in fractions
        ]
    )
    mapping["actual_fraction"] = mapping["top_m"] / float(K)
    mapping = mapping.drop_duplicates(["fraction_requested", "top_m"]).reset_index(drop=True)

    output_keys = [
        "fit_fingerprint",
        "condition",
        "K",
        "outer_rep",
        "direction",
        "fraction_requested",
        "top_m",
    ]

    for outer_rep in range(int(args.outer_reps)):
        for direction_offset, direction in enumerate(["A_to_B", "B_to_A"]):
            expected = {
                (float(row.fraction_requested), int(row.top_m))
                for row in mapping.itertuples(index=False)
            }
            existing_checkpoint = set()
            if len(checkpoint):
                current_checkpoint = checkpoint[
                    (checkpoint["outer_rep"].astype(int) == outer_rep)
                    & (checkpoint["direction"].astype(str) == direction)
                ]
                existing_checkpoint = {
                    (float(row.fraction_requested), int(row.top_m))
                    for row in current_checkpoint.itertuples(index=False)
                }
            missing = expected - existing_checkpoint
            if not missing and args.resume:
                continue

            ranked, ranking_outer_seed = ranking_for_direction(
                ranking_rows,
                outer_rep=outer_rep,
                direction=direction,
                K=K,
            )
            outer_seed = (
                ranking_outer_seed
                if ranking_outer_seed is not None
                else outer_seed_for(int(args.base_seed), str(task["condition"]), outer_rep, K)
            )
            group_a, group_b = split_outer_subject_indices(len(subjects), outer_seed)
            evaluation_idx = group_b if direction == "A_to_B" else group_a
            evaluation_subjects = [subjects[int(i)] for i in evaluation_idx]
            evaluation_seed = evaluation_seed_for(int(args.base_seed), K, outer_rep, direction_offset)

            direction_existing = existing_direction_rows(existing_rows, outer_rep, direction)
            full_override = full_summary_from_existing(direction_existing)
            rows = []
            missing_top_m = sorted({int(top_m) for _, top_m in missing if int(top_m) < K})
            computed = evaluate_missing_m(
                evaluation_subjects,
                ranked,
                missing_top_m,
                nfolds=int(args.nfolds),
                nreps=int(args.nreps),
                seed=evaluation_seed,
                full_summary_override=full_override,
            ) if missing_top_m else pd.DataFrame()
            computed_by_m = {
                int(row.top_m): pd.Series(row._asdict())
                for row in computed.itertuples(index=False)
            }

            for map_row in mapping.itertuples(index=False):
                pair = (float(map_row.fraction_requested), int(map_row.top_m))
                if pair not in missing:
                    continue
                top_m = int(map_row.top_m)
                exact = direction_existing[
                    direction_existing["top_m"].astype(int) == top_m
                ] if len(direction_existing) else pd.DataFrame()

                if len(exact):
                    source_row = exact.iloc[0]
                    summary = summary_from_existing_row(source_row)
                    full_summary = full_summary_from_existing(direction_existing)
                    chance = float(source_row["chance"]) if "chance" in source_row.index and pd.notna(source_row["chance"]) else 1.0 / np.asarray(evaluation_subjects[0]["S"]).shape[1]
                    if full_summary is None:
                        full_summary = batch_decode_stack(
                            full_temporal_reconstruction_stack(evaluation_subjects),
                            nfolds=int(args.nfolds),
                            nreps=int(args.nreps),
                            seed=evaluation_seed,
                        )
                    retention = (
                        float(source_row["retention_above_chance"])
                        if "retention_above_chance" in source_row.index and pd.notna(source_row["retention_above_chance"])
                        else (summary["mean"] - chance) / (full_summary["mean"] - chance)
                    )
                    selected = (
                        source_row["selected_archetypes"]
                        if "selected_archetypes" in source_row.index and pd.notna(source_row["selected_archetypes"])
                        else str(ranked[:top_m])
                    )
                    result_source = "reused_exact"
                elif top_m == K:
                    full_summary = full_override
                    if full_summary is None:
                        full_summary = batch_decode_stack(
                            full_temporal_reconstruction_stack(evaluation_subjects),
                            nfolds=int(args.nfolds),
                            nreps=int(args.nreps),
                            seed=evaluation_seed,
                        )
                    summary = dict(full_summary)
                    chance = 1.0 / np.asarray(evaluation_subjects[0]["S"]).shape[1]
                    retention = 1.0
                    selected = str([int(x) for x in ranked])
                    result_source = "reused_or_computed_full"
                else:
                    if top_m not in computed_by_m:
                        raise RuntimeError(f"Missing computed K={K}, m={top_m}")
                    source_row = computed_by_m[top_m]
                    summary = {
                        key: float(source_row[key]) if key in source_row.index and pd.notna(source_row[key]) else float("nan")
                        for key in ["mean", "err", "std", "rank_mean", "error_mean"]
                    }
                    summary["count"] = int(source_row["count"])
                    full_summary = {
                        "mean": float(source_row["full_mean_same_eval"]),
                        "err": float(source_row["full_err_same_eval"]),
                    }
                    chance = float(source_row["chance"])
                    retention = float(source_row["retention_above_chance"])
                    selected = str(source_row["selected_archetypes"])
                    result_source = "computed"

                rows.append(
                    {
                        "task_id": int(task["task_id"]),
                        "fit_label": str(task["fit_label"]),
                        "fit_seed": task["fit_seed"],
                        "fit_fingerprint": str(task["fit_fingerprint"]),
                        "fit_path": str(task["fit_path"]),
                        "condition": str(task["condition"]),
                        "K": K,
                        "outer_rep": outer_rep,
                        "direction": direction,
                        "outer_seed": int(outer_seed),
                        "evaluation_seed": int(evaluation_seed),
                        "fraction_requested": float(map_row.fraction_requested),
                        "top_m": top_m,
                        "actual_fraction": float(map_row.actual_fraction),
                        "selected_archetypes": selected,
                        **summary,
                        "full_mean_same_eval": float(full_summary["mean"]),
                        "full_err_same_eval": float(full_summary.get("err", float("nan"))),
                        "retention_above_chance": float(retention),
                        "chance": float(chance),
                        "n_evaluation_subjects": len(evaluation_subjects),
                        "result_source": result_source,
                        "script_version": VERSION,
                    }
                )

            if rows:
                new = pd.DataFrame(rows)
                checkpoint = pd.concat([checkpoint, new], ignore_index=True, sort=False) if len(checkpoint) else new
                checkpoint = checkpoint.drop_duplicates(output_keys, keep="last").sort_values(
                    ["outer_rep", "direction", "fraction_requested", "top_m"]
                ).reset_index(drop=True)
                atomic_write_csv(checkpoint, result_path)
                print(
                    f"task={int(task['task_id'])} {task['fit_label']} K={K} {task['condition']} "
                    f"outer={outer_rep + 1}/{args.outer_reps} {direction} "
                    f"completed fractions={sorted(float(x[0]) for x in missing)}",
                    flush=True,
                )

    expected_rows = int(args.outer_reps) * 2 * len(mapping)
    if len(checkpoint) != expected_rows:
        raise RuntimeError(
            f"Task completed with {len(checkpoint)} rows; expected {expected_rows}. Result: {result_path}"
        )
    print(f"COMPLETE: {result_path} ({len(checkpoint)} rows)")


def command_aggregate(args: argparse.Namespace) -> None:
    tasks = pd.read_csv(args.task_manifest)
    output_dir = Path(args.output_dir)
    frames = []
    coverage_rows = []
    fractions = parse_csv_floats(args.fractions)
    expected_per_task = int(args.outer_reps) * 2 * len(
        pd.DataFrame({"fraction": fractions}).drop_duplicates()
    )

    for task in tasks.itertuples(index=False):
        task_series = pd.Series(task._asdict())
        path = result_file_for_task(output_dir, task_series)
        if not path.exists():
            coverage_rows.append({**task._asdict(), "result_path": str(path), "rows": 0, "complete": False})
            continue
        df = pd.read_csv(path)
        frames.append(df)
        coverage_rows.append(
            {
                **task._asdict(),
                "result_path": str(path),
                "rows": len(df),
                "complete": len(df) == expected_per_task,
            }
        )

    coverage = pd.DataFrame(coverage_rows)
    aggregate_dir = output_dir / "aggregate"
    aggregate_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_csv(coverage, aggregate_dir / "matched_fraction_coverage.csv")
    if not frames:
        raise RuntimeError("No task outputs were found.")
    runs = pd.concat(frames, ignore_index=True, sort=False)
    key_cols = [
        "fit_fingerprint",
        "condition",
        "K",
        "outer_rep",
        "direction",
        "fraction_requested",
        "top_m",
    ]
    runs = runs.drop_duplicates(key_cols, keep="last").reset_index(drop=True)
    atomic_write_csv(runs, aggregate_dir / "matched_fraction_all_runs.csv")

    outer = (
        runs.groupby(
            [
                "fit_label",
                "fit_seed",
                "fit_fingerprint",
                "condition",
                "K",
                "fraction_requested",
                "top_m",
                "actual_fraction",
                "outer_rep",
            ],
            as_index=False,
            dropna=False,
        )
        .agg(
            accuracy=("mean", "mean"),
            retention=("retention_above_chance", "mean"),
            full_accuracy=("full_mean_same_eval", "mean"),
        )
    )
    atomic_write_csv(outer, aggregate_dir / "matched_fraction_outer_split.csv")

    seed_summary = (
        outer.groupby(
            [
                "fit_label",
                "fit_seed",
                "fit_fingerprint",
                "condition",
                "K",
                "fraction_requested",
                "top_m",
                "actual_fraction",
            ],
            as_index=False,
            dropna=False,
        )
        .agg(
            mean_accuracy=("accuracy", "mean"),
            outer_sem_accuracy=("accuracy", sem),
            mean_retention=("retention", "mean"),
            outer_sem_retention=("retention", sem),
            mean_full_accuracy=("full_accuracy", "mean"),
            n_outer=("outer_rep", "nunique"),
        )
    )
    atomic_write_csv(seed_summary, aggregate_dir / "matched_fraction_seed_summary.csv")

    summary = (
        seed_summary.groupby(
            ["condition", "K", "fraction_requested", "top_m", "actual_fraction"],
            as_index=False,
        )
        .agg(
            mean_accuracy=("mean_accuracy", "mean"),
            fit_seed_sd_accuracy=("mean_accuracy", "std"),
            fit_seed_sem_accuracy=("mean_accuracy", sem),
            mean_outer_sem_accuracy=("outer_sem_accuracy", "mean"),
            mean_retention=("mean_retention", "mean"),
            fit_seed_sd_retention=("mean_retention", "std"),
            fit_seed_sem_retention=("mean_retention", sem),
            mean_outer_sem_retention=("outer_sem_retention", "mean"),
            mean_full_accuracy=("mean_full_accuracy", "mean"),
            n_fit_seeds=("fit_label", "nunique"),
        )
    )
    atomic_write_csv(summary, aggregate_dir / "matched_fraction_summary.csv")

    peak_rows = []
    for (condition, fraction), current in summary.groupby(["condition", "fraction_requested"]):
        current = current.dropna(subset=["mean_retention"]).copy()
        if len(current) == 0:
            continue
        best = current.loc[current["mean_retention"].idxmax()]
        peak_rows.append(
            {
                "condition": condition,
                "fraction_requested": float(fraction),
                "best_K": int(best["K"]),
                "best_top_m": int(best["top_m"]),
                "best_actual_fraction": float(best["actual_fraction"]),
                "best_mean_retention": float(best["mean_retention"]),
                "n_fit_seeds": int(best["n_fit_seeds"]),
            }
        )
    atomic_write_csv(pd.DataFrame(peak_rows), aggregate_dir / "matched_fraction_peak_by_condition.csv")

    if not args.no_plots:
        import matplotlib.pyplot as plt

        for condition in sorted(summary["condition"].astype(str).unique()):
            current_condition = summary[summary["condition"].astype(str) == condition]
            plt.figure(figsize=(8.5, 5.2))
            for fraction in sorted(current_condition["fraction_requested"].unique()):
                current = current_condition[
                    np.isclose(current_condition["fraction_requested"], float(fraction))
                ].sort_values("K")
                yerr = current["fit_seed_sd_retention"].fillna(current["mean_outer_sem_retention"]).fillna(0)
                plt.plot(current["K"], current["mean_retention"], marker="o", label=f"{fraction:.0%} retained")
                plt.fill_between(
                    current["K"],
                    current["mean_retention"] - yerr,
                    current["mean_retention"] + yerr,
                    alpha=0.15,
                )
            plt.axhline(1.0, linestyle="--", linewidth=1)
            plt.xlabel("Temporal model order (K)")
            plt.ylabel("Fraction of above-chance full-model performance retained")
            plt.title(f"Temporal MS-AA matched retained fractions | {condition}")
            plt.legend()
            plt.tight_layout()
            plt.savefig(aggregate_dir / f"matched_fraction_{condition}.png", dpi=300, bbox_inches="tight")
            plt.close()

        for fraction in sorted(summary["fraction_requested"].unique()):
            plt.figure(figsize=(8.5, 5.2))
            for condition in sorted(summary["condition"].astype(str).unique()):
                current = summary[
                    (summary["condition"].astype(str) == condition)
                    & np.isclose(summary["fraction_requested"], float(fraction))
                ].sort_values("K")
                yerr = current["fit_seed_sd_retention"].fillna(current["mean_outer_sem_retention"]).fillna(0)
                plt.plot(current["K"], current["mean_retention"], marker="o", label=condition)
                plt.fill_between(
                    current["K"],
                    current["mean_retention"] - yerr,
                    current["mean_retention"] + yerr,
                    alpha=0.15,
                )
            plt.axhline(1.0, linestyle="--", linewidth=1)
            plt.xlabel("Temporal model order (K)")
            plt.ylabel("Fraction of above-chance full-model performance retained")
            plt.title(f"Temporal MS-AA | {fraction:.0%} retained")
            plt.legend()
            plt.tight_layout()
            plt.savefig(aggregate_dir / f"matched_fraction_{fraction:.2f}.png", dpi=300, bbox_inches="tight")
            plt.close()

    complete_count = int(coverage["complete"].sum())
    print(f"Coverage: {complete_count}/{len(coverage)} tasks complete")
    print(f"Wrote aggregate outputs to {aggregate_dir}")
    if complete_count != len(coverage) and args.require_complete:
        raise SystemExit(2)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build-manifests", help="Discover fits and create fit/task manifests")
    build.add_argument("--fit-root", action="append", required=True, help="Fit root; repeat for multiple roots")
    build.add_argument("--k-values", default=",".join(map(str, DEFAULT_K_VALUES)))
    build.add_argument("--conditions", default=",".join(DEFAULT_CONDITIONS))
    build.add_argument("--fit-manifest", required=True)
    build.add_argument("--task-manifest", required=True)
    build.set_defaults(func=command_build_manifests)

    audit = subparsers.add_parser("audit", help="Audit ranking coverage for every task")
    audit.add_argument("--task-manifest", required=True)
    audit.add_argument("--ranking-csv", action="append", required=True, help="CSV path or glob; repeat as needed")
    audit.add_argument("--outer-reps", type=int, default=20)
    audit.add_argument("--output")
    audit.set_defaults(func=command_audit)

    run = subparsers.add_parser("run-task", help="Run/resume one task (one row of the task manifest)")
    run.add_argument("--task-manifest", required=True)
    run.add_argument("--task-id", type=int, required=True)
    run.add_argument("--ranking-csv", action="append", required=True, help="CSV path or glob; repeat as needed")
    run.add_argument("--existing-run-csv", action="append", help="Optional existing held-out run CSV path/glob. Do not use with several fit seeds: every seed's run file has the same name, so rows would be matched to the wrong seed.")
    run.add_argument("--output-dir", required=True)
    run.add_argument("--fractions", default=",".join(map(str, DEFAULT_FRACTIONS)))
    run.add_argument("--outer-reps", type=int, default=20)
    run.add_argument("--nfolds", type=int, default=2)
    run.add_argument("--nreps", type=int, default=20)
    run.add_argument("--base-seed", type=int, default=42)
    run.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    run.set_defaults(func=command_run_task)

    aggregate = subparsers.add_parser("aggregate", help="Aggregate completed task CSVs")
    aggregate.add_argument("--task-manifest", required=True)
    aggregate.add_argument("--output-dir", required=True)
    aggregate.add_argument("--fractions", default=",".join(map(str, DEFAULT_FRACTIONS)))
    aggregate.add_argument("--outer-reps", type=int, default=20)
    aggregate.add_argument("--no-plots", action="store_true")
    aggregate.add_argument("--require-complete", action="store_true")
    aggregate.set_defaults(func=command_aggregate)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
