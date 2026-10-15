from __future__ import annotations
import json
from pathlib import Path


def load_profiles(repo_root: Path) -> dict:
    path = repo_root / "configs" / "profiles.json"
    return json.loads(path.read_text())


def load_profile(repo_root: Path, name: str) -> dict:
    profiles = load_profiles(repo_root)
    if name not in profiles:
        raise KeyError(f"Unknown profile {name!r}. Choose from {sorted(profiles)}")
    return profiles[name]


def parse_job_key(key: str) -> tuple[str, str]:
    try:
        analysis_type, fit_scope = key.split("_", 1)
    except ValueError as e:
        raise ValueError(f"Invalid job key {key!r}; expected e.g. spatial_across") from e
    if analysis_type not in {"spatial", "temporal"}:
        raise ValueError(key)
    if fit_scope not in {"within", "across"}:
        raise ValueError(key)
    return analysis_type, fit_scope
