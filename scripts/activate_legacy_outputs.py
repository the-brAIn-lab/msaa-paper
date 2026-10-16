#!/usr/bin/env python3
"""Expose portable legacy outputs at historical top-level paths via relative links."""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def manifest_roots(path: Path) -> list[str]:
    if not path.is_file():
        raise SystemExit(
            f"Manifest not found: {path}\n"
            "Run `make import-legacy OLD_REPO=/path/to/original/MS_AA` first."
        )
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows or "root" not in rows[0]:
        raise SystemExit(f"Manifest has no root column: {path}")
    return sorted({row["root"] for row in rows})


def same_target(link: Path, target: Path) -> bool:
    if not link.is_symlink():
        return False
    try:
        return link.resolve() == target.resolve()
    except OSError:
        return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--store", type=Path, default=Path("reference_outputs/legacy_full"))
    ap.add_argument("--manifest", type=Path, default=Path("configs/legacy_outputs_manifest.csv"))
    ap.add_argument("--deactivate", action="store_true")
    args = ap.parse_args()

    repo = args.repo_root.resolve()
    store = (repo / args.store).resolve() if not args.store.is_absolute() else args.store.resolve()
    manifest = repo / args.manifest if not args.manifest.is_absolute() else args.manifest
    roots = manifest_roots(manifest)

    if args.deactivate:
        removed = 0
        for root in roots:
            dest = repo / root
            target = store / root
            if dest.is_symlink() and same_target(dest, target):
                dest.unlink()
                print(f"Removed: {dest}")
                removed += 1
            elif dest.exists() or dest.is_symlink():
                print(f"Left untouched (not a matching legacy link): {dest}")
        print(f"Deactivated {removed} legacy output link(s).")
        return

    if not store.is_dir():
        raise SystemExit(f"Legacy store not found: {store}")

    created = 0
    for root in roots:
        target = store / root
        dest = repo / root
        if not target.exists():
            raise SystemExit(f"Manifest root missing from portable store: {target}")
        if dest.exists() or dest.is_symlink():
            if same_target(dest, target):
                print(f"Already active: {dest}")
                continue
            raise SystemExit(
                f"Refusing to overwrite existing path: {dest}\n"
                "Move/remove it explicitly, or run with a clean repo root."
            )
        relative_target = os.path.relpath(target, start=dest.parent)
        os.symlink(relative_target, dest, target_is_directory=target.is_dir())
        print(f"Activated: {dest} -> {relative_target}")
        created += 1

    print(f"Activated {created} legacy output link(s).")
    print("These links are repo-relative and remain valid when the whole folder is moved.")


if __name__ == "__main__":
    main()
