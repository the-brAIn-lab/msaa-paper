#!/usr/bin/env python3
"""Validate the portable legacy store against its tracked manifest."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

CHUNK = 8 * 1024 * 1024


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--store", type=Path, default=Path("reference_outputs/legacy_full")
    )
    ap.add_argument(
        "--manifest", type=Path, default=Path("configs/legacy_outputs_manifest.csv")
    )
    ap.add_argument("--mode", choices=["quick", "hash"], default="quick")
    args = ap.parse_args()

    if not args.manifest.is_file():
        raise SystemExit(
            f"Manifest not found: {args.manifest}\n"
            "Run `make import-legacy OLD_REPO=/path/to/original/MS_AA` first."
        )
    if not args.store.is_dir():
        raise SystemExit(f"Portable legacy store not found: {args.store}")

    with args.manifest.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    required = {"root", "relative_path", "size_bytes", "sha256"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(f"Manifest is invalid: {args.manifest}")

    failures = []
    expected = set()
    for i, row in enumerate(rows, start=1):
        rel = Path(row["relative_path"])
        expected.add(rel.as_posix())
        path = args.store / rel
        if not path.is_file():
            failures.append((rel.as_posix(), "missing"))
            continue
        actual_size = path.stat().st_size
        expected_size = int(row["size_bytes"])
        if actual_size != expected_size:
            failures.append((rel.as_posix(), f"size {actual_size} != {expected_size}"))
            continue
        if args.mode == "hash":
            digest = sha256_file(path)
            if digest != row["sha256"]:
                failures.append((rel.as_posix(), "SHA-256 mismatch"))
        if i % 100 == 0:
            print(f"Checked {i}/{len(rows)} files...")

    extras = []
    for p in args.store.rglob("*"):
        if p.is_file() and p.relative_to(args.store).as_posix() not in expected:
            extras.append(p.relative_to(args.store).as_posix())

    if failures:
        print("Validation FAILED:")
        for rel, reason in failures[:50]:
            print(f"  {rel}: {reason}")
        if len(failures) > 50:
            print(f"  ... {len(failures)-50} additional failures")
        raise SystemExit(2)

    print(f"Validation passed: {len(rows)} manifest files ({args.mode} mode).")
    if extras:
        print(f"WARNING: {len(extras)} unmanifested file(s) are also present in the store.")
        for rel in extras[:20]:
            print(f"  extra: {rel}")


if __name__ == "__main__":
    main()
