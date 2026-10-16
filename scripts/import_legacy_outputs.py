#!/usr/bin/env python3
"""Copy expensive original-analysis outputs into a portable ignored store.

The copy is resumable. Files already present with identical SHA-256 hashes are
left in place; changed/incomplete files are replaced. The resulting tracked
manifest contains portable paths, byte sizes, and SHA-256 hashes, never the
machine-specific path to the original repository.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
from pathlib import Path

CHUNK = 8 * 1024 * 1024


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_roots(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    required = {"root", "required", "purpose"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(f"Invalid roots config: {path}; expected columns {sorted(required)}")
    return rows


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_and_hash(src: Path, dst: Path) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".partial")
    h = hashlib.sha256()
    with src.open("rb") as fin, tmp.open("wb") as fout:
        for chunk in iter(lambda: fin.read(CHUNK), b""):
            fout.write(chunk)
            h.update(chunk)
    shutil.copystat(src, tmp, follow_symlinks=True)
    os.replace(tmp, dst)
    return h.hexdigest()


def iter_files(root: Path):
    if root.is_file():
        yield root
        return
    for p in sorted(root.rglob("*")):
        if p.is_file():
            yield p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("old_repo", type=Path)
    ap.add_argument(
        "--destination",
        type=Path,
        default=Path("reference_outputs/legacy_full"),
    )
    ap.add_argument(
        "--roots-config",
        type=Path,
        default=Path("configs/legacy_output_roots.csv"),
    )
    ap.add_argument(
        "--manifest",
        type=Path,
        default=Path("configs/legacy_outputs_manifest.csv"),
    )
    args = ap.parse_args()

    old_repo = args.old_repo.expanduser().resolve()
    destination = args.destination.resolve()
    if not old_repo.is_dir():
        raise SystemExit(f"Original repository not found: {old_repo}")

    roots = load_roots(args.roots_config)
    missing_required = [
        r["root"] for r in roots
        if truthy(r["required"]) and not (old_repo / r["root"]).exists()
    ]
    if missing_required:
        raise SystemExit(
            "Required legacy output roots are missing from the original repo:\n  "
            + "\n  ".join(missing_required)
        )

    destination.mkdir(parents=True, exist_ok=True)
    records = []
    copied = reused = 0
    copied_bytes = 0

    for spec in roots:
        root_name = spec["root"]
        src_root = old_repo / root_name
        if not src_root.exists():
            print(f"Optional root absent; skipping: {root_name}")
            continue

        print(f"\n[{root_name}]")
        for src in iter_files(src_root):
            rel_within = src.relative_to(src_root) if src_root.is_dir() else Path(src.name)
            manifest_rel = Path(root_name) / rel_within
            dst = destination / manifest_rel
            size = src.stat().st_size

            digest = None
            if dst.is_file() and dst.stat().st_size == size:
                src_hash = sha256_file(src)
                dst_hash = sha256_file(dst)
                if src_hash == dst_hash:
                    digest = dst_hash
                    reused += 1

            if digest is None:
                digest = copy_and_hash(src, dst)
                copied += 1
                copied_bytes += size

            records.append({
                "root": root_name,
                "relative_path": manifest_rel.as_posix(),
                "size_bytes": size,
                "sha256": digest,
            })

            n = copied + reused
            if n % 100 == 0:
                print(f"  processed {n} files ({copied} copied, {reused} reused)")

    records.sort(key=lambda r: r["relative_path"])
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    tmp_manifest = args.manifest.with_suffix(args.manifest.suffix + ".partial")
    with tmp_manifest.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["root", "relative_path", "size_bytes", "sha256"]
        )
        writer.writeheader()
        writer.writerows(records)
    os.replace(tmp_manifest, args.manifest)

    total_bytes = sum(int(r["size_bytes"]) for r in records)
    print("\nImport complete.")
    print(f"  files in manifest: {len(records)}")
    print(f"  total stored bytes: {total_bytes:,}")
    print(f"  copied this run: {copied} files / {copied_bytes:,} bytes")
    print(f"  reused this run: {reused} files")
    print(f"  store: {destination}")
    print(f"  manifest: {args.manifest}")
    print("\nNext: make validate-legacy && make activate-legacy")


if __name__ == "__main__":
    main()
