#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

DEFAULT_STORE = Path('reference_outputs/legacy_full')
DEFAULT_MANIFEST = Path('configs/legacy_outputs_manifest.csv')
DEFAULT_GITIGNORE = Path('.gitignore')
FIELDS = ['root', 'relative_path', 'size_bytes', 'sha256']
CHUNK = 8 * 1024 * 1024


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(CHUNK), b''):
            h.update(chunk)
    return h.hexdigest()


def iter_files(path: Path):
    if path.is_file():
        yield path
        return
    for p in sorted(path.rglob('*')):
        if p.is_file() and not p.is_symlink():
            yield p


def normalize_requested(store: Path, raw: str) -> tuple[Path, str]:
    p = Path(raw).expanduser()
    if p.is_absolute():
        resolved = p.resolve()
    else:
        resolved = (store / p).resolve()
    try:
        rel = resolved.relative_to(store.resolve())
    except ValueError as exc:
        raise SystemExit(
            f'Refusing path outside legacy store:\n  {resolved}\nStore is:\n  {store.resolve()}'
        ) from exc
    return resolved, rel.as_posix()


def requested_matches(row_rel: str, requested_rel: str, requested_path: Path) -> bool:
    if requested_path.is_file():
        return row_rel == requested_rel
    prefix = requested_rel.rstrip('/') + '/'
    return row_rel == requested_rel or row_rel.startswith(prefix)


def ensure_gitignore(roots: set[str], gitignore: Path) -> None:
    lines = gitignore.read_text(encoding='utf-8').splitlines() if gitignore.exists() else []
    normalized = {line.strip() for line in lines}
    to_add = []
    for root in sorted(roots):
        entry = f'/{root}'
        if entry not in normalized and root not in normalized:
            to_add.append(entry)
    if not to_add:
        print('.gitignore already covers registered root(s).')
        return
    with gitignore.open('a', encoding='utf-8') as f:
        if lines and lines[-1].strip():
            f.write('\n')
        f.write('# Activated portable legacy-output links\n')
        for entry in to_add:
            f.write(entry + '\n')
    for entry in to_add:
        print(f'Added to .gitignore: {entry}')


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Register already-copied portable legacy files in the tracked SHA-256 manifest.'
    )
    ap.add_argument(
        'paths', nargs='+',
        help='Files/directories relative to reference_outputs/legacy_full (or absolute paths inside it).'
    )
    ap.add_argument('--store', type=Path, default=DEFAULT_STORE)
    ap.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument('--gitignore', type=Path, default=DEFAULT_GITIGNORE)
    ap.add_argument('--no-gitignore', action='store_true')
    args = ap.parse_args()

    store = args.store.expanduser().resolve()
    manifest = args.manifest
    if not store.is_dir():
        raise SystemExit(f'Legacy store not found: {store}')
    if not manifest.is_file():
        raise SystemExit(f'Manifest not found: {manifest}\nRun the main legacy import first.')

    requested = []
    for raw in args.paths:
        path, rel = normalize_requested(store, raw)
        if not path.exists():
            raise SystemExit(f'Requested legacy path does not exist:\n  {path}')
        requested.append((path, rel))

    with manifest.open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or not set(FIELDS).issubset(reader.fieldnames):
            raise SystemExit(
                f'Unexpected manifest format in {manifest}; expected at least columns {FIELDS}'
            )
        rows = list(reader)

    kept_rows = []
    removed = 0
    for row in rows:
        row_rel = row.get('relative_path', '')
        replace = any(requested_matches(row_rel, rel, path) for path, rel in requested)
        if replace:
            removed += 1
        else:
            kept_rows.append(row)

    new_rows = []
    roots = set()
    for requested_path, requested_rel in requested:
        files = list(iter_files(requested_path))
        if not files:
            raise SystemExit(f'No files found under requested path:\n  {requested_path}')
        print(f'\nRegistering: {requested_rel}')
        for file_path in files:
            rel = file_path.relative_to(store).as_posix()
            root = Path(rel).parts[0]
            roots.add(root)
            size = file_path.stat().st_size
            digest = sha256_file(file_path)
            new_rows.append({
                'root': root,
                'relative_path': rel,
                'size_bytes': str(size),
                'sha256': digest,
            })
            print(f'  {rel}')
            print(f'    size:   {size} bytes')
            print(f'    sha256: {digest}')

    combined = kept_rows + new_rows
    deduped = {}
    for row in combined:
        deduped[row['relative_path']] = {field: row.get(field, '') for field in FIELDS}
    final_rows = sorted(deduped.values(), key=lambda r: r['relative_path'])

    tmp = manifest.with_suffix(manifest.suffix + '.partial')
    with tmp.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(final_rows)
    tmp.replace(manifest)

    if not args.no_gitignore:
        ensure_gitignore(roots, args.gitignore)

    print()
    print(f'Removed/replaced prior manifest rows: {removed}')
    print(f'Registered files this run: {len(new_rows)}')
    print(f'Manifest rows now: {len(final_rows)}')
    print(f'Updated manifest: {manifest}')
    print()
    print('Next:')
    print('  make validate-legacy')
    print('  make activate-legacy')


if __name__ == '__main__':
    main()
