# Portable validated legacy outputs

`reference_outputs/legacy_full/` is the local, portable store for the expensive
full-analysis outputs produced by the original repository.

The directory is intentionally ignored by Git because it is several gigabytes.
A tracked SHA-256 manifest in `configs/legacy_outputs_manifest.csv` records the
exact files that were imported.

Import once from the original repository:

```bash
make import-legacy OLD_REPO=/path/to/original/MS_AA
```

Then expose the stored directories at the historical repo-relative paths used by
older notebooks:

```bash
make activate-legacy
```

Activation creates **relative symlinks within this repository**, so the project
folder can be moved to another machine. If an archive/copy operation does not
preserve symlinks, simply run `make activate-legacy` again after unpacking.

A Git clone alone will not contain `reference_outputs/legacy_full/`. To hand the
complete working project to another person, copy/archive the whole project
folder, including ignored files, and have them run:

```bash
make handoff-check
make activate-legacy
```
