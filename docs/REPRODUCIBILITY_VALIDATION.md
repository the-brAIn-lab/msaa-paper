# Reproducibility validation and portable reuse of the original full sweep

## Policy

The cleaned repository supports two complementary modes.

1. **Regeneration mode** starts from the two HTFA MAT inputs and recomputes
   selected MS-AA fits/decoding using `smoke`, `mini-paper`, or the optional full
   profile.
2. **Validated legacy mode** reuses the original expensive full-sweep derived
   outputs after representative recomputation established parity. The bytes are
   copied into `reference_outputs/legacy_full/` inside the project folder, while
   a tracked manifest records relative paths, byte sizes, and SHA-256 hashes.

This keeps the working project self-contained for a student handoff without
putting multi-gigabyte outputs into Git.

## Migration validation completed

The `mini-paper` profile was rerun from the two MAT inputs using all participants
and the original fit/decoding settings. The migration checks established:

- core full and single-archetype decoding matched the original results within
  the predefined numerical tolerances;
- the spatial across-condition K=46 solution reproduced essentially exactly:
  all 46 archetypes matched one-to-one, mean matched spatial-map correlation was
  approximately 0.99999965, and the minimum matched correlation was
  approximately 0.99999332;
- the five manuscript-selected K=46 archetypes retained the same component IDs
  (4, 5, 6, 12, 8);
- the original participant-held-out outer split seeds are frozen in
  `configs/heldout_outer_seeds.csv` so the cleaned held-out analysis uses the
  same participant halves as the manuscript analysis;
- aggregated held-out top-m results across the K=5/K=46 validation differed from
  the original by at most approximately 0.000458 in mean accuracy and 0.000196
  in SEM.

Small split-level differences remain because tiny numerical changes in
single-archetype scores can reorder near-tied components at a top-m boundary.
They do not reflect a different K=46 decomposition or a materially different
aggregate result.

## Import the original derived outputs into this project

From the cleaned repository root:

```bash
make import-legacy OLD_REPO=/Users/lowen/repos/arch_analysis/MS_AA
```

The import is resumable. It copies configured paper-supporting output roots into:

```text
reference_outputs/legacy_full/
```

and writes:

```text
configs/legacy_outputs_manifest.csv
```

The large store is ignored by Git. The manifest is intentionally tracked and
contains only portable repo-relative paths, byte sizes, and SHA-256 hashes.

The roots to import are documented in `configs/legacy_output_roots.csv`. Core
fit/decoding outputs are required; smaller downstream/supporting roots are marked
optional so historical repos that lack one can still be imported.

## Validate the portable copy

A fast existence/size check:

```bash
make validate-legacy
```

A full cryptographic verification (slower because it rereads all legacy bytes):

```bash
make validate-legacy-hash
```

## Activate historical paths

Older notebooks expect output folders such as `msaa_flexible_outputs_npz/` at
the repository root. Expose the internal portable store with:

```bash
make activate-legacy
```

This creates relative symlinks only. They point into
`reference_outputs/legacy_full/`, not to the original repository on Lucy's
machine. Thus no absolute machine-specific path is embedded in the analysis.

Before intentionally running a full from-scratch fit/decode, deactivate them:

```bash
make deactivate-legacy
```

The Makefile also guards `full-fit`, `full-core`, and `full-heldout` and refuses
to run while these links are active.

## Student handoff

The complete working folder contains three classes of material:

- **tracked Git content:** source code, notebooks, configurations, documentation,
  frozen split seeds, and the legacy-output manifest;
- **ignored scientific inputs:** `data/input/pieman_data.mat` and
  `data/input/pieman_posterior_K700.mat`;
- **ignored validated derived outputs:** `reference_outputs/legacy_full/`.

A Git clone contains only the first class. For a student who should receive the
ready-to-work project, copy/archive the **whole project folder including ignored
files**. On the receiving machine:

```bash
make handoff-check
make activate-legacy
```

If the transfer mechanism does not preserve symlinks, that is fine: the legacy
bytes remain in the project and `make activate-legacy` recreates the relative
links.

## Provenance rule

Files under `reference_outputs/legacy_full/` must be treated as archived derived
artifacts. Do not edit them in place. Fresh recomputation belongs in the normal
current pipeline outputs. If the archived legacy bytes ever change, re-import
from the authoritative original source and commit the newly generated manifest
only after reviewing why the hashes changed.
