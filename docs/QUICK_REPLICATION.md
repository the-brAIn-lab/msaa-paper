# Quick replication modes

The cleaned repository contains two deliberately small verification profiles.

## 1. `smoke`: structural end-to-end check

```bash
make smoke
```

This profile uses:

- spatial across-condition MS-AA: `K=5`
- temporal across-condition MS-AA: `K=5`
- first 4 participants per condition
- 20 fitting iterations
- 2 decoding repetitions
- top-m = 1, 2, 3

It runs:

`MAT inputs -> normalization -> MS-AA fit -> reconstruction -> per-archetype decoding -> condition-specific ranking -> top-m decoding -> QC figure`

This is intentionally fast and **is not expected to reproduce manuscript numbers** because it uses fewer participants and fewer optimization/decoding iterations.

Outputs are written under `outputs/smoke/` and cannot overwrite full-paper results.

## 2. `mini-paper`: small exact-replication check

```bash
make mini
```

This profile keeps the original computational settings and all participants, but
fits only a small subset of model orders:

- spatial across: K = 5, 46
- temporal across: K = 5
- spatial within: K = 5 for each condition
- temporal within: K = 5 for each condition
- 200 fit iterations
- 20 decoding repetitions
- seed = 42

K=5 gives a common checkpoint across spatial/temporal and within/across models.
K=46 additionally exercises the representative spatial across-condition fit used
later for manuscript interpretation.

To compare the new small run directly with the corresponding rows/files in the
old repository:

```bash
make compare OLD_REPO=/path/to/original/MS_AA
```

The comparison writes `outputs/mini_paper/comparison_to_original.json`. For fit
files, it compares reconstructed data rather than only raw archetype ordering,
which makes the check robust to a potential component permutation.

## Why both modes?

`smoke` answers: **Does the reorganized code still execute end-to-end?**

`mini-paper` answers: **When we use the original settings on a few selected K
values, do we recover the old analysis?**
