# 08 Multi-seed stability

Reruns every MS-AA fit with nine more FurthestSum starting picks (seeds 1-9, next to the
paper's seed 0), then asks what stays the same: archetypes matched one-to-one to seed 0,
decoding and held-out top-m accuracy averaged over the 10 fits, and the stability of the
selected archetypes, their maps, networks and time courses.

Everything runs on one machine from the repo root. All paths are relative to the repo
(see `paths.py`), so no script needs a path argument.

## 1. Setup

From the repo root, as in the main README:

```bash
uv sync && source .venv/bin/activate        # or: pip install -r requirements.txt (Python 3.12)
# put pieman_data.mat and pieman_posterior_K700.mat in data/input/
make check
```

The node-to-network table `configs/node_yeo7.csv` is committed. To rebuild it from the
atlas (this needs network access once):

```bash
python scripts/fetch_reference_atlas.py
python scripts/build_node_yeo7.py
```

## 2. Run the fits and held-out top-m

`run_local.py` runs the whole chain in order. Every step is resumable: rerunning it
skips finished fits, decodes and held-out splits.

| step | what it does |
|---|---|
| `seed0` | the paper's reference fits, using the frozen picks in `configs/furthestsum_init.csv.gz`; then decode, held-out top-m at the 30 paper K, and merge |
| `picks` | FurthestSum starting picks for seeds 1-9 |
| `fits` | one fit per seed x analysis x scope x K (x condition) |
| `decode` | decode each seed's fits (exact low-rank decoder) |
| `topm` | same-data top-m for each seed (diagnostic only; the paper uses held-out top-m) |
| `merge` | cross-seed tables: Hungarian matching to seed 0, per-K mean and SD, ... |
| `heldout30` | participant-held-out top-m at the 30 paper K, every seed |
| `heldoutK` | participant-held-out top-m at every K, every seed |

```bash
# small test, about 30 min on 4 cores: one K, all 10 seeds
python analysis/08_multiseed_stability/run_local.py --ks 5 --seeds 1,2,3,4,5,6,7,8,9 --workers 4

# what the figures below need: every K, all seeds (very large, see below)
python analysis/08_multiseed_stability/run_local.py --workers 32 --threads 1
```

Options:
- `--steps`: a comma list of the steps above.
- `--ks`: limits K (the default is every K: spatial 2-700, temporal 2-300).
- `--seeds`: the extra seeds (seed 0 always runs).
- `--jobs`: for example `spatial_across,temporal_across`.
- `--workers` and `--threads`: the number of processes, and the BLAS threads per process.

**Cost.** The full run is about 36,000 fits plus held-out splits at every K. Held-out at
every K alone is roughly 430 CPU-hours per seed for spatial and 90 for temporal. The paper
run took several days on about 100 cluster cores. Locally it is weeks of CPU time.

**Outputs** (in the repo root, gitignored):

```text
msaa_flexible_outputs_npz/                       seed-0 fits
msaa_condrank_decoding_outputs_<analysis>_across/ seed-0 decoding
seeds/full/seedNN/                               one self-contained run root per seed (its slim/ holds compact fits)
seeds/slim/seed00/                               compact copies of the seed-0 fits
seeds/results/                                   per_K.csv, matching_*.csv, jumpiness.csv, ...
heldout-seeds-results/                           held-out top-m, 30 paper K
heldout-everyk-<spatial|temporal>-<seed0|all>-results/   held-out top-m, every K
```

## 3. Make the figures and tables

Each script writes to `outputs/multiseed_stability/<script name>/` (gitignored). The
archetype IDs below are 0-based seed-0 IDs, as used in the paper. Every input and output
can be overridden with flags (`--help` lists them).

| figure / table | command |
|---|---|
| Peak-K bootstrap CIs, used by Fig 3 | `for a in spatial temporal; do for w in seed0 all; do python plots/peakk_ci.py $a --which $w --W 1; done; done` |
| Fig 3 by analysis: seed 0 vs 10-fit mean, decoding and held-out top-m | `python plots/fig3_by_analysis.py spatial` and `... temporal` |
| K=51 Yeo-7 composition bars, all seeds | `python plots/yeo7_bars_all_seeds.py 51 28,40,16,50,37` |
| K=51 glass brains, all seeds | `python plots/glassbrains_all_seeds.py 51 28,40,16,50,37` |
| K=51 per-seed network figures (`K51_seed{0-9}_networks`) and cross-seed similarity | `python plots/loadx.py` (once), then `python plots/fig5_per_seed.py 51 28,40,16,50,37` |
| K=14 temporal motifs, all seeds | `python plots/temporal_motifs_by_seed.py 14 9,4,5,7,0` |
| Top-5 vs top-10 selection stability across the K neighbourhood; M0-M4 top-5 vs top-10 | `python plots/topm_stability_neighbourhood.py 45,46,47,48,49,50,51,52,53,54,55` |
| M0-M4 crosswalk and K=46 supplementary table (run after the line above) | `python k46_stability/make_supp_tables.py` |
| K=46 paper archetypes across seeds; top-5 overlap (K46, K51, temporal K14) | `python k46_stability/k46_cross_seed_stability.py` |
| Raw top-5 vs top-10 nesting check | `python k46_stability/nesting_check.py` |
| Leave-one-seed-out check of the good-K bands on held-out top-m accuracy | `python plots/goodk_loso.py spatial` (and `temporal`; `--W 9` default, `--W 21` as in peak-K CI) |
| Time-tolerant matching of temporal archetypes between two runs | `python temporal_drift/time_tolerant_match.py --a RUN_A --b RUN_B` |

Run them from this folder, or prefix each path with `analysis/08_multiseed_stability/`.

`DEFINITIONS.md` defines every selection-stability metric in these tables: recovery, the Hungarian
cost, recovery-aware labels, consensus sets, chance-corrected Jaccard and the averaging order.
- The paper's K=46 archetypes are M0-M4 = #28, 35, 16, 39, 45 (seed-0 ids;
  `src/msaa_pipeline/paper_archetypes.py`). `k46_cross_seed_stability.py` uses them by
  default; set `PAPER_SET=28,35,37,31,39` (or any list) to check another set, e.g. the
  submitted one.
- `fig5_per_seed.py` at K=46 also writes the matching table that
  `k46_cross_seed_stability.py` checks itself against when it is present.

To make figures from a finished run somewhere else (for example results copied from a
cluster), point the scripts at that run root:

```bash
MSAA_ROOT=/path/to/run_root python plots/fig3_by_analysis.py spatial
```

## 4. What to expect when you rerun

- **Spatial results reproduce across machines:** archetype numbering is identical, and
  matched maps have r >= 0.9999.
- **Summaries reproduce across machines:** decoding and held-out summaries agree within
  about 0.005.
- **Temporal archetypes at high K may differ:** a few can come out slightly different on
  other hardware or BLAS builds (local optima). Most such differences are shifts of 1-3 TRs.
  `temporal_drift/time_tolerant_match.py` measures this.
- **Nothing forces bit-identical results:** the BLAS is whatever the numpy wheel ships.

## Supporting scripts (in `scripts/`)

- `run_parallel.py`: the parallel fit / decode / held-out / merge driver that `run_local.py` calls.
- `freeze_init.py`, `freeze_everyk.py`: compute the frozen FurthestSum picks in
  `configs/furthestsum_init.csv.gz`. Existing rows are never changed, so extending the
  K grid does not renumber old fits.
- `verify_lowrank.py`: checks the low-rank decoder against the dense one
  (`MSAA_REF_RUN=/path/to/reference_run python scripts/verify_lowrank.py --quick`).
- `build_node_yeo7.py`: builds `configs/node_yeo7.csv`.
