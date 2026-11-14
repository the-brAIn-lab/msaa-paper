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
| `merge` | cross-seed tables: Hungarian matching to seed 0, per-K mean and SD, ... |
| `heldout30` | participant-held-out top-m at the 30 paper K, every seed |
| `heldoutK` | participant-held-out top-m at every K, every seed |
| `highm` | extra m (7-100) at the 15 matched spatial/temporal K pairs, reusing the every-K rankings (Fig 4C, S01) |
| `landscape` | extra m (7-50) at spatial K = 10, 20, ..., 700, same reuse (S12) |

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
seeds/results/                                   per_K.csv, matching_*.csv, paper_K46.csv, ...
heldout-seeds-results/                           held-out top-m, 30 paper K
heldout-everyk-<spatial|temporal>-<seed0|all>-results/   held-out top-m, every K
heldout-highm-results/, heldout-highm-landscape-results/   held-out top-m at high m (checked against the every-K run)
```

## 3. Make the figures and tables

Each script writes to `outputs/multiseed_stability/<script name>/` (gitignored). The
archetype IDs below are 0-based seed-0 IDs, as used in the paper. Every input and output
can be overridden with flags (`--help` lists them).

| figure / table | command |
|---|---|
| Fig 4A/B (10-seed held-out top-m, per condition and in the one-axis layout), Fig 4D (temporal matched retained fraction, 25/50/70 %), full-model held-out accuracy vs K (run after the matched-fraction section below) | `python plots/fig4_heldout_panels.py` |
| Fig 4C and S01: spatial - temporal held-out top-m at the 15 exactly matched K/D ratios (m = 1-100), mean of seeds and number of seeds with spatial > temporal; needs `heldout-highm-results/` (step `highm`) | `python plots/fig4c_matched_ratio_heatmaps.py` |
| Supplement: complete held-out spatial top-m landscape over K x m (accuracy and retained fraction; m = 1-50, K = 10-700 in steps of 10); needs `heldout-highm-landscape-results/` (step `landscape`) | `python plots/spatial_topm_landscape.py` |
| K=46 per-seed matching to seed 0, map and time-course r for every archetype (needed by the K=46 stability script); also draws one M0-M4 page per seed | `python plots/loadx.py` (once), then `python plots/fig5_per_seed.py 46 28,35,16,39,45` |
| Top-5 vs top-10 selection stability across the K neighbourhood; M0-M4 top-5 vs top-10 | `python plots/topm_stability_neighbourhood.py 45,46,47,48,49,50,51,52,53,54,55` |
| M0-M4 crosswalk and K=46 supplementary table (run after the line above) | `python k46_stability/make_supp_tables.py` |
| K=46 paper archetypes (M0-M4) across seeds: map r, time-course r, Yeo-7 r, recovered; top-5 overlap (run after `fig5_per_seed.py 46`) | `python k46_stability/k46_cross_seed_stability.py` |
| Raw top-5 vs top-10 nesting check (checks the S05/S06 selection tables; see `DEFINITIONS.md`) | `python k46_stability/nesting_check.py` |
| Compact M0-M4 summary across seeds (map r, temporal profile r, Yeo-7 r, recovered at r >= 0.9; run after the K=46 stability script) | `python k46_stability/m0m4_summary.py` |
| Temporal K=14 supplement: accumulation, generator-anchor masking, motifs at higher K (run the revision_10a and revision_10c notebooks first) | `python plots/k14_supplement.py` |

Run them from this folder, or prefix each path with `analysis/08_multiseed_stability/`.

`DEFINITIONS.md` defines every selection-stability metric in these tables: recovery, the Hungarian
cost, recovery-aware labels, consensus sets, chance-corrected Jaccard and the averaging order.
- The paper's K=46 archetypes are M0-M4 = #28, 35, 16, 39, 45 (seed-0 ids;
  `src/msaa_pipeline/paper_archetypes.py`). `k46_cross_seed_stability.py` uses them by
  default; set `PAPER_SET=28,35,37,31,39` (or any list) to check another set, e.g. the
  submitted one.
- `fig5_per_seed.py` at K=46 writes the per-seed similarity table (including the
  time-course r) that `k46_cross_seed_stability.py` reads and checks its matching against.

To make figures from a finished run somewhere else (for example results copied from a
cluster), point the scripts at that run root:

```bash
MSAA_ROOT=/path/to/run_root python plots/fig4_heldout_panels.py
```

## 3b. Temporal matched retained fraction and the K=14 contrast

Fixed top-m favours small K because the same m is a larger share of the basis. This asks
the question at a fixed fraction instead: with m = round(fraction x K) temporal archetypes,
selected on one participant half and scored on the other (the same splits, rankings and
decoder seeds as the held-out top-m above), what fraction of the full model's above-chance
accuracy is retained? It needs the `heldoutK` step for temporal (its per-split rankings in
`heldout-everyk/parts/temporal/`) and every seed's temporal across fits; nothing is refit.
Run from the repo root (`R=.`), or set `R` to another run root:

```bash
cd analysis/08_multiseed_stability/matched_fraction
R=../../..; B=$R/outputs/multiseed_stability/matched_fraction; mkdir -p $B
KS=5,8,10,12,14,16,20,25,30,40,50; CONDS=intact,word; SEEDS=0:1:2:3:4:5:6:7:8:9
python collect_rankings.py $R $B/rankings.csv ${KS//,/:} ${CONDS//,/:} $SEEDS 20
roots=(--fit-root $R/msaa_flexible_outputs_npz)
for s in 01 02 03 04 05 06 07 08 09; do roots+=(--fit-root $R/seeds/full/seed$s/msaa_flexible_outputs_npz); done
python temporal_matched_fraction.py build-manifests "${roots[@]}" --k-values $KS --conditions $CONDS \
  --fit-manifest $B/fit_manifest.csv --task-manifest $B/task_manifest.csv      # 220 tasks
python temporal_matched_fraction.py audit --task-manifest $B/task_manifest.csv \
  --ranking-csv $B/rankings.csv --outer-reps 20 --output $B/rankings_audit.csv
N=$(($(wc -l < $B/task_manifest.csv) - 1))
for i in $(seq 0 $((N - 1))); do        # independent and resumable: run them in parallel if you like
  python temporal_matched_fraction.py run-task --task-manifest $B/task_manifest.csv --task-id $i \
    --ranking-csv $B/rankings.csv --output-dir $B --fractions 0.25,0.50,0.70,1.00 \
    --outer-reps 20 --nfolds 2 --nreps 20 --base-seed 42
done
python validate.py $R $B            # must say PASS: identical to held-out top-m where m is 1/2/3/5/10
python temporal_matched_fraction.py aggregate --task-manifest $B/task_manifest.csv --output-dir $B \
  --fractions 0.25,0.50,0.70,1.00 --outer-reps 20 --require-complete
python k14_local_contrast.py $B/aggregate $B/k14_local_contrast
```

`aggregate/` holds the retention curves and tables, averaged direction -> outer split -> fit
seed (`matched_fraction_summary.csv`, `matched_fraction_seed_summary.csv`, one plot per
condition and per fraction). `k14_local_contrast.py` writes the seed-level table (seed,
condition, K, requested fraction, selected m, actual fraction m/K, raw subset and
full-model accuracy, above-chance retention), the paired contrast
Delta14 = R(14) - [R(12) + R(16)] / 2 with a bootstrap 95% CI and the number of seeds with
Delta14 > 0, and plots of raw subset accuracy, full-model accuracy and per-seed K=12/14/16.
At 50% the fraction is exactly matched there (6/12 = 7/14 = 8/16). Each task takes a minute
or two on 2 cores; the 220 tasks are a few CPU-hours.

## 4. What to expect when you rerun

- **Spatial results reproduce across machines:** archetype numbering is identical, and
  matched maps have r >= 0.9999.
- **Summaries reproduce across machines:** decoding and held-out summaries agree within
  about 0.005.
- **Temporal archetypes at high K may differ:** a few can come out slightly different on
  other hardware or BLAS builds (local optima). Most such differences are shifts of 1-3 TRs.
- **Nothing forces bit-identical results:** the BLAS is whatever the numpy wheel ships.

## Supporting scripts (in `scripts/`)

- `run_parallel.py`: the parallel fit / decode / held-out / merge driver that `run_local.py` calls.
- `freeze_init.py`, `freeze_everyk.py`: compute the frozen FurthestSum picks in
  `configs/furthestsum_init.csv.gz`. Existing rows are never changed, so extending the
  K grid does not renumber old fits.
- `build_node_yeo7.py`: builds `configs/node_yeo7.csv`.
