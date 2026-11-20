# Selection-stability metrics: definitions and checks

This file defines the quantities in the multi-seed selection-stability outputs: the K=50 tables,
the top-5 vs top-10 neighbourhood summary and the M0-M4 supplementary table. Everything here is
spatial participant-held-out selection on the every-K, 10-seed run. Fit seeds are 0-9, and seed 0
is the reference (canonical) fit.

Scripts (in this folder):

| Script | Output folder (under `outputs/multiseed_stability/`) |
|---|---|
| `k50_stability/nesting_check.py` | `nesting_check/` |
| `k50_stability/k50_cross_seed_stability.py` | `k50_cross_seed_stability/` |
| `plots/topm_stability_neighbourhood.py` | `topm_stability_neighbourhood/` |
| `k50_stability/make_supp_tables.py` | `k50_supp_tables/` |

## 1. Raw nesting check (before consensus and before alignment)

`nesting_check.py` reads the raw per-ranking rows (`heldout_topm_runs_all_seeds.csv.gz`). These
use each seed's own component IDs, with no Hungarian matching and no aggregation. For every
orientation x K x condition x fit seed x outer split x direction, it asserts that
`set(top5_raw_ids) <= set(top10_raw_ids)`, and that the top-5 is the ordered prefix of the top-10.
m is clamped to K, as in the pipeline. It also checks that every K x condition x seed has
40 rankings, and that the per-seed selection-stability tables match a recount from the raw rows.

**Consensus sets are not always nested.** The consensus top-5 and consensus top-10 are built
separately (Sec. 2d). So a seed's consensus top-5 can include an archetype that is missing from its
consensus top-10, even though every raw top-5 is nested in its raw top-10. The cases are listed
in `nesting_check/consensus_nesting_by_seed.csv`.

## 2. Metric definitions

**a. Per-ranking selection (raw).** Each outer split r (r = 0...19) divides the 36 participants
into halves A and B (18/18). The same split is used for every fit seed at a given K and
condition. The direction A->B ranks archetypes on A: each archetype is decoded on its own, then
archetypes are sorted by mean accuracy, with ties broken by mean decoder rank and then by index.
The top-m of that ranking is its first m entries, and they are scored on B. B->A does the reverse.
This gives 20 splits x 2 directions = 40 rankings per K x condition x seed.

**b. Hungarian matching cost.** For each K and each seed s, the cost matrix is -r, where r is the
Pearson correlation between the seed-0 and seed-s spatial maps. Each map is the subject-mean map
(Smean) of one archetype, a 700-node vector. The match is a one-to-one assignment over all K
archetypes (`scipy.optimize.linear_sum_assignment`). For the M0-M4 motifs at K != 50, all 50 K=50
seed-0 maps are matched against that fit's K maps, with the matrix padded to square.

**c. Recovery.** An archetype counts as recovered in seed s when its Hungarian partner's map
correlation is at least 0.90. Seed 0 is the reference, so it always counts as recovered. Matched
r is bimodal (near 1 or below about 0.4), so 0.90 is a clean cut.

**d. Consensus set per seed.** Selection frequency is the number of the seed's 40 rankings whose
top-m contains the archetype, divided by 40. The consensus top-m is the m archetypes with the
highest frequency, with ties broken by lower mean rank position. Top-5 and top-10 use the top-5
and top-10 rows respectively, so they are formed independently (see the caveat in Sec. 1).

**e. Recovery-aware labels.** Each pick in a seed's set is relabelled. It gets the seed-0 ID of
its Hungarian partner if that match is recovered (r >= 0.90). Otherwise it gets a unique
seed-specific label (`u{seed}_{id}`), which cannot overlap with anything. The *forced* variant
always uses the Hungarian partner, whatever r is, and is an upper bound. The recovery-aware
Jaccard is the main stability metric.

**f. Jaccard and chance correction.** J(A,B) = |intersect(A, B)| / |union(A, B)|. The seed-set Jaccard is the mean
over all 45 pairs of seeds' consensus sets. The chance value J_rand is the exact expected Jaccard
of two independent uniform random m-subsets of K archetypes:
J_rand = sum_i P(|intersect(A, B)| = i) * i / (2m - i), where P is hypergeometric. The chance-corrected value is
(mean J - J_rand) / (1 - J_rand), so the correction is applied to the pair-averaged J.

**g. Directions (A->B vs B->A).** The two directions are not merged into one ranking before the
seed summaries. Each one is its own ranking and adds 1/40 to the selection frequency. This is the
same as averaging the two direction indicators within a split and then averaging over the 20
splits. In "same split, across seeds", each direction is compared only with the same direction
in the other seed. In "across splits, within seed", all 780 pairs of the 40 rankings are
compared, including the A->B / B->A pair from the same split.

**h. Per-motif quantities (M0-M4).**
- *Recovered*: the number of seeds 1-9 in which the motif is recovered (Sec. 2c).
- *Selected & recovered*: the number of seeds (out of 10) in which the motif is recovered and
  appears in the seed's consensus top-m.
- *Selected but not recovered* (`selected_not_recovered`): the number of seeds in which the
  motif's forced partner (r < 0.90) appears in the consensus top-m. These are least-bad pairings,
  not genuine recoveries. In the K=50 summary across all archetypes, it counts (seed 1-9, pick)
  occurrences whose Hungarian partner has r < 0.90, so its maximum is 9 x m.
- *Selection given recovery*: selected & recovered divided by recovered seeds, with seed 0
  included (out of 10).
- *Split-level frequency*: the motif's 40-ranking selection frequency, averaged over recovered
  seeds.

**i. Aggregate Yeo-7 profile.** For each seed's consensus set, the profile is computed from the
selected maps themselves, not from their assigned seed-0 labels. Each map's |weight| is summed
per Yeo-7 network (`configs/node_yeo7.csv`) and normalised to shares, and the shares are averaged
over the set. The seed-to-seed agreement is the Pearson r between these 7-vectors.

## 3. Averaging hierarchy

- **Selection frequency:** directions within outer split, then outer splits within fit seed.
  This holds exactly, because every K x condition x seed has 2 directions x 20 splits (asserted).
- **Per-seed metrics** (consensus sets, M0-M4 counts, within-seed Jaccard): computed within each
  seed, then summarised across the 10 seeds.
- **Across-seed metrics** (seed-set Jaccard): the last level is the 45 seed pairs, because these
  compare seeds. "Same split, across seeds" is computed as seed pairs within split x direction,
  then directions within split, then splits. "Across splits, within seed" is computed per seed,
  then averaged across seeds. The design is balanced, so these hierarchical versions (the `*_hier`
  columns) equal the pooled means to rounding error.

## 4. M0-M4 crosswalk

M0-M4 are the only named motifs. They come from the representative K=50 spatial across-condition
fit. They were chosen by participant-held-out top-5 selection frequency pooled across the 10 fit seeds, as stated
in the paper's Methods:

- For each archetype and condition, selection frequency is the fraction of the 40 cross-fitted
  rankings (20 outer splits x 2 directions) in which it appears in the top 5.
- Frequencies are averaged over fit seeds 0-9 after aligning each seed to seed 0 by Hungarian
  matching of spatial maps (Sec. 2). A seed that does not recover a motif contributes the
  frequency of its forced match.
- Each archetype is scored by its highest 10-seed mean frequency over the three conditions, and the
  five highest are kept after excluding archetypes not reproducibly recovered across seeds. At
  K=50 that excludes canonical 31 (recovered, matched map r >= 0.90, in only 5 of 9 seeds).
  K=50 ranking: 28 .995 (word), 16 .885 (intact), 39 .698 (word), 49 .622 (intact),
  37 .508 (intact); next 36 .460, 31 .393 (excluded). Source:
  `heldout-everyk-spatial-all-results/selection_by_seed.csv`.
- The motifs are displayed from, and numbered by, the seed-0 fit.

The ids are defined once in `src/msaa_pipeline/paper_archetypes.py`.

| Motif | K=50 canonical id (seed 0, 0-based) |
|---|---|
| M0 | 28 |
| M1 | 39 |
| M2 | 16 |
| M3 | 37 |
| M4 | 49 |

`make_supp_tables.py` writes the crosswalk with each motif's dominant Yeo-7 networks
(`k50_supp_tables/M0M4_crosswalk.csv`), along with the K=50 per-motif table defined in Sec. 2h.
