# Canonical sparse analysis: participant-held-out top-m only

The cleaned repository intentionally **does not regenerate the historical same-data top-m analysis**.

The only canonical top-m procedure is participant-held-out cross-fitting:

1. For a condition and model order K, split participants into outer halves A and B.
2. Rank individual archetypes using A only.
3. Select top-m from the A-derived ranking and evaluate only on B.
4. Reverse the direction (B ranks, A evaluates).
5. Repeat across deterministic outer participant splits.
6. Average A->B and B->A within each outer split before computing SEM across outer splits.

Canonical outputs are written to:

```text
msaa_heldout_topm_outputs_spatial_across/
  heldout_topm_crossfit_runs.csv
  heldout_per_archetype_ranking_runs.csv
  heldout_topm_summary.csv
  heldout_archetype_selection_stability.csv
```

The historical same-data top-m notebooks remain only under `legacy/original_notebooks/` for provenance.
They are not part of the reproducible paper pipeline.

The core decoding directories named `msaa_condrank_decoding_outputs_*` are retained for compatibility
with downstream full-reconstruction and single-archetype diagnostics, but the cleaned runner no longer
creates `topm_summary.csv` or rankings for same-data top-m selection there.
