# Build validation

The cleaned core runner was tested during repository construction with temporary
synthetic MATLAB inputs that matched the required interface:

- cell-array conditions `intact`, `word`, and `rest`;
- a rest series 100 time points longer so the historical rest trim yielded the
  same analysis length;
- a MATLAB `posterior` struct with `centers` and `widths`.

The test successfully ran:

1. input validation;
2. spatial across-condition K=5 MS-AA fitting;
3. temporal across-condition K=5 MS-AA fitting;
4. full-reconstruction decoding;
5. per-archetype decoding;
6. condition-specific ranking;
7. top-m decoding;
8. QC-table and QC-figure generation;
9. a second run to verify resumable behavior.

The synthetic MAT files and generated outputs were removed afterward and are not
included in the distributed repository.

This test validates code/path mechanics only. Scientific numerical equivalence
should be checked with the real MAT inputs via the `mini-paper` profile and the
old-repository comparison command.
