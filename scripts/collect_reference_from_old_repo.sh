#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 /path/to/original/MS_AA" >&2
  exit 2
fi
OLD="$(cd "$1" && pwd)"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$HERE/reference/original_outputs"
mkdir -p "$DEST"
for d in \
  msaa_condrank_decoding_outputs_spatial_across \
  msaa_condrank_decoding_outputs_temporal_across \
  msaa_condrank_decoding_outputs_spatial_within \
  msaa_condrank_decoding_outputs_temporal_within; do
  if [ -d "$OLD/$d" ]; then
    mkdir -p "$DEST/$d"
    rsync -a --include='*.csv' --include='rankings.npy' --exclude='*' "$OLD/$d/" "$DEST/$d/"
  fi
done
# Copy only compact fitted models used by mini-paper comparison.
mkdir -p "$DEST/msaa_flexible_outputs_npz"
for f in \
  spatialAA_across_acrossCond_K5.npz \
  spatialAA_across_acrossCond_K46.npz \
  temporalAA_across_acrossCond_K5.npz \
  spatialAA_within_intact_K5.npz spatialAA_within_word_K5.npz spatialAA_within_rest_K5.npz \
  temporalAA_within_intact_K5.npz temporalAA_within_word_K5.npz temporalAA_within_rest_K5.npz; do
  [ -f "$OLD/msaa_flexible_outputs_npz/$f" ] && cp "$OLD/msaa_flexible_outputs_npz/$f" "$DEST/msaa_flexible_outputs_npz/"
done
echo "Reference subset copied to $DEST"
