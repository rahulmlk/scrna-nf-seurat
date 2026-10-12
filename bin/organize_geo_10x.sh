#!/usr/bin/env bash
# Sort flat GEO supplementary downloads into one folder per sample, with the
# file names Seurat's Read10X expects (barcodes / features / matrix).
#
# GEO usually stores files as <GSM>_<label>_barcodes.tsv.gz, <GSM>_<label>_features.tsv.gz,
# <GSM>_<label>_matrix.mtx.gz. This groups them by the GSM prefix.
#
# Usage: organize_geo_10x.sh <download_dir> <output_dir>
# Check the printed summary: every sample should end up with 3 files.

set -euo pipefail

src="${1:?usage: $0 <download_dir> <output_dir>}"
dst="${2:?usage: $0 <download_dir> <output_dir>}"
mkdir -p "$dst"

shopt -s nullglob
for f in "$src"/GSM*; do
    base=$(basename "$f")
    gsm="${base%%_*}"

    case "$base" in
        *barcodes*) name="barcodes.tsv" ;;
        *features*|*genes*) name="features.tsv" ;;
        *matrix*)   name="matrix.mtx" ;;
        *) echo "skipping unrecognised file: $base" >&2; continue ;;
    esac

    # keep .gz if the source is gzipped (Read10X handles both)
    [[ "$base" == *.gz ]] && name="$name.gz"

    mkdir -p "$dst/$gsm"
    cp "$f" "$dst/$gsm/$name"
done

rd() { if [[ "$1" == *.gz ]]; then zcat "$1"; else cat "$1"; fi; }

echo "Per-sample summary (expect 3 files each; check species in the feature names):"
for d in "$dst"/*/; do
    id=$(basename "$d")
    n=$(find "$d" -type f | wc -l)
    bc=$(ls "$d"/barcodes.tsv* 2>/dev/null | head -1 || true)
    ft=$(ls "$d"/features.tsv* 2>/dev/null | head -1 || true)
    [[ "$n" -ne 3 ]] && echo "WARNING: $id has $n files, expected 3" >&2
    cells="NA"; genes="NA"
    [[ -n "$bc" ]] && cells=$(rd "$bc" | wc -l)
    [[ -n "$ft" ]] && genes=$(rd "$ft" | head -3 | cut -f2 | paste -sd, -)
    printf '%s\t%s files\t%s cells\tfirst genes: %s\n' "$id" "$n" "$cells" "$genes"
done
