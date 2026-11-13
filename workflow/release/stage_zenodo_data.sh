#!/usr/bin/env bash
# Collect the data files of the archive into one directory, with checksums and
# a README, ready to upload to Zenodo by hand.
# Usage: workflow/release/stage_zenodo_data.sh <data_root> <out_dir> [version]
# Run from the repository root on the machine that holds the data (the server).
set -euo pipefail
root=${1:?data root}; out=${2:?output directory}; version=${3:-0.1.0}
mkdir -p "$out/reference" "$out/spatial" "$out/deconvolution" "$out/genetics" "$out/results"

echo "REFERENCE"
python workflow/release/reference_for_release.py \
  --merged "$root/interim/reference_merged.h5ad" --embeddings "$root/interim/reference_embeddings.h5ad" \
  --annotation "$root/interim/reference_annotation.tsv.gz" --exclude bhuiyan2025_hdrg \
  --out "$out/reference/pnatlas_reference_v$version.h5ad" --cells-out "$out/reference/pnatlas_reference_cells.tsv.gz"
cp "$root/interim/genetics/specificity.tsv.gz" "$out/reference/state_specificity.tsv.gz"
cp "$root/interim/genetics/expression_cpm.tsv.gz" "$out/reference/state_expression_cpm.tsv.gz"
cp "$root/interim/genetics/pseudobulk_cells.tsv" "$out/reference/state_pseudobulk_cells.tsv"

echo "SPATIAL"
rsync -a --exclude .claim "$root/interim/spatial/c2l/" "$out/spatial/cell2location/"
rsync -a "$root/interim/spatial/transfer/" "$out/spatial/panel_transfer/"
rsync -a --exclude model "$root/interim/spatial/signatures/" "$out/spatial/signatures/"

echo "DECONVOLUTION"
for f in reference_counts.mtx reference_cells.tsv genes.txt donors.tsv mixtures_counts.tsv.gz mixtures_truth.tsv; do
  cp "$root/interim/deconv/$f" "$out/deconvolution/"
done
gzip -f "$out/deconvolution/reference_counts.mtx"

echo "RESULTS"
rsync -a --exclude 'anatomy' results/ "$out/results/"
cp data/MANIFEST.tsv data/README.md "$out/"

echo "CHECKSUMS"
(cd "$out" && find . -type f ! -name MD5SUMS -print0 | sort -z | xargs -0 md5sum > MD5SUMS)
du -sh "$out"
echo "STAGED"
