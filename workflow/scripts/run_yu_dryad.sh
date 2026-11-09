#!/usr/bin/env bash
# Yu 2024 from the Dryad bundles (runbook 3.9): take the members the analysis
# needs out of each rar (cells table, cell by gene matrix, post-Xenium H and E
# image with its alignment, run metadata), build the per-section objects of
# both Yu datasets (the hand-segmented neurons with corrected coordinates and
# every segmented cell), rebuild the section summary, compute the panel
# markers and run the signature transfer for the all-cell dataset, collect
# the checks and draw the H and E anatomy panels. Each step prints a
# "<STEP> <name> <exit code>" line; the run ends with YUDONE.
# Runs in the project checkout with the pnatlas environment active.
# Usage: workflow/scripts/run_yu_dryad.sh <data_root> [results dir]
set -uo pipefail
root=${1:?data root}
res=${2:-results/spatial}
dryad=$root/raw/yu2024_hdrg_xenium/dryad
mkdir -p logs
for n in 1 2 3 4; do
  nice -n 15 python -I workflow/scripts/extract_members.py \
    "$dryad/output-XETG00171__0018220__Region_${n}__20240207__003742.rar" "$dryad/extracted/Region_$n" --strip 1 \
    '*/cells.csv.gz' '*/cells.parquet' '*/cell_feature_matrix.h5' "*/DRG_area_$n.ome.tif" \
    "*/DRG_area_${n}_alignment_files/*" '*/experiment.xenium' '*/gene_panel.json' '*/metrics_summary.csv' \
    > "logs/extract_dryad_Region_$n.log" 2>&1
  echo "EXTRACT Region_$n $?"
done
for d in yu2024_hdrg_xenium yu2024_hdrg_xenium_cells; do
  nice -n 15 python workflow/scripts/build_spatial.py "$d" --data-root "$root" --out-dir "$root/interim/spatial/$d" \
    > "$res/sections_$d.tsv" 2> "logs/build_spatial_$d.log"
  echo "BUILD $d $?"
done
cat "$res"/sections_*.tsv | awk 'NR==1||!/^dataset/' > "$res/section_summary.tsv"
d=yu2024_hdrg_xenium_cells
mkdir -p "$res/$d"
python -c "import anndata as ad, glob, sys; a = ad.read_h5ad(sorted(glob.glob(sys.argv[1] + '/*.h5ad'))[0], backed='r'); print('\n'.join(a.var_names))" \
  "$root/interim/spatial/$d" > "$res/$d/panel_genes.txt"
nice -n 15 python workflow/scripts/reference_markers.py --counts "$root/interim/reference_merged.h5ad" \
  --annotation "$root/interim/reference_annotation.tsv.gz" --genes "$res/$d/panel_genes.txt" \
  --out "$res/$d/panel_markers.tsv" > "logs/panel_markers_$d.log" 2>&1
echo "MARKERS $d $?"
nice -n 15 python workflow/scripts/panel_transfer.py "$d" --signatures "$root/interim/spatial/signatures/reference_signatures.csv.gz" \
  --annotation "$root/interim/reference_annotation.tsv.gz" --markers "$res/$d/panel_markers.tsv" \
  --sections-dir "$root/interim/spatial/$d" --out-dir "$root/interim/spatial/transfer/$d" \
  --results-dir "$res/$d" > "logs/panel_transfer_$d.log" 2>&1
echo "TRANSFER $d $?"
python workflow/scripts/collect_spatial_checks.py --c2l-dir "$root/interim/spatial/c2l" --results-dir "$res" \
  --sections "$res/section_summary.tsv" > logs/collect_spatial_checks.log 2>&1
echo "COLLECT $?"
nice -n 15 python workflow/scripts/anatomy_panels.py "$d" --data-root "$root" --out-dir "$res/anatomy/$d" --width-px 1600 \
  > "logs/anatomy_panels_$d.log" 2>&1
echo "PANELS $?"
echo YUDONE
