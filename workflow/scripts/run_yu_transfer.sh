#!/usr/bin/env bash
# Rerun the signature transfer of the two Yu 2024 datasets (runbook 3.10)
# after a change of their candidate states, then collect the checks again
# and redraw the H and E anatomy panels of the all-cell dataset. The
# section objects, panel genes and panel markers of 3.9 are reused. Each
# step prints a "<STEP> <name> <exit code>" line; the run ends with YUTRANSFERDONE.
# Runs in the project checkout with the pnatlas environment active.
# Usage: workflow/scripts/run_yu_transfer.sh <data_root> [results dir]
set -uo pipefail
root=${1:?data root}
res=${2:-results/spatial}
mkdir -p logs
for d in yu2024_hdrg_xenium yu2024_hdrg_xenium_cells; do
  nice -n 15 python workflow/scripts/panel_transfer.py "$d" --signatures "$root/interim/spatial/signatures/reference_signatures.csv.gz" \
    --annotation "$root/interim/reference_annotation.tsv.gz" --markers "$res/$d/panel_markers.tsv" \
    --sections-dir "$root/interim/spatial/$d" --out-dir "$root/interim/spatial/transfer/$d" \
    --results-dir "$res/$d" > "logs/panel_transfer_$d.log" 2>&1
  echo "TRANSFER $d $?"
done
python workflow/scripts/collect_spatial_checks.py --c2l-dir "$root/interim/spatial/c2l" --results-dir "$res" \
  --sections "$res/section_summary.tsv" > logs/collect_spatial_checks.log 2>&1
echo "COLLECT $?"
d=yu2024_hdrg_xenium_cells
nice -n 15 python workflow/scripts/anatomy_panels.py "$d" --data-root "$root" --out-dir "$res/anatomy/$d" --width-px 1600 \
  > "logs/anatomy_panels_$d.log" 2>&1
echo "PANELS $?"
echo YUTRANSFERDONE
