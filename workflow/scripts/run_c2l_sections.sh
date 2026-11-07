#!/usr/bin/env bash
# Map every Visium section of the given datasets with cell2location, one
# section after another (one model on the GPU at a time). A section whose
# checks.tsv exists is skipped, so the loop can be relaunched.
# Usage: workflow/scripts/run_c2l_sections.sh <data_root> <signatures.csv.gz> <markers.tsv> <dataset> [<dataset> ...]
set -uo pipefail
root=$1; sig=$2; markers=$3; shift 3
mkdir -p logs
for ds in "$@"; do
  for sec in "$root"/interim/spatial/"$ds"/*.h5ad; do
    name=$(basename "$sec" .h5ad)
    out="$root/interim/spatial/c2l/$ds/$name"
    [[ -f "$out/checks.tsv" ]] && { echo "skip $ds $name"; continue; }
    nice -n 15 python workflow/scripts/c2l_map.py --section "$sec" --signatures "$sig" --markers "$markers" \
      --out-dir "$out" > "logs/c2l_map_${ds}_${name}.log" 2>&1
    echo "EXIT $ds $name $?"
  done
done
echo ALLDONE
