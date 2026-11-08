#!/usr/bin/env bash
# Map every Visium section of the given datasets with cell2location, one
# section after another. A section whose checks.tsv exists is skipped, and
# a section another loop has claimed (its .claim directory exists) is
# skipped too, so several loops can run side by side on disjoint or
# overlapping dataset lists and a loop can be relaunched. A section that
# fails loses its claim so that a relaunch retries it.
# Usage: workflow/scripts/run_c2l_sections.sh <data_root> <signatures.csv.gz> <markers.tsv> <dataset> [<dataset> ...]
set -uo pipefail
root=$1; sig=$2; markers=$3; shift 3
mkdir -p logs
for ds in "$@"; do
  for sec in "$root"/interim/spatial/"$ds"/*.h5ad; do
    name=$(basename "$sec" .h5ad)
    out="$root/interim/spatial/c2l/$ds/$name"
    [[ -f "$out/checks.tsv" ]] && { echo "skip $ds $name (done)"; continue; }
    mkdir -p "$out"
    mkdir "$out/.claim" 2>/dev/null || { echo "skip $ds $name (claimed)"; continue; }
    nice -n 15 python workflow/scripts/c2l_map.py --section "$sec" --signatures "$sig" --markers "$markers" \
      --out-dir "$out" > "logs/c2l_map_${ds}_${name}.log" 2>&1
    rc=$?
    echo "EXIT $ds $name $rc"
    (( rc != 0 )) && rmdir "$out/.claim" 2>/dev/null
  done
done
echo ALLDONE
