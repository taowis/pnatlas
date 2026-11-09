#!/usr/bin/env bash
# Move the project's large data directories (data/raw, data/interim) from the
# project checkout on the server to a directory on the server's data array and
# leave symbolic links behind, so every command keeps using --data-root data.
# Runs on the server. Each directory is copied, the copy is compared with the
# source by checksum, and only then is the source removed and replaced by the
# link; a README is written in the target directory.
# Usage (on the server, in the project checkout):
#   nice -n 15 workflow/remote/relocate_data.sh <target data dir> [dir ...]
# e.g. workflow/remote/relocate_data.sh /path/on/array/pnatlas-dev/data raw interim
set -euo pipefail
target=${1:?target data directory on the array}
shift
dirs=("$@")
[[ ${#dirs[@]} -eq 0 ]] && dirs=(raw interim)
mkdir -p "$target"
if [[ ! -f "$target/../README.md" ]]; then
  cat > "$target/../README.md" <<EOF
# pnatlas-dev data

Raw downloads and intermediate files of the pnatlas-dev analysis (human
peripheral nervous system reference), moved here from the project checkout
in the home directory because the system disk is small. The checkout keeps
symbolic links data/raw and data/interim pointing here, so the code and the
runbook are unchanged. Every raw file is listed with its source URL and md5
in the checkout's data/MANIFEST.tsv and data/README.md and can be fetched
again; the intermediate files are rebuilt by the workflow. Nothing here is
the only copy of anything. State and what to run next: the project README
and docs/runbook/ in the checkout.
EOF
fi
for d in "${dirs[@]}"; do
  src="data/$d"
  if [[ -L "$src" ]]; then echo "$src is already a link to $(readlink "$src")"; continue; fi
  [[ -d "$src" ]] || { echo "$src: no such directory"; continue; }
  echo "copying $src to $target/$d"
  rsync -a "$src/" "$target/$d/"
  echo "comparing by checksum"
  diff=$(rsync -anc --delete --out-format='%n' "$src/" "$target/$d/" | grep -v '/$' || true)
  if [[ -n "$diff" ]]; then echo "DIFFERENCES after copy of $src:"; echo "$diff"; exit 1; fi
  rm -rf "$src"
  ln -s "$target/$d" "$src"
  echo "done: $src -> $(readlink "$src")"
done
df -h . "$target" | tail -2
