#!/usr/bin/env bash
# Pack the staged archive into the files uploaded to Zenodo, which takes at
# most 100 files per record: the reference object alone, one tar.gz per
# component, and the README, data README, manifest and inner MD5 list as
# they are. Writes MD5SUMS_bundles next to the bundles.
# Usage: workflow/release/bundle_zenodo.sh <staged_dir> <out_dir> [version]
set -euo pipefail
staged=${1:?staged directory}; out=${2:?output directory}; version=${3:-0.1.0}
mkdir -p "$out"
tarball() {  # name, then the paths inside the staged directory
  local name=$1; shift
  tar -C "$staged" --no-xattrs -cf - "$@" | gzip -n -6 > "$out/$name.tar.gz"
  echo "$name.tar.gz: $(tar -tzf "$out/$name.tar.gz" | grep -vc '/$') files"
}
cp "$staged/reference/pnatlas_reference_v$version.h5ad" "$out/"
tarball reference_tables reference/pnatlas_reference_cells.tsv.gz reference/state_specificity.tsv.gz \
  reference/state_expression_cpm.tsv.gz reference/state_pseudobulk_cells.tsv
tarball spatial spatial
tarball deconvolution deconvolution
tarball results results
for f in README.md data_README.md MANIFEST.tsv MD5SUMS; do cp "$staged/$f" "$out/"; done
(cd "$out" && md5sum pnatlas_reference_v"$version".h5ad *.tar.gz README.md data_README.md MANIFEST.tsv MD5SUMS > MD5SUMS_bundles)
ls -l "$out"
echo "BUNDLED $(ls "$out" | wc -l) files"
