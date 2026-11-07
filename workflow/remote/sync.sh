#!/usr/bin/env bash
# Copy the committed code, config and metadata to the server and record the
# git commit there. Large data never goes through this script.
# Usage: REMOTE=user@host REMOTE_DIR='~/Projects/pnatlas-dev' workflow/remote/sync.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${REMOTE:?set REMOTE=user@host}"
REMOTE_DIR=${REMOTE_DIR:-'~/Projects/pnatlas-dev'}
if [[ -n "$(git status --porcelain -- src workflow config/config.yaml config/datasets.yaml pyproject.toml environment.yml)" ]]; then
  echo "uncommitted changes in code; commit first so COMMIT is meaningful" >&2; exit 1
fi
ssh -o BatchMode=yes "$REMOTE" "mkdir -p $REMOTE_DIR/data/metadata $REMOTE_DIR/results $REMOTE_DIR/logs"
rsync -a --delete src workflow tests pyproject.toml environment.yml "$REMOTE:$REMOTE_DIR/"
rsync -a config/config.yaml config/datasets.yaml config/paths.example.yaml "$REMOTE:$REMOTE_DIR/config/"
rsync -a data/metadata/ "$REMOTE:$REMOTE_DIR/data/metadata/"
rsync -a data/README.md "$REMOTE:$REMOTE_DIR/data/"
[[ -f data/MANIFEST.tsv ]] && rsync -a data/MANIFEST.tsv "$REMOTE:$REMOTE_DIR/data/"
git rev-parse HEAD | ssh -o BatchMode=yes "$REMOTE" "cat > $REMOTE_DIR/COMMIT"
echo "synced $(git rev-parse --short HEAD) to $REMOTE:$REMOTE_DIR"
