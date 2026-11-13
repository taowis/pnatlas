#!/usr/bin/env bash
# Build the public copy of this repository: a rewritten clone of main (and its
# tags) without the manuscript, with no attribution trailers in any commit
# message, checked for private names. Without --push it only builds and
# reports; with --push it creates the public repository if needed and pushes.
# Usage: PRIVATE_PATTERN='<regex of names that must not appear>' \
#        workflow/release/build_public.sh <work_dir> [--push]
# PUBLIC_REPO defaults to taowis/pnatlas.
set -euo pipefail
work=${1:?work directory}; push=${2:-}
: "${PRIVATE_PATTERN:?set PRIVATE_PATTERN to the regex of names that must not appear}"
PUBLIC_REPO=${PUBLIC_REPO:-taowis/pnatlas}
repo=$(git rev-parse --show-toplevel)
mkdir -p "$work"
[[ -x "$work/fr/bin/git-filter-repo" ]] || { python3 -m venv "$work/fr" && "$work/fr/bin/pip" -q install git-filter-repo; }
rm -rf "$work/pub"
git clone -q --no-local --single-branch -b main "$repo" "$work/pub"
cd "$work/pub"
git fetch -q origin 'refs/tags/*:refs/tags/*'
"$work/fr/bin/git-filter-repo" --force --quiet --invert-paths --path manuscript/ --message-callback '
import re
lines = [l for l in message.split(b"\n") if not re.match(rb"Co-Authored-By:", l, re.I)]
return b"\n".join(lines).rstrip(b"\n") + b"\n"
'
echo "commits: $(git rev-list --count HEAD)"
echo "attribution lines in history: $(git log --format=%B | grep -c -i -E 'co-authored-by|generated with' || true)"
echo "private names in any tree: $(git log --all -p | grep -c -i -E "$PRIVATE_PATTERN" || true)"
echo "manuscript files in any tree: $(git log --all --name-only --format= | grep -c '^manuscript/' || true)"
git ls-files | grep -E '^(private_notes/|CLAUDE.md|\.claude/)' && { echo "private paths tracked" >&2; exit 1; }
if [[ "$push" == "--push" ]]; then
  gh repo view "$PUBLIC_REPO" >/dev/null 2>&1 || gh repo create "$PUBLIC_REPO" --public \
    --description "Multiscale reference of human peripheral nervous system cell states: ganglia and nerve single-cell, spatial, GTEx bulk and neuropathy genetics"
  git remote add public "https://github.com/$PUBLIC_REPO.git"
  git -c credential.helper='!gh auth git-credential' push -q public main
  git -c credential.helper='!gh auth git-credential' push -q public --tags
  echo "pushed to $PUBLIC_REPO"
else
  echo "built in $work/pub (not pushed)"
fi
