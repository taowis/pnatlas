# Stage 7: portal, data archive and public copy

## 7.1 Static gene lookup site

`portal/build_portal.py` writes the site from the committed result tables
(no other input); `portal/index.html` is the page. The build is checked by
`tests/test_portal.py`, which builds into a temporary directory and reads
MPZ back. The site itself is not committed: the pages workflow
(`.github/workflows/pages.yml`) rebuilds and deploys it; it is run by
hand (`gh workflow run pages`) until Pages is set to "GitHub Actions" in
the repository settings, after which a push trigger on `main` is added.

```bash
conda activate pnatlas
python portal/build_portal.py                          # portal/site/
(cd portal/site && python -m http.server 8931 --bind 127.0.0.1)   # open http://127.0.0.1:8931/#gene/MPZ
ruff check src tests portal && pytest -q tests/test_portal.py
```

Check: the build prints `17154 genes in 557 shards, 31 states`; the page
for MPZ lists myelinating Schwann first (specificity 0.608) and the GTEx
column has myelinating Schwann at 0.88; the States and About views open
from the header; no console errors. A port that another local
application has registered a service worker on shows that application
instead of the site, so use a port no other project has used.

## 7.2 Citation and archive metadata

`CITATION.cff` and `.zenodo.json` carry the title, author with ORCID and
affiliation, MIT licence and keywords, in the layout of the PNPrank
release; `results/SNAPSHOT.txt` names the commit the result tables, the
portal and the data archive are built from.

## 7.3 Data archive

The cell-level reference object is written without the Bhuiyan 2025
cells (their portal forbids redistributing cell-level derivatives,
`config/datasets.yaml`); the state-level tables keep every dataset. The
staging script collects the reference, the spatial outputs (cell2location
abundances per section, panel transfers, the signature model's tables),
the deconvolution reference and mixtures, the result tables and the data
records into one directory with an MD5 list, on the server, for the
Zenodo record (7.5).

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && free -g | sed -n 2p && who | wc -l && nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader && \
  ~/bin/memguard.sh start logs/memguard.log "reference_for_release" && \
  setsid -f nohup nice -n 15 bash -c "source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && \
    workflow/release/stage_zenodo_data.sh data data/release/zenodo_v0.1.0 0.1.0" > logs/stage_zenodo.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && tail -3 logs/stage_zenodo.log && ~/bin/memguard.sh stop logs/memguard.log && \
  cat data/release/zenodo_v0.1.0/MD5SUMS | wc -l && du -sh data/release/zenodo_v0.1.0'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/release/zenodo_v0.1.0/MD5SUMS" data/metadata/zenodo_v0.1.0_md5sums.txt
```

Check: the log ends with `STAGED`; `reference_for_release.py` reports the
kept cells (644,681 minus the 47,226 Bhuiyan cells); the MD5 list is
copied into `data/metadata/` so the archive's content is recorded in the
repository before the upload.

## 7.4 Public copy

The build is rehearsed first so that the rewrite is known to be clean:

```bash
PRIVATE_PATTERN='<regex of the server name and address>' workflow/release/build_public.sh /tmp/pnatlas_public
```

Check: `attribution lines in history: 0`, `private names in any tree: 0`,
`manuscript files in any tree: 0`, and `built in ... (not pushed)`. The
same command with `--push` creates `taowis/pnatlas` and pushes `main`
and the tags; rerun after every merge to `main` that should reach the
public copy (the rewrite is deterministic, so the push is a fast
forward):

```bash
PRIVATE_PATTERN='<regex of the server name and address>' workflow/release/build_public.sh /tmp/pnatlas_public --push
gh repo view taowis/pnatlas --json visibility,url
```

## 7.5 Zenodo record

The staged archive (7.3) is copied from the server to the local data
root, verified against the MD5 list, and packed into the nine upload
files, because Zenodo takes at most 100 files per record:

```bash
rsync -a --stats "$REMOTE:~/Projects/pnatlas-dev/data/release/zenodo_v0.1.0/" "$DATA_ROOT/release/zenodo_v0.1.0/"
(cd "$DATA_ROOT/release/zenodo_v0.1.0" && md5sum -c MD5SUMS | grep -vc ': OK$')   # 0
workflow/release/bundle_zenodo.sh "$DATA_ROOT/release/zenodo_v0.1.0" "$DATA_ROOT/release/zenodo_v0.1.0_upload" 0.1.0
cp "$DATA_ROOT/release/zenodo_v0.1.0_upload/MD5SUMS_bundles" data/metadata/zenodo_v0.1.0_bundles_md5sums.txt
```

The draft record is created and filled with a personal access token
(scope `deposit:write`, made at
https://zenodo.org/account/settings/applications/tokens/new/ and saved
in `~/.config/zenodo/token`, which is outside the repository). The
metadata is `workflow/release/zenodo_metadata.json`; the DOI is reserved
with the draft and printed:

```bash
python workflow/release/zenodo_upload.py --dir "$DATA_ROOT/release/zenodo_v0.1.0_upload"
```

Check: every file line says `uploaded` or `present`, the last line is
`DRAFT READY <record id> <DOI>`, and the draft page lists the nine files
with the sizes of `MD5SUMS_bundles`. The reserved DOI is then written
into `CITATION.cff` (`identifiers`), `data/README.md` (stage 7 section)
and the README data availability line, and the record is published from
its draft page in the browser (the publish step registers the DOI and
cannot be undone, so it is the last step). The server copy of the
staged archive is deleted once the record is published and the local
copy verified:

```bash
ssh "$REMOTE" 'rm -rf ~/Projects/pnatlas-dev/data/release/zenodo_v0.1.0'
```

## 7.6 GitHub Pages

Pages is enabled on the public repository with the "GitHub Actions"
source, then the pages workflow is run there once by hand; after that
every push to `main` of the public copy that touches the site's inputs
redeploys it (the deploy job runs only in `taowis/pnatlas`, the build job
also runs in the private repository as a check):

```bash
gh api -X POST repos/taowis/pnatlas/pages -f build_type=workflow
gh workflow run pages -R taowis/pnatlas
gh run list -R taowis/pnatlas --workflow pages --limit 1
gh api repos/taowis/pnatlas/pages --jq '.html_url+" "+.status'
```

Check: the run ends with both jobs green and the site answers at the
`html_url` with the gene page for `#gene/MPZ`.
