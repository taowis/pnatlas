# Runbook, stage 2: joint ganglia and nerve reference

Commands in run order. `$REPO` and `$DATA_ROOT` are set as in
[README.md](README.md). Every step is run with the `pnatlas` environment
active unless stated otherwise.

### 2.1 Environment and machine config

The environment build is in the Environment section of [README.md](README.md).
Machine config for this run:

```bash
cd "$REPO"
cp config/paths.example.yaml config/paths.yaml
# data_root: a directory on the external drive with at least 300 GB free
# pnprank_root: the pnprank-dev clone (tag submitted-human-genetics)
mkdir -p "$DATA_ROOT"
conda activate pnatlas
python -c "import torch; print('mps', torch.backends.mps.is_available())"
```

Prints `mps True` on the laptop.

### 2.2 Download the reference layer on the server

The server is reached as `$REMOTE` (user@host, set in the shell, never in a
file). Code goes over with the sync script, which refuses uncommitted code and
writes the commit hash to `COMMIT` on the server.

Build the file list on the laptop (GEO supplementary listings, the Pennsieve
file browser for SPARC 465, the painseq Dropbox link) and commit it:

```bash
cd "$REPO"
python3 -I workflow/scripts/stage2_list_files.py > data/metadata/stage2_files.tsv
awk -F'\t' 'NR>1{s+=$5} END{print "total GB", s/1073741824}' data/metadata/stage2_files.tsv
```

19 files, 11.3 GB: six GEO series (RAW tars, author metadata, the Jung 2023
annotated object), the SPARC 465 filtered matrix with its SDS sheets (ATAC
fragments left out, 17.8 GB), and the painseq Seurat object.

Sync and launch on the server (the job survives the ssh session; its log is
`logs/stage2_download.log` on the server):

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && mkdir -p data/raw logs && \
  setsid -f nohup python3 -I workflow/scripts/download_files.py data/metadata/stage2_files.tsv \
    --data-root data --manifest data/MANIFEST.tsv \
    --fetch "workflow/scripts/download_files.py data/metadata/stage2_files.tsv" \
    > logs/stage2_download.log 2>&1 < /dev/null'
```

Wait for it, then verify the checksums and bring the manifest back:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && grep -cE "^done" logs/stage2_download.log && \
  python3 -I workflow/scripts/download_files.py data/metadata/stage2_files.tsv --data-root data --verify'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/MANIFEST.tsv" data/MANIFEST.tsv
```

The verify run prints `OK` for all 19 paths. The manifest rows (url, md5,
size, access date) are then copied into the dataset sections of
`data/README.md`.
