# Runbook, stage 3: spatial mapping

Commands in run order. `$REPO`, `$DATA_ROOT` and `$REMOTE` are set as in
[README.md](README.md) and [02_reference.md](02_reference.md). Every step
is run with the `pnatlas` environment active unless stated otherwise.

### 3.1 Environment additions

cell2location and its dependencies are added to the `pip:` section of
`environment.yml` and installed into the existing environments on both
machines. cell2location 0.1.5 does not import with the pyro that
scvi-tools pulls in (1.9.2), so pyro is pinned to 1.9.0:

```bash
pip install --no-cache-dir "cell2location==0.1.5" openpyxl
pip install --no-cache-dir "pyro-ppl==1.9.0"
python -c "import cell2location, scvi, torch, pyro; print(cell2location.__version__, scvi.__version__, torch.__version__, pyro.__version__)"
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && \
  pip install --no-cache-dir "cell2location==0.1.5" openpyxl && pip install --no-cache-dir "pyro-ppl==1.9.0"'
conda env export -n pnatlas --no-builds | grep -v "^prefix:" > docs/runbook/environment_lock.yml
```

Check: the import line prints `0.1.5 1.4.2 2.12.0 1.9.0` (laptop; the
server has the same cell2location, scvi-tools and pyro).

### 3.2 Download the spatial layer on the server

Build the file list on the laptop (GEO supplementary listings and the two
public Box folders of the UT Dallas apps, whose file identifiers are
parsed from the folder pages) and the sample sheets of the GEO series,
and commit them:

```bash
cd "$REPO"
python3 -I workflow/scripts/stage3_list_files.py > data/metadata/stage3_files.tsv
awk -F'\t' 'NR>1{s+=$5} END{print "total GB", s/1073741824}' data/metadata/stage3_files.tsv
python3 -I workflow/scripts/geo_samples.py GSE285984 GSE295206 GSE243073 GSE215994 GSE273557 > data/metadata/geo_samples_spatial.tsv
```

29 files, 21.7 GB: the RAW tars of the five GEO series (Heming 2025 Xenium
9.9 GB, Shiers 2025 DRG Visium 9.0 GB, Yu 2024 Xenium transcripts 0.8 GB,
Vroman 2023 and Zhang 2024 DRG Visium), the Yu 2024 neuron matrix and
its metadata, the Shiers sample metadata, and from Box the 16 Space
Ranger H5 files, the sample sheet and the combined Seurat object of the
Tavares-Ferreira 2022 DRG Visium data and the annotated Seurat object and
sample sheet of the DPN sural nerve Visium data. The Box H and E images
and Loupe files (37 GB) are not taken.

Sync and launch on the server (log `logs/stage3_download.log`):

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && \
  setsid -f nohup python3 -I workflow/scripts/download_files.py data/metadata/stage3_files.tsv \
    --data-root data --manifest data/MANIFEST.tsv \
    --fetch "workflow/scripts/download_files.py data/metadata/stage3_files.tsv" \
    > logs/stage3_download.log 2>&1 < /dev/null'
```

Wait for it, verify the checksums and bring the manifest back:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && grep -cE "^done" logs/stage3_download.log && \
  python3 -I workflow/scripts/download_files.py data/metadata/stage3_files.tsv --data-root data --verify'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/MANIFEST.tsv" data/MANIFEST.tsv
```

Check: the verify run prints `OK` for all 29 paths; the manifest rows
(url, md5, size, access date) are then copied into the dataset sections
of `data/README.md`.
