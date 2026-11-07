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

### 2.3 Read every dataset into an AnnData of raw counts

Environments on the server (the Python one from the same `environment.yml`
as the laptop, the R one from `environment-r.yml`, used only to export the
authors' Seurat and RDS objects). Run them one after the other: two conda
solves at once share one repodata cache and the second fails with
"database is locked".

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && cd ~/Projects/pnatlas-dev && \
  conda env create -f environment.yml -n pnatlas -y && \
  conda activate pnatlas && pip install -e . && \
  python -c "import torch, scvi, jax; print(torch.__version__, torch.cuda.is_available(), scvi.__version__, jax.devices())" && \
  conda env create -f environment-r.yml -n pnatlas-r -y'
```

Prints `2.12.0 True 1.4.2 [CudaDevice(id=0)]`: the CUDA build of PyTorch is
picked on the server, the CPU build on the laptop.

Unpack the GEO RAW tars (the per-sample files inside are what the builder
reads; the tars stay as the downloaded copy):

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev/data/raw && \
  for d in nguyen2021_hdrg/GSE168243 jung2023_hdrg/GSE201586 heming2025_sural_snrna/GSE285983; do \
    ds=${d%/*}; acc=${d#*/}; mkdir -p $ds/extracted && tar -xf $ds/${acc}_RAW.tar -C $ds/extracted; done'
```

Sample to donor for SPARC 465 comes from its SDS `samples.xlsx` (library
`sam-hSG-lib-1` pools three donors); it is written into
`config/reference_inputs.yaml` as `donor_map`. The sheet is read with the
same stdlib parser as in step 1.5, with inline strings handled.

Export the three R objects (Seurat or dgCMatrix) as matrix market plus CSV,
with the memory guard running (the painseq object is the largest):

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && cd ~/Projects/pnatlas-dev && \
  ~/bin/memguard.sh start logs/memguard.log "build_anndata|rds_to_mtx|Rscript"; \
  conda activate pnatlas-r && \
  for p in bhuiyan2025_hdrg/seurat_object_preprint_release.Rds \
           jung2023_hdrg/GSE201586_humanallcelltypes.rds.gz \
           yang2022_htg/GSE197289_snRNA-seq_human_raw_counts.RDS.gz; do \
    d=${p%%/*}; nice -n 15 Rscript workflow/scripts/rds_to_mtx.R data/raw/$p data/interim/rds/$d > logs/rds_to_mtx_$d.log 2>&1; \
    cat data/interim/rds/$d/object_info.txt; done'
```

Build one h5ad per dataset (raw counts, obs columns dataset, sample, donor,
author_label, author_label_fine and the authors' metadata as `author_*`),
then the summary table:

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  for d in nguyen2021_hdrg jung2023_hdrg yu2024_singlesoma zhang_hdrg_glia yang2022_htg \
           yang2025_hsg_multiome heming2025_sural_snrna bhuiyan2025_hdrg; do \
    nice -n 15 python workflow/scripts/build_anndata.py $d --data-root data \
      --out data/interim/reference/$d.h5ad | tail -1; done && \
  mkdir -p results/reference && \
  python workflow/scripts/anndata_summary.py data/interim/reference > results/reference/anndata_summary.tsv'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/reference/anndata_summary.tsv" results/reference/
```

The donor map for Jung 2023 comes from the GEO sample characteristics:

```bash
python3 -I workflow/scripts/geo_samples.py GSE168243 GSE201586 GSE249746 GSE189501 GSE197289 GSE285983 \
  > data/metadata/geo_samples.tsv
```

Check: `results/reference/anndata_summary.tsv` has eight rows and the cell
counts match the analysis log (2.3).

### 2.4 QC and the cell type map

Thresholds are fixed in the `qc` block of `config/config.yaml` before the
first run: 200 genes minimum, 5 percent mitochondrial counts for nuclei and
20 percent for whole cells, scrublet per sample with expected rate 0.06 on
samples of at least 500 cells. Author-annotated datasets keep every cell and
carry the flags; the three unannotated datasets are filtered.

`config/celltype_map.tsv` is written by hand from the author label lists
(one row per dataset and author label, with class, state and fine level);
the label lists come from the h5ad files:

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  python - <<EOF
import anndata as ad, glob
for f in sorted(glob.glob("data/interim/reference/*.h5ad")):
    o = ad.read_h5ad(f, backed="r").obs
    print("===", o["dataset"].iloc[0]); print(o["author_label_fine"].value_counts().to_string())
EOF'
```

Run the QC, one dataset at a time, and collect the summary lines:

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  ~/bin/memguard.sh start logs/memguard.log "qc_dataset|build_anndata|Rscript"; \
  echo -e "dataset\tcells_in\tcells_out\tflag_low_genes\tflag_high_mt\tflag_doublet\tmedian_counts\tmedian_genes\tmedian_pct_mt\tpolicy" \
    > results/reference/qc_summary.tsv; \
  for d in nguyen2021_hdrg yu2024_singlesoma zhang_hdrg_glia yang2022_htg yang2025_hsg_multiome \
           jung2023_hdrg bhuiyan2025_hdrg heming2025_sural_snrna; do \
    nice -n 15 python workflow/scripts/qc_dataset.py $d --in data/interim/reference/$d.h5ad \
      --out data/interim/reference_qc/$d.h5ad 2> logs/qc_$d.log | tail -1 >> results/reference/qc_summary.tsv; done'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/reference/qc_summary.tsv" results/reference/
```

Check: eight rows in `results/reference/qc_summary.tsv`; every author label
in the h5ad files has a row in `config/celltype_map.tsv` (checked in 2.5
when the map is applied).
