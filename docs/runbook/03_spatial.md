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

28 files, 21.7 GB: the RAW tars of the five GEO series (Heming 2025 Xenium
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

Check: the verify run prints `OK` for all 28 paths; the manifest rows
(url, md5, size, access date) are then copied into the dataset sections
of `data/README.md`.

## 3.3 Section inputs

List the tar members, extract the matrices, cell tables and positions
(not the morphology images, the Loupe files or the raw matrices) and
look at what each format holds:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev/data/raw && \
  for t in heming2025_sural_xenium/GSE285984_RAW.tar yu2024_hdrg_xenium/GSE273557_RAW.tar \
           zhang2024_hdrg_visium/GSE243073_RAW.tar vroman2023_hdrg_visium/GSE215994_RAW.tar \
           shiers2025_dpn_drg/GSE295206_RAW.tar; do echo "== $t"; tar -tvf $t | awk "{print \$3, \$6}"; done'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev/data/raw && \
  mkdir -p heming2025_sural_xenium/extracted zhang2024_hdrg_visium/extracted \
           vroman2023_hdrg_visium/extracted shiers2025_dpn_drg/extracted && \
  tar -xf heming2025_sural_xenium/GSE285984_RAW.tar -C heming2025_sural_xenium/extracted \
      --wildcards "*_cell_feature_matrix.h5" "*_cells.parquet.gz" && \
  gunzip -f heming2025_sural_xenium/extracted/*.parquet.gz && \
  tar -xf zhang2024_hdrg_visium/GSE243073_RAW.tar -C zhang2024_hdrg_visium/extracted \
      --wildcards "*_filtered_feature_bc_matrix.tar.gz" && \
  (cd zhang2024_hdrg_visium/extracted && for f in *.tar.gz; do \
     mkdir -p ${f%.tar.gz} && tar -xzf $f -C ${f%.tar.gz} --strip-components=1; done) && \
  tar -xf vroman2023_hdrg_visium/GSE215994_RAW.tar -C vroman2023_hdrg_visium/extracted \
      --wildcards "*_filtered_feature_bc_matrix.h5" "*_tissue_hires_image.png.gz" && \
  tar -xf shiers2025_dpn_drg/GSE295206_RAW.tar -C shiers2025_dpn_drg/extracted \
      --wildcards "*_filtered_feature_bc_matrix.h5" "*tissue_positions*" "*NgnNeuSelections*" "*sample[0-9]*.csv.gz" && \
  gunzip -f shiers2025_dpn_drg/extracted/*.csv.gz vroman2023_hdrg_visium/extracted/*.png.gz && \
  for f in zhang2024_hdrg_visium/extracted/*.tar.gz; do tar -tzvf $f | awk "{print \$3, \$6}"; done && \
  md5sum shiers2025_dpn_drg/extracted/*filtered*.h5 | awk "{print \$1}" | sort > /tmp/s.md5 && \
  md5sum tavaresferreira2022_hdrg_visium/*.h5 | awk "{print \$1}" | sort > /tmp/t.md5 && \
  comm -12 /tmp/s.md5 /tmp/t.md5 | wc -l'
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev/data/raw && \
  python -I -c "
import pandas as pd, h5py, numpy as np, scanpy as sc, glob
d = pd.read_parquet(\"heming2025_sural_xenium/extracted/GSM8714492_S01_cells.parquet\"); print(d.shape); print(d.dtypes.to_string())
f = h5py.File(\"heming2025_sural_xenium/extracted/GSM8714492_S01_cell_feature_matrix.h5\")
print(np.unique(f[\"matrix/features/feature_type\"][:].astype(str), return_counts=True))
for p in sorted(glob.glob(\"*/extracted/*.h5\")) + sorted(glob.glob(\"tavaresferreira2022_hdrg_visium/*.h5\")):
    a = sc.read_10x_h5(p, gex_only=False); print(p.split(\"/\")[-1], a.shape, a.var[\"feature_types\"].unique().tolist(), int(a.X.sum()))
"'
```

The sample sheets that are not on GEO: the Shiers 2025 docx (samples,
donors, batches), the Tavares-Ferreira sample sheet and the DPN nerve
sample sheet are small and read on the laptop:

```bash
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/raw/shiers2025_dpn_drg/GSE295206_metadata_nageotteSamples.docx.gz" \
  "$REMOTE:~/Projects/pnatlas-dev/data/raw/tavaresferreira2022_hdrg_visium/README_sampleIDs.xlsx" \
  "$REMOTE:~/Projects/pnatlas-dev/data/raw/utdallas_dpn_nerves/sample_metadata.csv" /tmp/stage3_small/
gunzip -f /tmp/stage3_small/GSE295206_metadata_nageotteSamples.docx.gz
python3 -I -c "
import zipfile, re
x = zipfile.ZipFile('/tmp/stage3_small/GSE295206_metadata_nageotteSamples.docx').read('word/document.xml').decode()
print(re.sub(r'<[^>]+>', '', re.sub(r'</w:tc>', '\t', re.sub(r'</w:p>', '\n', x))))"
python3 -I -c "
import pandas as pd
for n, d in pd.read_excel('/tmp/stage3_small/README_sampleIDs.xlsx', sheet_name=None, header=None).items(): print(n); print(d.to_string())"
cat /tmp/stage3_small/sample_metadata.csv
```

The per-section inputs are written to `config/spatial_inputs.yaml` and
read by `workflow/scripts/build_spatial.py`, one h5ad per section.

Export the two Seurat objects (counts, metadata, coordinates with scale
factors per image) in the R environment, under the memory guard, and
build the sections of the five GEO datasets:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && cd ~/Projects/pnatlas-dev && \
  ~/bin/memguard.sh stop; ~/bin/memguard.sh start logs/memguard.log "rds_spatial_to_mtx|Rscript|build_spatial" && \
  conda activate pnatlas-r && mkdir -p data/interim/rds && \
  setsid -f nohup bash -c "for p in tavaresferreira2022_hdrg_visium/drg.combined.rds \
      utdallas_dpn_nerves/merged_obj_harmony_layerintegration_resol0.5_annotated_final.rds; do d=\${p%%/*}; \
    nice -n 15 Rscript workflow/scripts/rds_spatial_to_mtx.R data/raw/\$p data/interim/rds/\$d \
      > logs/rds_spatial_\$d.log 2>&1; done; echo ALLDONE >> logs/rds_spatial_utdallas_dpn_nerves.log" \
    > /dev/null 2>&1 < /dev/null'
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  mkdir -p data/interim/spatial results/spatial && \
  setsid -f nohup bash -c "for d in heming2025_sural_xenium yu2024_hdrg_xenium zhang2024_hdrg_visium \
      vroman2023_hdrg_visium shiers2025_dpn_drg; do \
    nice -n 15 python workflow/scripts/build_spatial.py \$d --data-root data --out-dir data/interim/spatial/\$d \
      > results/spatial/sections_\$d.tsv 2> logs/build_spatial_\$d.log; echo \"EXIT \$d \$?\" >> logs/build_spatial.log; done" \
    > /dev/null 2>&1 < /dev/null'
```

Check: `logs/build_spatial.log` has five `EXIT <dataset> 0` lines and both
`logs/rds_spatial_<dataset>.log` end with the object summary
(`spots_with_coordinates` equal to `spots`).

Add the two Box datasets to the config from the export, rebuild Shiers
(conditions added) and build the two, then collect one summary table:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  for d in shiers2025_dpn_drg tavaresferreira2022_hdrg_visium utdallas_dpn_nerves; do \
    nice -n 15 python workflow/scripts/build_spatial.py $d --data-root data --out-dir data/interim/spatial/$d \
      > results/spatial/sections_$d.tsv 2> logs/build_spatial_$d.log; echo "EXIT $d $?"; done && \
  cat results/spatial/sections_*.tsv | awk "NR==1||!/^dataset/" > results/spatial/section_summary.tsv && \
  rm results/spatial/sections_*.tsv && du -sh data/interim/spatial'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/spatial/section_summary.tsv" results/spatial/
```

Check: 56 sections (8 Heming, 4 Yu, 6 Zhang, 2 Vroman, 16 Shiers, 14
Tavares-Ferreira, 6 DPN nerve); coordinates present for all but Zhang,
Vroman and Tavares-Ferreira.


## 3.4 Reference signatures for cell2location

Parameters are the `cell2location` block of `config/config.yaml`. Probe
first (20,000 cells, 5 epochs, under the memory guard, peak memory from
`time -v`), then the full run detached on the GPU:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader && who && \
  ~/bin/memguard.sh start logs/memguard.log "c2l_reference_signatures|c2l_map" && \
  mkdir -p data/interim/spatial/probe && \
  nice -n 15 /usr/bin/time -v python workflow/scripts/c2l_reference_signatures.py \
    --counts data/interim/reference_merged.h5ad --annotation data/interim/reference_annotation.tsv.gz \
    --out-dir data/interim/spatial/probe --subsample 20000 --max-epochs 5 > logs/c2l_signatures_probe.log 2>&1; \
  grep -E "^cells|^loaded|^genes|^elbo|^signatures|Maximum resident|Elapsed" logs/c2l_signatures_probe.log'
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  setsid -f nohup nice -n 15 /usr/bin/time -v python workflow/scripts/c2l_reference_signatures.py \
    --counts data/interim/reference_merged.h5ad --annotation data/interim/reference_annotation.tsv.gz \
    --out-dir data/interim/spatial/signatures > logs/c2l_signatures.log 2>&1 < /dev/null'
```

When it ends: copy the small outputs back, remove the probe directory.

```bash
mkdir -p results/spatial/signatures
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/interim/spatial/signatures/{signature_cells.tsv,training_history.tsv}" results/spatial/signatures/
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/interim/spatial/signatures/reference_signatures.csv.gz" "$DATA_ROOT/interim/spatial/"
ssh "$REMOTE" 'rm -r ~/Projects/pnatlas-dev/data/interim/spatial/probe'
```

Check: `logs/c2l_signatures.log` ends with `signatures (<genes>, 31) written`;
the ELBO in `training_history.tsv` has flattened over the last 50 epochs.
