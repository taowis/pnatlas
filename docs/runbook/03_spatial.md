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

The reference markers for check (a) (parameters: the `spatial_checks`
block of the config) are computed on the server while the signatures
train (CPU only, about 10 GB):

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  setsid -f nohup nice -n 15 /usr/bin/time -v python workflow/scripts/reference_markers.py \
    --counts data/interim/reference_merged.h5ad --annotation data/interim/reference_annotation.tsv.gz \
    --out results/spatial/reference_markers.tsv > logs/reference_markers.log 2>&1 < /dev/null'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/spatial/reference_markers.tsv" results/spatial/
```

Check: 20 rows per state in `results/spatial/reference_markers.tsv` (a
state with fewer passes the filters with fewer genes).

## 3.5 Mapping the sections

Visium sections, one cell2location model per section on the GPU, in a
loop that skips finished sections (`workflow/scripts/run_c2l_sections.sh`),
after the signatures are in place (the memory guard from 3.4 is still
running):

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader && who && \
  setsid -f nohup workflow/scripts/run_c2l_sections.sh data data/interim/spatial/signatures/reference_signatures.csv.gz \
    results/spatial/reference_markers.tsv shiers2025_dpn_drg utdallas_dpn_nerves tavaresferreira2022_hdrg_visium \
    zhang2024_hdrg_visium vroman2023_hdrg_visium > logs/c2l_sections.log 2>&1 < /dev/null'
```

The two targeted panels run on the CPU while the sections map. The panel
markers for check (a) are the reference markers of 3.4 computed within
the panel's genes (the gene list is taken from the first section of each
dataset); then the signature-likelihood transfer of plan amendment 2
(parameters: `panel_transfer` in the config; check (b) rotates over ten
folds of the shared genes, `heldout_gene_fraction` in the config). The
first attempt, a scANVI and scArches transfer, is described in the
analysis log; its script version is in the history of
`workflow/scripts/panel_transfer.py` and its outputs were removed.

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && rm -rf data/interim/spatial/transfer && \
  setsid -f nohup bash -c "source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && \
    for d in heming2025_sural_xenium yu2024_hdrg_xenium; do \
      python -c \"import anndata as ad, glob; a = ad.read_h5ad(sorted(glob.glob(\\\"data/interim/spatial/\$d/*.h5ad\\\"))[0], backed=\\\"r\\\"); print(chr(10).join(a.var_names))\" > results/spatial/\$d/panel_genes.txt; \
      nice -n 15 /usr/bin/time -v python workflow/scripts/reference_markers.py --counts data/interim/reference_merged.h5ad \
        --annotation data/interim/reference_annotation.tsv.gz --genes results/spatial/\$d/panel_genes.txt \
        --out results/spatial/\$d/panel_markers.tsv > logs/panel_markers_\$d.log 2>&1; echo \"MARKERS \$d \$?\"; \
      nice -n 15 python workflow/scripts/panel_transfer.py \$d --signatures data/interim/spatial/signatures/reference_signatures.csv.gz \
        --annotation data/interim/reference_annotation.tsv.gz --markers results/spatial/\$d/panel_markers.tsv \
        --sections-dir data/interim/spatial/\$d --out-dir data/interim/spatial/transfer/\$d \
        --results-dir results/spatial/\$d > logs/panel_transfer_\$d.log 2>&1; echo \"TRANSFER \$d \$?\"; \
    done; echo PANELDONE" > logs/panel_rerun.log 2>&1 < /dev/null'
```

Check: `logs/panel_rerun.log` has `MARKERS <dataset> 0`, `TRANSFER <dataset> 0`
for both datasets and ends with `PANELDONE`; `logs/panel_transfer_<dataset>.log`
reports 96 shared genes for Heming 2025 and 85 for Yu 2024, both in ten
folds; `results/spatial/<dataset>/transfer_summary.tsv` has one row per
section. The same two commands run on the laptop against `$DATA_ROOT`
give the same tables up to the last digits of the AUROC values.

When the section loops end with `ALLDONE`, collect the checks and bring
back the small outputs:

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  python workflow/scripts/collect_spatial_checks.py --c2l-dir data/interim/spatial/c2l \
    --results-dir results/spatial --sections results/spatial/section_summary.tsv'
rsync -a --include "*/" --include "*.tsv" --exclude "*" "$REMOTE:~/Projects/pnatlas-dev/results/spatial/" results/spatial/
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/interim/spatial/c2l" "$REMOTE:~/Projects/pnatlas-dev/data/interim/spatial/transfer" "$DATA_ROOT/interim/spatial/"
```

Check: 44 `EXIT <dataset> <section> 0` lines; `results/spatial/state_mapping_status.tsv`
lists every state with the sections checked and its mapped / not mapped status.

Running the sections side by side. One loop takes about 17 minutes per
section (30,000 epochs on the whole section, about 28 epochs per second,
5 GB of GPU memory) while the GPU is only two thirds busy, so after the
first 17 sections the single loop was stopped (its running section left
to finish, with a claim directory so that no loop restarts it) and two
loops were started on disjoint dataset lists; the loop script skips
finished and claimed sections. Alongside them runs the GPU yield guard
(`workflow/remote/gpu_yield.sh`), which pauses our newest mapping job
whenever another person's compute process uses the GPU or free memory
falls below 25 percent, down to one running job, and resumes the paused
jobs once the GPU has been free of others for three checks.

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && \
  cur=$(ps -eo args | grep "[c]2l_map.py" | sed -E "s/.*--out-dir ([^ ]+).*/\1/"); mkdir -p "$cur/.claim"; \
  kill $(pgrep -f "run_c2l_sections.sh data")'
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  setsid -f nohup workflow/scripts/run_c2l_sections.sh data data/interim/spatial/signatures/reference_signatures.csv.gz \
    results/spatial/reference_markers.tsv utdallas_dpn_nerves zhang2024_hdrg_visium vroman2023_hdrg_visium \
    > logs/c2l_sections_A.log 2>&1 < /dev/null; \
  setsid -f nohup workflow/scripts/run_c2l_sections.sh data data/interim/spatial/signatures/reference_signatures.csv.gz \
    results/spatial/reference_markers.tsv tavaresferreira2022_hdrg_visium > logs/c2l_sections_B.log 2>&1 < /dev/null; \
  workflow/remote/gpu_yield.sh start logs/gpu_yield.log "[c]2l_map.py"'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && workflow/remote/gpu_yield.sh status logs/gpu_yield.log'
```

Check: both loop logs end with `ALLDONE`; `logs/gpu_yield.log` lists every
PAUSE and RESUME with the process that caused it. When the loops are
done, stop both guards (the GPU guard resumes anything still paused),
confirm nothing of ours is left, then run the collection block above:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && workflow/remote/gpu_yield.sh stop logs/gpu_yield.log && \
  ~/bin/memguard.sh stop logs/memguard.log; pgrep -af "c2l|run_c2l|gpu_yield|memguard" | grep -v pgrep; \
  grep -c -E "PAUSE|RESUME" logs/gpu_yield.log; grep -E "WARN|PAUSE" logs/memguard.log | tail'
```

Check: `pgrep` prints nothing; 44 `checks.tsv` files under the c2l
directory; `results/spatial/section_checks.tsv` has one row per section
and state (1,540 rows) and `results/spatial/state_mapping_status.tsv`
28 state and compartment rows.

## 3.6 Anatomy check inputs and panels

Check (c) needs a tissue image under the state map. Images are taken
from four places: the Tavares-Ferreira H and E images on Box (16 TIFF
files), the Space Ranger high-resolution tissue image and scale factors
of every Shiers 2025 section from the SPARC copy of that dataset
(discover dataset 490, public S3 bucket), the DAPI morphology images
inside the Heming 2025 GEO tar, and the low-resolution image kept in the
DPN nerve Seurat object. List and download the first two (the download
script appends to the manifest, so nothing else may write it meanwhile;
`sync.sh` copies the laptop's manifest to the server, so it is not run
while a download is going):

```bash
cd "$REPO"
python3 -I workflow/scripts/stage3_list_files.py --images > data/metadata/stage3_images.tsv
awk -F'\t' 'NR>1{s+=$5} END{print NR-1, "files, GB", s/1073741824}' data/metadata/stage3_images.tsv
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && df -h / | tail -1 && \
  setsid -f nohup python3 -I workflow/scripts/download_files.py data/metadata/stage3_images.tsv \
    --data-root data --manifest data/MANIFEST.tsv \
    --fetch "workflow/scripts/download_files.py data/metadata/stage3_images.tsv" \
    > logs/stage3_images_download.log 2>&1 < /dev/null'
```

Meanwhile take the morphology images out of the Heming tar (gzipped
OME-TIFF, 10.5 GB; they land next to the tar and are moved into the
extracted directory once decompressed):

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev/data/raw/heming2025_sural_xenium && \
  setsid -f nohup bash -c "nice -n 15 tar -xf GSE285984_RAW.tar --wildcards \"*morphology.ome.tif.gz\" && \
    nice -n 15 gunzip -f GSM87144*_morphology.ome.tif.gz && mv GSM87144*_morphology.ome.tif extracted/ && echo MORPHDONE" \
    > ../../../logs/heming_morphology_extract.log 2>&1 < /dev/null'
```

When the download has ended (`grep -c ^done logs/stage3_images_download.log`
equals the number of rows in the list), verify, bring the manifest
back, and re-export the DPN nerve object with the image export added to
`rds_spatial_to_mtx.R` (same counts, coordinates and metadata as before,
plus `images/<section>.png`):

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && python3 -I workflow/scripts/download_files.py data/metadata/stage3_images.tsv --data-root data --verify | sort | uniq -c -w 2'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/MANIFEST.tsv" data/MANIFEST.tsv
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas-r && cd ~/Projects/pnatlas-dev && \
  nice -n 15 Rscript workflow/scripts/rds_spatial_to_mtx.R \
    data/raw/utdallas_dpn_nerves/merged_obj_harmony_layerintegration_resol0.5_annotated_final.rds \
    data/interim/rds/utdallas_dpn_nerves > logs/rds_export_utdallas_dpn_nerves.log 2>&1; \
  grep -E "^image_|spots_with" logs/rds_export_utdallas_dpn_nerves.log'
```

Then draw the panels (image sources per dataset in the `anatomy` block
of `config/spatial_inputs.yaml`: for a Visium section the image and one
map per state of the estimated cells of that state in each spot; for a
Xenium section the image and the cells of the three states in colour)
and bring them back with the reader's sheet. The Xenium image is JPEG
2000 compressed, so the environment needs `imagecodecs` next to
`tifffile` (both in `environment.yml`; on the existing server
environment `pip install imagecodecs`):

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  pip install -q imagecodecs && \
  for d in heming2025_sural_xenium shiers2025_dpn_drg utdallas_dpn_nerves; do \
    nice -n 15 python workflow/scripts/anatomy_panels.py $d --data-root data --out-dir results/spatial/anatomy/$d; done'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/spatial/anatomy/" results/spatial/anatomy/
```

Check: one JPEG per section with coordinates and an image (8 Heming, 16
Shiers, 6 DPN nerve, 26 MB in all); `results/spatial/anatomy/anatomy_check.tsv`
has three rows per section (one per state) with the count drawn and an
empty `placed_correctly` column for the reader to fill with yes, no or
not present. The Tavares-Ferreira images have no spot positions to
overlay (the positions are only in the Loupe files, which no open tool
reads), so they are kept for viewing and no panel is drawn for them.

## 3.7 Yu 2024 full Xenium bundles from Dryad

The four Xenium output bundles of Yu 2024 (all cells, the morphology
images and the post-Xenium H and E images; 12.4 GB in five files) are on
Dryad under doi 10.5061/dryad.gf1vhhmxq (version 3, CC0). Dryad serves
file downloads only to a browser session or with an API token (section
1.7), so the files are saved by hand and then verified against the
sha-256 digests that the Dryad API publishes. Manual steps, in a browser
on the laptop:

1. Open https://datadryad.org/dataset/doi:10.5061/dryad.gf1vhhmxq.
2. Under "Oct 25, 2024 version files" click each file name in turn
   (`README.md`, then `output-XETG00171__0018220__Region_1__20240207__003742.rar`
   to `Region_4`); each click opens a download tab and the browser saves
   the file to its download folder. With little free space on the laptop
   system disk, take one bundle at a time and move it on before the next.
3. Move each finished file into the repo data directory (on the external
   drive) under `data/raw/yu2024_hdrg_xenium/dryad/`, keeping its name.

Then verify and record (the script reads the file list of the dataset
version from the Dryad API, compares size and sha-256, computes the md5
and appends one manifest row per verified file; rerun it as the files
come in):

```bash
cd "$REPO"
mkdir -p data/raw/yu2024_hdrg_xenium/dryad
mv ~/Downloads/README.md ~/Downloads/output-XETG00171__0018220__Region_*.rar data/raw/yu2024_hdrg_xenium/dryad/
python3 -I workflow/scripts/dryad_verify.py --list
python3 -I workflow/scripts/dryad_verify.py --data-root data --dir raw/yu2024_hdrg_xenium/dryad \
  --dataset yu2024_hdrg_xenium --manifest data/MANIFEST.tsv
bsdtar -tf data/raw/yu2024_hdrg_xenium/dryad/output-XETG00171__0018220__Region_1__20240207__003742.rar | head -60
```

Check: five `ok` lines, `data/MANIFEST.tsv` gains five rows with
`source` Dryad, and the listing shows a Xenium output folder
(`experiment.xenium`, `cells.csv.gz`, `cell_feature_matrix.h5`,
`transcripts.parquet`, the `morphology*.ome.tif` images and the
`DRG_area_<n>` H and E image with its alignment file).

Copy the bundles to the server (the laptop's manifest is now the one
with the Dryad rows, so `sync.sh` may run again) and keep the laptop
copy on the external drive until the server copy has verified:

```bash
rsync -a --partial --info=progress2 data/raw/yu2024_hdrg_xenium/dryad/ "$REMOTE:~/Projects/pnatlas-dev/data/raw/yu2024_hdrg_xenium/dryad/"
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && python3 -I workflow/scripts/dryad_verify.py --data-root data --dir raw/yu2024_hdrg_xenium/dryad --dataset yu2024_hdrg_xenium --manifest /dev/null'
```

Check: five `ok` lines on the server too.

## 3.8 Project data on the server's data array

The server's system disk (907 GB, shared by all home directories) was
91 percent full with the project's `data/raw` (72 GB) and `data/interim`
(13 GB) in the checkout, and the Dryad bundles unpack to a further
15 GB or so. The two large directories are moved to the server's data
array and replaced by symbolic links in the checkout, so `--data-root
data` and every command above are unchanged. The target directory is
passed at run time as `REMOTE_DATA` (a directory of the array owned by
the login user); the script copies each directory, compares the copy
with the source by checksum, and only then removes the source and
links it, writing a README next to the target. Run with nothing of
the project running on the server:

```bash
cd "$REPO"
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" "cd ~/Projects/pnatlas-dev && setsid -f nohup nice -n 15 workflow/remote/relocate_data.sh $REMOTE_DATA raw interim > logs/relocate_data.log 2>&1 < /dev/null"
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/relocate_data.log && ls -l data/ && python3 -I workflow/scripts/dryad_verify.py --data-root data --dir raw/yu2024_hdrg_xenium/dryad --dataset yu2024_hdrg_xenium --manifest /dev/null'
```

Check: the log ends with `done: data/raw -> ...` and `done: data/interim
-> ...` and no `DIFFERENCES` line; `ls -l data/` shows the two links;
the Dryad files verify through the link. Before this step the Dryad
bundles were copied into the checkout with `rsync -a --partial` from
the laptop (section 3.7), resumed after an interrupted connection, and
verified on the server with the same script.

## 3.9 Yu 2024: every cell of the Dryad bundles, and H and E panels

The Dryad bundles (3.7) hold, per region, the cells table and the cell by
gene matrix of every segmented cell (27,000 to 46,000 cells per region,
against the 300 to 440 hand-segmented neurons of the GEO table), the
post-Xenium H and E image of the section and its alignment to the
morphology image from Xenium Explorer. The nine members the analysis
needs (0.1 to 0.3 GB per region) are taken out of each rar with libarchive
(`workflow/scripts/extract_members.py`; `python-libarchive-c` in
`environment.yml`, on the existing server environment `pip install
libarchive-c`) into `data/raw/yu2024_hdrg_xenium/dryad/extracted/Region_<n>/`.
The all-cell data are a second dataset, `yu2024_hdrg_xenium_cells` in
`config/spatial_inputs.yaml` (format `xenium`, same donors, slide and
levels as the neuron table), transferred with the signature classifier of
3.5 with the ganglion states plus the three nerve sheath states as
candidates (`extra_states` in `config/config.yaml`), and drawn for check
(c) on the H and E through the inverse of the alignment matrix (`anatomy`
block of kind `he_aligned`). The neuron table's coordinates are corrected
at the same time (`position_order: yx`, `pixel_size_um` in the config):
its cell ids are ImageJ ROI names, row-column of the morphology image.

All steps run in one script on the server, after the usual checks and
with the memory guard running:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  pip install -q libarchive-c && nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader && who && free -g | sed -n 2p && \
  ~/bin/memguard.sh start logs/memguard.log "extract_members|build_spatial|reference_markers|panel_transfer|anatomy_panels" && \
  setsid -f nohup workflow/scripts/run_yu_dryad.sh data > logs/yu_dryad.log 2>&1 < /dev/null'
```

Check: `logs/yu_dryad.log` has four `EXTRACT Region_<n> 0` lines, two
`BUILD <dataset> 0` lines, `MARKERS`, `TRANSFER`, `COLLECT` and `PANELS`
with exit 0, and ends with `YUDONE`; each `logs/extract_dryad_Region_<n>.log`
lists nine members; `results/spatial/section_summary.tsv` has the rows
of the two Yu datasets replaced and every other dataset's rows kept (60
sections). Then stop the guard and bring the
tables, the panels and the per-cell labels back:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/yu_dryad.log && ~/bin/memguard.sh stop logs/memguard.log; grep -E "WARN|PAUSE" logs/memguard.log | tail -3'
rsync -a --include "*/" --include "*.tsv" --include "*.txt" --exclude "*" "$REMOTE:~/Projects/pnatlas-dev/results/spatial/" results/spatial/
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/spatial/anatomy/" results/spatial/anatomy/
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/interim/spatial/transfer/yu2024_hdrg_xenium_cells" "$DATA_ROOT/interim/spatial/transfer/"
```

Check: `results/spatial/yu2024_hdrg_xenium_cells/transfer_summary.tsv` has
four rows; `results/spatial/anatomy/yu2024_hdrg_xenium_cells/` has four
JPEGs with the H and E under the cells, and `anatomy_check.tsv` twelve
new rows for the reader. The same script runs on the laptop against a
data root that holds the bundles, the reference and the signatures.

## 3.10 Yu 2024: the sympathetic neuron state is not a DRG candidate

A DRG holds no sympathetic neurons, and in 3.9 the sympathetic neuron
state of the ganglion compartment took 15,593 of the all-cell dataset's
cells. The state is removed from the candidates of both Yu datasets
(`exclude_states` under `panel_transfer.datasets` in `config/config.yaml`,
read by `workflow/scripts/panel_transfer.py`), the two transfers are
rerun on the section objects, panel genes and panel markers of 3.9, the
checks are collected again and the H and E panels of the all-cell
dataset redrawn, all in one script on the server:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader && who && free -g | sed -n 2p && \
  ~/bin/memguard.sh start logs/memguard.log "panel_transfer|anatomy_panels" && \
  setsid -f nohup workflow/scripts/run_yu_transfer.sh data > logs/yu_transfer.log 2>&1 < /dev/null'
```

Check: `logs/yu_transfer.log` has two `TRANSFER <dataset> 0` lines,
`COLLECT 0`, `PANELS 0` and ends with `YUTRANSFERDONE`; the transfer logs
report 23 candidate states for the all-cell dataset and 11 for the
neurons. Then stop the guard and bring the tables, panels and per-cell
labels back as in 3.9:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/yu_transfer.log && ~/bin/memguard.sh stop logs/memguard.log; grep -E "WARN|PAUSE" logs/memguard.log | tail -3'
rsync -a --include "*/" --include "*.tsv" --include "*.txt" --exclude "*" "$REMOTE:~/Projects/pnatlas-dev/results/spatial/" results/spatial/
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/spatial/anatomy/yu2024_hdrg_xenium_cells/" results/spatial/anatomy/yu2024_hdrg_xenium_cells/
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/interim/spatial/transfer/yu2024_hdrg_xenium" "$REMOTE:~/Projects/pnatlas-dev/data/interim/spatial/transfer/yu2024_hdrg_xenium_cells" "$DATA_ROOT/interim/spatial/transfer/"
```

Check: `results/spatial/state_mapping_status.tsv` has no sympathetic
neuron row in the ganglion, and the per-cell files of the transfer
directories carry no sympathetic neuron label.

## 3.11 Anatomy check (c): the panels read and the sheet filled

Two things had to change before the panels could be read. The DPN nerve
object is Seurat 5, whose `GetTissueCoordinates` returns `x` as the image
row and `y` as the image column, so the spots of the six UT Dallas
sections fell beside the tissue; `position_order: yx` in the dataset's
config block makes `build_spatial.py` read the exported `x,y` as
row,column (the eight orientations of the spots were scored against the
tissue mask of the low-resolution image and only this one puts 92 to 99
percent of the spots on tissue). For a targeted panel section the single
cell panel, on which a few hundred coloured cells among tens of thousands
of grey ones could not be seen, is replaced by one panel per sheath state
with that state's cells drawn at twice the marker width over the other
cells in grey, as for the Visium sections. Rebuild the UT Dallas sections
and redraw the panels of the three datasets concerned:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  d=utdallas_dpn_nerves && nice -n 15 python workflow/scripts/build_spatial.py $d --data-root data \
    --out-dir data/interim/spatial/$d > results/spatial/sections_$d.tsv 2> logs/build_spatial_$d.log && \
  for d in heming2025_sural_xenium yu2024_hdrg_xenium_cells utdallas_dpn_nerves; do \
    nice -n 15 python workflow/scripts/anatomy_panels.py $d --data-root data --out-dir results/spatial/anatomy/$d --width-px 1100; done'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/spatial/anatomy/" results/spatial/anatomy/
```

Reading. One reader goes through the 34 sections (8 Heming, 16 Shiers,
6 DPN nerve, 4 Yu) and writes yes, no or not present with a short note
into the `placed_correctly` and `note` columns of
`results/spatial/anatomy/anatomy_check.tsv`, one row per section and
state, against the image in the first panel: perineurial cells or spots
must lie along the fascicle borders (nerve) or the capsule and root
sheaths (ganglion); epineurial fibroblasts in the outer and
interfascicular epineurium, the capsule or the loose connective and
adipose tissue outside it; endoneurial fibroblasts inside the fascicles
or the parenchyma and not at the margins. "Yes" needs the highest spots
or most of the cells in the expected place, "no" is a state spread over
the section regardless of its structure or concentrated in the wrong
place, "not present" a state with too few cells to judge (fewer than
about 20). The reading is done on the panels as drawn; for the Yu
sections, whose sheath cells are few, the same cells were also viewed
on the H and E at pyramid level 1 with larger markers. Then collect the
checks again (the collector now reads the sheet and reports the sections
read and those read as misplaced per state and compartment, without
changing the mapped verdict):

```bash
python workflow/scripts/collect_spatial_checks.py --c2l-dir "$DATA_ROOT/interim/spatial/c2l" \
  --results-dir results/spatial --sections results/spatial/section_summary.tsv
```

Check: the sheet has 102 filled rows; `state_mapping_status.tsv` has the
`read_c` and `misplaced_c` columns filled for the three sheath states in
both compartments.
