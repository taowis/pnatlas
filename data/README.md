# Data

Nothing under `data/` is committed except this file, `MANIFEST.tsv` and
`data/metadata/`. Large files live under `data_root` from `config/paths.yaml`
(on the server, `data/` inside the project directory). `MANIFEST.tsv` has one
row per downloaded file: path relative to `data_root`, source, exact URL,
access date, md5, size and the fetch command. Re-download everything with

```bash
python3 -I workflow/scripts/download_files.py data/metadata/stage2_files.tsv --data-root "$DATA_ROOT"
python3 -I workflow/scripts/download_files.py data/metadata/stage2_files.tsv --data-root "$DATA_ROOT" --verify
```

Candidate datasets and their access terms are in `docs/datasets.md`; the
entries used by the analysis are in `config/datasets.yaml`.

## Stage 2: reference layer (single-cell)

### bhuiyan2025_hdrg

DRG, snRNA-seq, scRNA-seq. DRG backbone of the reference (126 donors, nomenclature source for the neuron and non-neuron labels). Source: humandrg.painseq.com (Seurat object on Dropbox), PMID 41279342. Access: portal.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/bhuiyan2025_hdrg/seurat_object_preprint_release.Rds` | 1915388302 | b45ecd7a86901febf7b4d331d54cb3fb | 2026-11-07 |

Download link: the "Download Seurat object" button on humandrg.painseq.com (a public Dropbox link, recorded in `data/metadata/stage2_files.tsv`). No licence is stated on the portal; the object is used with citation and nothing cell-level derived from it is redistributed.

### nguyen2021_hdrg

DRG, snRNA-seq. DRG neurons, snRNA-seq, 6 donors. Source: GSE168243, PMID 34825887. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/nguyen2021_hdrg/GSE168243_RAW.tar` | 9349120 | 8cc74d3afbba048cf4f57c187948fa2a | 2026-11-07 |

### jung2023_hdrg

DRG, snRNA-seq. DRG snRNA-seq, 7 donors; the human sub-series of the cross-species series GSE201654, with the authors' annotated Seurat object. Source: GSE201586, PMID 36690629. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/jung2023_hdrg/GSE201586_RAW.tar` | 2115338240 | b55b572095de66a7028631cbe217fc94 | 2026-11-07 |
| `raw/jung2023_hdrg/GSE201586_humanallcelltypes.rds.gz` | 4223086818 | a1074dc8312d8ffe5519dd83f313da60 | 2026-11-07 |

### yu2024_singlesoma

DRG, single-soma Smart-seq2. DRG single-soma Smart-seq2 counts; neuron subtype anchor. Source: GSE249746, PMID 39496796. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/yu2024_singlesoma/GSE249746_Expression_matrix_raw_counts.csv.gz` | 11830344 | 9f230acb1c49bbc30ffeb916f15b9df3 | 2026-11-07 |

### zhang_hdrg_glia

DRG, spinal cord, snRNA-seq. DRG glia, 3 donors; optional, no linked paper, cited by accession. Source: GSE189501. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/zhang_hdrg_glia/GSE189501_DRG_10X_raw_UMI_count_matrix.txt.gz` | 73029906 | 7a6a58868f1394aa3cbbb62254487beb | 2026-11-07 |

### yang2022_htg

trigeminal ganglion, snRNA-seq, snATAC-seq. Trigeminal ganglion snRNA-seq (human files only; the snATAC and mouse files are not used). Source: GSE197289, PMID 35349784. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/yang2022_htg/GSE197289_snRNA-seq_human_barcode_meta.csv.gz` | 842041 | 9ac5f3037e8d6406720c688f8ab03671 | 2026-11-07 |
| `raw/yang2022_htg/GSE197289_snRNA-seq_human_raw_counts.RDS.gz` | 182486123 | 3b534f263fc38e6deafc06e8721f9725 | 2026-11-07 |

### yang2025_hsg_multiome

sympathetic ganglia, snMultiome (RNA + ATAC). Sympathetic ganglia snMultiome, RNA part only (filtered Cell Ranger ARC matrix, author metadata, SDS sheets); the ATAC fragments (17.8 GB) are not downloaded. Source: SPARC 465, doi 10.26275/oqxx-si9h, PMID 41000833. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/yang2025_hsg_multiome/README.txt` | 808 | ce9236c4f6db98d27c036850b4f79adf | 2026-11-07 |
| `raw/yang2025_hsg_multiome/barcodes.tsv.gz` | 422581 | a7887feff1b59fb46c34681fc92ce0fa | 2026-11-07 |
| `raw/yang2025_hsg_multiome/dataset_description.xlsx` | 15236 | b3558ac18218df4307f1587b2ddc621b | 2026-11-07 |
| `raw/yang2025_hsg_multiome/features.tsv.gz` | 3400616 | 9e5c16b305b4a6150afb93a0fb81543a | 2026-11-07 |
| `raw/yang2025_hsg_multiome/matrix.mtx.gz` | 2168008189 | 67d1a4fd76ff1c6a1db0ccefcf4c5331 | 2026-11-07 |
| `raw/yang2025_hsg_multiome/meta.data.csv` | 11996682 | 201de11a44ceb49402ee55946326f351 | 2026-11-07 |
| `raw/yang2025_hsg_multiome/samples.xlsx` | 9795 | 34e3083898ceac6b9345246058d75d22 | 2026-11-07 |
| `raw/yang2025_hsg_multiome/subjects.xlsx` | 15962 | a4c36daf398537c70ef92894c8bf485b | 2026-11-07 |

Files are read from the public SPARC S3 bucket (`sparc-prod-aod-discover-publish50-use1`, dataset 465 version 1); CC BY 4.0.

### heming2025_sural_snrna

sural nerve, snRNA-seq. Sural nerve snRNA-seq, 33 polyneuropathy and 4 control donors; nerve backbone of the reference. Source: GSE285983, PMID 40849297. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/heming2025_sural_snrna/GSE285983_RAW.tar` | 1429237760 | 0b7aaf8f059425efe2d1ae4ca46a4049 | 2026-11-07 |
| `raw/heming2025_sural_snrna/GSE285983_metadata_all.csv.gz` | 22228346 | b740a4ff1dc27691ff84cbedb120e8d1 | 2026-11-07 |
| `raw/heming2025_sural_snrna/GSE285983_metadata_ic.csv.gz` | 810136 | df2f4fb84ffbdf6a976864207bbcb33b | 2026-11-07 |
