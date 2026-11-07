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

and the same with `data/metadata/stage3_files.tsv` for the spatial layer.

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

## Stage 3: spatial layer

Downloaded on the server with

```bash
python3 -I workflow/scripts/download_files.py data/metadata/stage3_files.tsv --data-root "$DATA_ROOT"
python3 -I workflow/scripts/download_files.py data/metadata/stage3_files.tsv --data-root "$DATA_ROOT" --verify
```

The Box files are fetched through the shared-link download endpoint of
the UT Dallas Box account (exact URL per file in `MANIFEST.tsv`); the
sample sheets are in `data/metadata/geo_samples_spatial.tsv` for the GEO
series and in the downloaded xlsx and csv for the Box folders.

### heming2025_sural_xenium

sural nerve, Xenium. Eight sural nerve sections (two CIDP, two CIAP, two vasculitic neuropathy, two controls), Xenium outputs of all sections in one RAW tar; main nerve spatial input. Source: GSE285984, PMID 40849297. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/heming2025_sural_xenium/GSE285984_RAW.tar` | 10614650880 | 9c32428ebcc9579fe56a9ed19f634882 | 2026-11-08 |

### yu2024_hdrg_xenium

DRG, Xenium (100-gene panel). Per-region transcript tables (RAW tar) and the matrix of 1,415 hand-segmented neurons by 87 genes with positions and donor metadata. The full Xenium bundles are on Dryad (doi 10.5061/dryad.gf1vhhmxq, browser download) and are not downloaded. Source: GSE273557, PMID 39496796. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/yu2024_hdrg_xenium/GSE273557_RAW.tar` | 862965760 | a116c8fc0aa9ec742100890ea3ec7532 | 2026-11-08 |
| `raw/yu2024_hdrg_xenium/GSE273557_hDRG_gene_expression_matrix.csv.gz` | 117846 | 6bc079ebb396f3c637c2e7ec0749ed98 | 2026-11-08 |
| `raw/yu2024_hdrg_xenium/GSE273557_hDRG_gene_expression_matrix_metadata.csv.gz` | 13370 | 0a1929fcfd51ea8c0330943a69f0ebc5 | 2026-11-08 |

### zhang2024_hdrg_visium

DRG, Visium. Six donors (DRG 2 to 7), one filtered matrix archive and one image per section. Source: GSE243073 (sub-series of GSE243077), PMID 38289829. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/zhang2024_hdrg_visium/GSE243073_RAW.tar` | 291328000 | 6507239167f930143d0b30b3705ccc77 | 2026-11-08 |

### vroman2023_hdrg_visium

DRG, Visium. Two sections of one donor, Space Ranger H5, Loupe file and high-resolution image per section. Source: GSE215994, PMID 37664243. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/vroman2023_hdrg_visium/GSE215994_RAW.tar` | 809451520 | 2f16ae7ed6f9d728f47c51103f86c6d1 | 2026-11-08 |

### shiers2025_dpn_drg

DRG, Visium, donors with and without diabetes. 16 sections from seven donors in three batches: filtered and raw H5, tissue positions, Loupe files and the authors' neuron selections per section; donor diabetes status in the docx. Source: GSE295206, PMID 40325011. Access: open.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/shiers2025_dpn_drg/GSE295206_RAW.tar` | 9698846720 | 121636444c268a13bcf0d8f5bbe877bd | 2026-11-08 |
| `raw/shiers2025_dpn_drg/GSE295206_metadata_nageotteSamples.docx.gz` | 19532 | c6c9c82ab58f4c8058301cd048c49021 | 2026-11-08 |

### tavaresferreira2022_hdrg_visium

DRG, Visium. 16 Space Ranger H5 files (two sections each of four female and four male donors), the sample sheet and the combined Seurat object; the H and E images (26.6 GB) and Loupe files (10.4 GB) in the same folder are not downloaded. Source: UT Dallas Box folder linked from sensoryomics.shinyapps.io/RNA-Data, PMID 35171654. Access: open, no licence stated.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/tavaresferreira2022_hdrg_visium/README_sampleIDs.xlsx` | 12511 | 6b32ae96f53c5f3be182642db7f95afa | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF1-1_filtered_feature_bc_matrix.h5` | 7856990 | 98cac2790c37ee303f075c094a23c42b | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF1-2_filtered_feature_bc_matrix.h5` | 9560354 | 2270b223d98b0b5977c9fe8b4cc61bd8 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF2-1_filtered_feature_bc_matrix.h5` | 8954905 | 9a540e32e86d3509313a6ac4ef2c72ce | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF2-2_filtered_feature_bc_matrix.h5` | 10579454 | 7870d39e816245e126dee548acb1b974 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF3-1_filtered_feature_bc_matrix.h5` | 8888168 | 1fe14e4074c6126590908a32313e645a | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF3-2_filtered_feature_bc_matrix.h5` | 8826816 | 1eb4770c983141d07a67cf04b3cf80d8 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF4-1_filtered_feature_bc_matrix.h5` | 6280805 | 9e4ce496ab9190691eab479a7fc272d0 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VF4-2_filtered_feature_bc_matrix.h5` | 5866551 | 9805068b9e24dd930246c3d35df62da9 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM1-1_filtered_feature_bc_matrix.h5` | 7086091 | e9dea4a1877b3ab4b85311ea2c98a344 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM1-2_filtered_feature_bc_matrix.h5` | 4848035 | 6f3a2138e80416081acec41d77805ff7 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM2-1_filtered_feature_bc_matrix.h5` | 4840565 | da00826879e19c46fe12acd29e1471be | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM2-2_filtered_feature_bc_matrix.h5` | 9698557 | 5cd9cb26f9c1a90e0f46970eeb1734f4 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM3-1_filtered_feature_bc_matrix.h5` | 8526709 | 6cc890d5b01600f81c2f71dcc61941b8 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM3-2_filtered_feature_bc_matrix.h5` | 8092430 | d94584875e86ea350da1f765c848e85d | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM4-1_filtered_feature_bc_matrix.h5` | 7681706 | 69786d2eb2e9757242443751038ac416 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/VM4-2_filtered_feature_bc_matrix.h5` | 6502306 | 96413f935f76e86e84b606c8e3c95614 | 2026-11-08 |
| `raw/tavaresferreira2022_hdrg_visium/drg.combined.rds` | 195534826 | 6e40a62c807bd82f46fde8435452df33 | 2026-11-08 |

### utdallas_dpn_nerves

sural nerve, Visium. The annotated, Harmony-integrated Seurat object of five sural nerve sections (four diabetic amputation cases, one control) and the sample sheet; the Loupe files are not downloaded. Source: UT Dallas Box folder linked from sensoryomics.shinyapps.io/human-peripheral-nerves, PMID 40828619. Access: open, no licence stated.

| File (under data_root) | Size (bytes) | md5 | Access date |
|---|---|---|---|
| `raw/utdallas_dpn_nerves/merged_obj_harmony_layerintegration_resol0.5_annotated_final.rds` | 708154858 | def955735ff1b683320cb412a85bf8c1 | 2026-11-08 |
| `raw/utdallas_dpn_nerves/sample_metadata.csv` | 797 | 7364cce3886747ac84132aeff5613ff5 | 2026-11-08 |
