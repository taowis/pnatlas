# Stage 1: data inventory and access

Variables: `$REPO` (repo root), set as in [README.md](README.md).
Results are summarised in `docs/datasets.md`; what was found at each step is in
`docs/analysis_log/01_datasets.md`. All steps below were run at commit 89ba816.

### 1.1 Candidate list

The candidate datasets come from the literature and database search in
`docs/novelty_check.md`, section 4. That search is described there, not repeated
here.

### 1.2 Verify candidate GEO series

Env: system python3 (stdlib only).

```bash
cd "$REPO"
python3 -I workflow/scripts/stage1_verify_geo.py \
  GSE168243 GSE201654 GSE201586 GSE249746 GSE189501 GSE169301 GSE310000 \
  GSE296776 GSE197289 GSE285983 GSE285984 GSE273557 GSE295206 GSE243073 \
  GSE215994 GSE78150 GSE78148 GSE77968 GSE148059 GSE24290 GSE255436 \
  GSE165569 GSE59739 GSE139088 GSE154659 GSE175421 GSE182098 GSE142541 \
  > data/metadata/geo_verification.tsv
wc -l data/metadata/geo_verification.tsv
cut -f2 data/metadata/geo_verification.tsv | sort | uniq -c
```

Output: `data/metadata/geo_verification.tsv`.
Check: 29 lines (header plus 28), 28 `found`.
Notes: GSE189501 and GSE165569 have an empty `pubmed_ids` column. GEO records can
be edited by submitters, so a rerun may differ in small ways (one rerun listed
the GSE142541 taxa in the other order).

### 1.3 Check access terms of the non-GEO sources

Env: system python3.

```bash
cd "$REPO"
python3 -I workflow/scripts/stage1_verify_nongeo.py > data/metadata/nongeo_access.tsv
wc -l data/metadata/nongeo_access.tsv
```

Output: `data/metadata/nongeo_access.tsv`, 5 rows (2 SPARC, 1 Dryad, 2 painseq
Dropbox links). Check: `wc -l` gives 6.
What it queries: DataCite (`https://api.datacite.org/dois/<doi>`), Pennsieve
(`https://api.pennsieve.io/discover/datasets/<id>`), Dryad
(`https://datadryad.org/api/v2/datasets/<doi>` and its `/files`), and a HEAD
request on each Dropbox link for the size.
Key values: SPARC 465 CC BY 4.0, public; SPARC 464 embargoed, release
2027-04-29; Dryad CC0, 12,412,536,831 bytes; painseq Seurat object
1,915,388,302 bytes; painseq harmonized folder 19,935,643,591 bytes.

### 1.4 Look for a data use agreement on the embargoed SPARC dataset

```bash
curl -s -o /dev/null -w "%{http_code}\n" \
  "https://api.pennsieve.io/discover/datasets/464/data-use-agreement"
```

Output: `204` (no agreement attached).

### 1.5 Count donors and list files in SPARC 465 (sympathetic ganglia)

Downloads one 16 KB file.

```bash
cd "$REPO"
WORK=$(mktemp -d)
curl -s -o "$WORK/subjects.xlsx" \
  "https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com/465/files/subjects.xlsx"
shasum -a 256 "$WORK/subjects.xlsx"
python3 -I workflow/scripts/sds_subjects_summary.py "$WORK/subjects.xlsx"
curl -s "https://api.pennsieve.io/discover/datasets/465/versions/1/files?path=files/subjects.xlsx" \
  | python3 -c 'import sys,json,base64; print(base64.b64decode(json.load(sys.stdin)["sha256"]).hex())'
B="https://api.pennsieve.io/discover/datasets/465/versions/1/files/browse"
for p in "" files files/primary files/primary/filtered_feature_bc_matrix; do
  curl -s "$B?path=$p&limit=100" | python3 -c '
import sys, json
for f in json.load(sys.stdin).get("files", []):
    print(f.get("type"), f.get("path"), f.get("size"))'
done
```

Check: both sha256 values are
`649b8a260f34b341477dd8b84519433448b378098f866bb6f67f217cddf5016d`.
Output: 7 donors, 5 male, 2 female, 23 to 59 years. Processed data are
`files/primary/filtered_feature_bc_matrix/` (2,171,831,386 bytes) and
`files/primary/meta.data.csv`; `atac_fragments.tsv.gz` is 17,753,218,897 bytes.

### 1.6 Read the painseq portals and the preprint licence

Manual, in a browser:

1. Open https://humandrg.painseq.com (redirects to
   https://painseq.shinyapps.io/u19humandrgatlas/). The page links the preprint
   (10.1101/2025.11.05.686654) and the Seurat object
   `seurat_object_preprint_release.Rds` on Dropbox. No licence or terms of use
   are shown.
2. Open https://harmonized.painseq.com (redirects to
   https://painseq.shinyapps.io/harmonized_painseq_v1/), "Data download" tab. It
   links a Dropbox folder with the full atlas objects and Table S4, and GEO
   GSE255436. No licence or terms of use are shown.

Licence of the preprint (covers the article only):

```bash
curl -s "https://api.biorxiv.org/details/biorxiv/10.1101/2025.11.05.686654" \
  | python3 -c 'import sys, json; [print(c["version"], c["license"]) for c in json.load(sys.stdin)["collection"]]'
curl -s "https://api.biorxiv.org/pubs/biorxiv/10.1101/2025.11.05.686654" \
  | python3 -c 'import sys, json; d=json.load(sys.stdin); print(d["messages"][0]["status"], len(d["collection"]))'
```

Output: version 1, `cc_by_nc_nd`; no published version.

### 1.7 Dryad file access

```bash
F=$(curl -s "https://datadryad.org/api/v2/versions/324492/files" | python3 -c '
import sys, json
for f in json.load(sys.stdin)["_embedded"]["stash:files"]:
    if f["path"] == "README.md": print(f["_links"]["stash:download"]["href"])')
curl -s -o /dev/null -w "api_download %{http_code}\n" "https://datadryad.org$F"
curl -sL -o /dev/null -w "web_download %{http_code}\n" "https://datadryad.org/downloads/file_stream/${F//[^0-9]/}"
```

Output: API download 401 (needs a Dryad API token), scripted web download 403.
The files have to be downloaded from
https://datadryad.org/dataset/doi:10.5061/dryad.gf1vhhmxq in a browser, or with
a Dryad API token. Not downloaded yet.

### 1.8 Zenodo record for Heming 2025

```bash
cd "$REPO"
python3 -I workflow/scripts/stage1_lookup_sources.py zenodo 15108216
```

Output: title "PNS Atlas", concept DOI 10.5281/zenodo.14226218, CC BY 4.0, open.
36 `.qs` files, all figure data objects (largest `perineurial_figure.qs`,
312,674,344 bytes); no merged count matrix or Seurat object. The concept DOI
redirects to the latest version, 15108216.

### 1.9 sensoryomics.com

The site is a single-page application; `curl` returns only the shell. The
list of applications was read from its bundled script:

```bash
curl -sL https://sensoryomics.com/static/js/main.b40bede6.js \
  | grep -o -E 'https?://[^"'"'"' )]{10,120}' | grep -i shinyapps | sort -u
```

Output: 11 shinyapps.io applications, including `human-peripheral-nerves`,
`Human-SCG`, `SC_Visium_Dashboard`, `human-surgical-c2drg`, `RNA-Data`,
`RNAIsoforms`, `Proteomics`. The applications themselves need a browser
(a scripted request returns a "Please Wait" page). Not yet opened.

Data availability statement of the DPN nerve paper (PMID 40828619):

```bash
python3 -I workflow/scripts/stage1_lookup_sources.py epmc-fulltext PMC12520680 'Data and code availability'
```

Output: raw RNA-seq in dbGaP phs001158; processed data at
`sensoryomics.shinyapps.io/human-peripheral-nerves`; code on the utdal GitHub.

### 1.10 Source papers of GSE189501 and GSE165569

```bash
cd "$REPO"
for g in GSE189501 GSE165569; do
  python3 -I workflow/scripts/stage1_lookup_sources.py geo $g
  python3 -I workflow/scripts/stage1_lookup_sources.py epmc-accession $g
done
python3 -I workflow/scripts/stage1_lookup_sources.py pubmed 33686078
python3 -I workflow/scripts/stage1_lookup_sources.py epmc-fulltext PMC10945563 'GSE[0-9]+'
```

Output: GSE165569 has the same title as Kupari 2021 Nat Commun (PMID 33686078)
and Europe PMC lists that paper for the accession. GSE189501 (Zhang D, West
China Hospital; 3 donors; DRG and spinal cord glia; raw UMI matrices as
supplementary files) is cited by one 2026 PLoS One paper (PMID 42546053) but no
source paper was found; the same group's eLife spinal cord atlas (PMID 38289829,
PMC10945563) lists GSE103892, GSE154659 and GSE243077, not GSE189501.

### 1.11 PMID of the sympathetic chain preprint

```bash
python3 -I workflow/scripts/stage1_lookup_sources.py pubmed 40667173
```

Output: PubMed's record of the bioRxiv preprint 10.1101/2025.06.24.661373 (no
journal version). Matches the `Human-SCG` application on sensoryomics.

### 1.12 Fill config/datasets.yaml

Written by hand from the results above; 22 entries. Check:

```bash
grep -c '^  - id:' config/datasets.yaml
```

### 1.13 sensoryomics applications, read in a browser

Manual. Each application is a Shiny app on shinyapps.io and takes 10 to 20
seconds to load.

1. https://sensoryomics.shinyapps.io/human-peripheral-nerves/ (DPN nerves,
   PMID 40828619). Home tab: buttons "Explore bulk rna-seq data", "Download
   spatial data", "Explore microRNA data". Gene Search tab: per-gene qnTPM plots
   for dataset 1 (6 control vs 6 DPN sural nerves, 100M reads per sample) and
   dataset 2 (5 paired sural and tibial DPN nerves); no table download.
   "Download spatial data" opens the public Box folder
   https://utdallas.app.box.com/s/apafnhru6dihoi2va20wxexi40lia9b0
   ("Visium_DPN_human_sural_nerve", 7 files, 3.6 GB):
   `merged_obj_harmony_layerintegration_resol0.5_annotated_final.rds` (675.3 MB),
   `sample_metadata.csv` (797 B), five Loupe files `cloupe8`, `cloupe11`,
   `cloupe14`, `cloupe24L`, `cloupe24T` (158.8 MB to 928.3 MB).
2. https://sensoryomics.shinyapps.io/Human-SCG/ (sympathetic chain ganglia
   Visium, PMID 40667173). Home, Gene Search and About tabs. Gene Search shows
   a CPM table and plots per gene; no download link anywhere.
3. https://sensoryomics.shinyapps.io/SC_Visium_Dashboard/ is the mouse spinal
   cord atlas of Loo 2026 (Nat Commun), not the human DRG Visium.
4. https://sensoryomics.shinyapps.io/RNA-Data/ is the human DRG app: Home tab
   links PMID 35171654 (Tavares-Ferreira 2022 DRG Visium) with "Download Loupe
   Browser Files" and "Download Data", and the Ray 2023 Brain paper
   (10.1093/brain/awac266) with a "Neuropathic Pain-TPM" tab. "Download Data"
   opens the public Box folder
   https://utdallas.app.box.com/s/dgs8a7vo9uu9cx6l81somaxb468hpgnz
   ("Visium_hDRG_shared_data_PriceLab", 83 files, 34.7 GB): `H5_files/` (17),
   `Loupe_browser_files/` (17), `H&E_images/` (16),
   `SUPPL_files_scitranslmed.abj8186_files_s1_to_s31/` (31),
   `README_sampleIDs.xlsx` (12.2 KB), `drg.combined.rds` (186.5 MB).

Nothing was downloaded. Box folder listings can be read without an account.

### 1.14 Verify genetics layer sources

Env: system python3 (stdlib only).

```bash
cd "$REPO"
python3 -I workflow/scripts/stage1_verify_genetics.py > data/metadata/genetics_sources.tsv
wc -l data/metadata/genetics_sources.tsv
grep -c '^  - id:' config/datasets.yaml
```

Output: `data/metadata/genetics_sources.tsv`, 38 rows (10 GWAS Catalog studies,
5 eQTL Catalogue datasets, 1 GTEx portal query, 16 GTEx bucket files, 6 PanelApp
panels). Check: `wc -l` gives 39; `config/datasets.yaml` has 35 entries after
the 13 genetics entries were added by hand.
What it queries: GWAS Catalog REST API v2 (`/studies/<GCST>`), the eQTL
Catalogue dataset table on GitHub (the API v2 `/datasets` endpoint returns
HTTP 410), GTEx portal API v2 (`/dataset/sample`), the `gs://adult-gtex` bucket
listing, PanelApp GE and AU (`/api/v1/panels/<id>/`; needs a User-Agent header).
PNPrank truth set: read from `pnprank-dev/data/truth/neuropathy_genes.tsv`, 242
genes.

### 1.15 Freeze the analysis plan

Manual. `docs/analysis_plan.md` rewritten from `config/datasets.yaml`,
`docs/datasets.md` and the stage 1 findings; the seven decisions are recorded
in the Stage 1 section of `docs/decisions.md`. Checks:

```bash
cd "$REPO"
grep -rn -P "\x{2014}|\x{2013}" docs config workflow
grep -rn -E "/Users/|/home/|iz89" docs config workflow
```

Both print nothing.
