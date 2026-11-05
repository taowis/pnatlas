# Stage 1: data inventory and access

Variables: `$REPO` (repo root), set as in [README.md](README.md).
Results of these steps are summarised in `docs/datasets.md`.

### 1.1 Candidate list

The candidate datasets come from the literature and database search in
`docs/novelty_check.md`, section 4 (searches run 2026-11-02). That search is
described there, not repeated here.

### 1.2 Verify candidate GEO series

Date: 2026-11-02. Repo commit before the step: 6cb6398. Env: system python3.
Reconstructed on 2026-11-06 from the script's usage line and the accession order
in its output file.

```bash
cd "$REPO"
python3 -I workflow/scripts/stage1_verify_geo.py \
  GSE168243 GSE201654 GSE201586 GSE249746 GSE189501 GSE169301 GSE310000 \
  GSE296776 GSE197289 GSE285983 GSE285984 GSE273557 GSE295206 GSE243073 \
  GSE215994 GSE78150 GSE78148 GSE77968 GSE148059 GSE24290 GSE255436 \
  GSE165569 GSE59739 GSE139088 GSE154659 GSE175421 GSE182098 GSE142541 \
  > data/metadata/geo_verification.tsv
```

Output: `data/metadata/geo_verification.tsv`.
Check: `wc -l data/metadata/geo_verification.tsv` gives 29 (header plus 28);
`cut -f2 data/metadata/geo_verification.tsv | sort | uniq -c` shows 28 `found`.
Notes: GSE189501 and GSE165569 have an empty `pubmed_ids` column. GEO sample
counts can change if submitters update a series, so a rerun may differ slightly.

### 1.3 Check access terms of the non-GEO sources

Date: 2026-11-06. Repo commit before the step: 35a2c7f. Env: system python3.

```bash
cd "$REPO"
python3 -I workflow/scripts/stage1_verify_nongeo.py > data/metadata/nongeo_access.tsv
```

Output: `data/metadata/nongeo_access.tsv`, 5 rows (2 SPARC, 1 Dryad, 2 painseq
Dropbox links). Check: `wc -l data/metadata/nongeo_access.tsv` gives 6.
What it queries: DataCite (`https://api.datacite.org/dois/<doi>`), Pennsieve
(`https://api.pennsieve.io/discover/datasets/<id>`), Dryad
(`https://datadryad.org/api/v2/datasets/<doi>` and its `/files`), and a HEAD
request on each Dropbox link for the size.
Key results on this date: SPARC 465 CC BY 4.0, public; SPARC 464 embargoed until
2027-04-29; Dryad CC0, 12,412,536,831 bytes; painseq Seurat object 1,915,388,302
bytes; painseq harmonized folder 19,935,643,591 bytes.

### 1.4 Look for a data use agreement on the embargoed SPARC dataset

Date: 2026-11-06.

```bash
curl -s -o /dev/null -w "%{http_code}\n" \
  "https://api.pennsieve.io/discover/datasets/464/data-use-agreement"
```

Output on this date: `204` (no agreement attached).

### 1.5 Count donors in SPARC 465 (sympathetic ganglia)

Date: 2026-11-06. Downloads one 16 KB file.

```bash
cd "$REPO"
WORK=$(mktemp -d)
curl -s -o "$WORK/subjects.xlsx" \
  "https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com/465/files/subjects.xlsx"
shasum -a 256 "$WORK/subjects.xlsx"
python3 -I workflow/scripts/sds_subjects_summary.py "$WORK/subjects.xlsx"
```

Check: sha256 is
`649b8a260f34b341477dd8b84519433448b378098f866bb6f67f217cddf5016d`, the same as
the `sha256` Pennsieve gives for `files/subjects.xlsx`
(`curl -s "https://api.pennsieve.io/discover/datasets/465/versions/1/files?path=files/subjects.xlsx"`,
base64 `ZJuKJg80s0FHfdi4RRlDNEizeAmPhmu29n8hfN31AW0=`).
Output: 7 donors, 5 male, 2 female, 23 to 59 years.

To list the files in the dataset:

```bash
B="https://api.pennsieve.io/discover/datasets/465/versions/1/files/browse"
for p in "" files files/primary files/primary/filtered_feature_bc_matrix; do
  curl -s "$B?path=$p&limit=100" | python3 -c '
import sys, json
for f in json.load(sys.stdin).get("files", []):
    print(f.get("type"), f.get("path"), f.get("size"))'
done
```

### 1.6 Read the painseq portals

Date: 2026-11-06. Manual, in a browser.

1. Open https://humandrg.painseq.com (redirects to
   https://painseq.shinyapps.io/u19humandrgatlas/). The page links the preprint
   (10.1101/2025.11.05.686654) and the Seurat object
   `seurat_object_preprint_release.Rds` on Dropbox. No licence or terms of use are
   shown.
2. Open https://harmonized.painseq.com (redirects to
   https://painseq.shinyapps.io/harmonized_painseq_v1/), "Data download" tab. It
   links a Dropbox folder with the full atlas objects and Table S4, and GEO
   GSE255436. No licence or terms of use are shown.

Licence of the preprint (article only):

```bash
curl -s "https://api.biorxiv.org/details/biorxiv/10.1101/2025.11.05.686654" \
  | python3 -c 'import sys, json; [print(c["version"], c["date"], c["license"]) for c in json.load(sys.stdin)["collection"]]'
```

Output on this date: `1 2025-11-06 cc_by_nc_nd`; the `pubs` endpoint reports no
published version yet.

### 1.7 Dryad file access

Date: 2026-11-06. File downloads through the Dryad API return 401 without an API
token, and a scripted request to the web download link returns 403. Download the
files from https://datadryad.org/dataset/doi:10.5061/dryad.gf1vhhmxq in a browser,
or with a Dryad API token. Not downloaded yet.
