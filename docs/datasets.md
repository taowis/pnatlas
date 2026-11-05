# Datasets: verification and access

Stage 1 record of what each candidate dataset is, where it lives and under what
terms it can be used. The candidate list is in `docs/novelty_check.md`, section 4.
Where this file and the novelty check disagree, this file is the later check.

## GEO series

Checked with NCBI E-utilities (`workflow/scripts/stage1_verify_geo.py`,
output `data/metadata/geo_verification.tsv`). All 28 candidate series exist and
their titles and sample counts match the novelty check. GSE189501 and GSE165569
have no linked PMID in GEO and still need a manual look.

## Non-GEO sources

Checked with DataCite, the Pennsieve discover API (SPARC), the Dryad
API v2 and HEAD requests to the Dropbox links (`workflow/scripts/stage1_verify_nongeo.py`,
output `data/metadata/nongeo_access.tsv`). Nothing was downloaded except the SPARC
`subjects.xlsx` (16 KB) used to count donors.

| Dataset | Location | Licence | Status | Size | What is there |
|---|---|---|---|---|---|
| Human sympathetic ganglia snMultiome (Yang 2025 preprint, PMID 41000833) | SPARC 465, DOI 10.26275/oqxx-si9h | CC BY 4.0 | Public; anonymous download from the SPARC S3 bucket works | 19.9 GB, 22 files | 7 donors (5 male, 2 female, 23 to 59 years). Cell Ranger ARC `filtered_feature_bc_matrix` (2.2 GB), `meta.data.csv` (12 MB), `atac_fragments.tsv.gz` (17.8 GB), SDS subject and sample sheets |
| Human DRG snMultiome (Zhao/Boyer 2025 preprint, PMID 41040282) | SPARC 464; DOI 10.26275/sojr-qhvm in the novelty check, Pennsieve now reports 10.26275/z6bl-slmq; both resolve to dataset 464 version 1 | CC BY 4.0 once released | Embargoed until 2027-04-29 | 237 GB, 1,038 files | Not available now. No data use agreement is attached in Pennsieve (API returns 204). Earlier access would need a request to the data owners through the SPARC portal; the terms of such a request were not checked |
| Human DRG Xenium (Yu 2024 Nat Neurosci, PMID 39496796) | Dryad, DOI 10.5061/dryad.gf1vhhmxq (version 3, 2024-10-25) | CC0 1.0 | Public. The novelty check marked this DOI as not resolved; it does resolve | 12.4 GB, 5 files | Four Xenium output bundles (`Region_1` to `Region_4`, `.rar`), 100-gene custom panel, 1,340 manually segmented neurons, plus `README.md`. Same four samples as GEO GSE273557. File download through the API needs a Dryad API token (401 without one) and scripted web download returned 403, so download in a browser or with a token |
| Human DRG reference atlas, preprint release (Bhuiyan 2025, PMID 41279342) | humanDRG.painseq.com (redirects to painseq.shinyapps.io/u19humandrgatlas), Seurat object `seurat_object_preprint_release.Rds` on Dropbox | None stated on the portal. The preprint itself is CC BY-NC-ND; that covers the article, not the object | Downloadable without login. No journal version yet (bioRxiv API). Raw data promised for dbGaP/SPARC | 1.92 GB | Integrated Seurat object with donor_ID, SpinalRegion, Level, Sex, Age, Platform and atlas annotation |
| Harmonized DRG/TG atlas (Bhuiyan 2024 Sci Adv, PMID 38905344) | harmonized.painseq.com (painseq.shinyapps.io/harmonized_painseq_v1), Dropbox folder from the "Data download" tab | None stated on the portal | Downloadable without login. GEO GSE255436 holds only the new axolotl data | 19.9 GB (folder) | Full neuronal and non-neuronal objects per atlas (for example 75,928 DRG neurons from 21 datasets) and Table S4 |

## What this means for the plan

- The DRG backbone for layer 1 can be the Bhuiyan 2025 preprint object. It is
  enough for analysis with citation, but because no licence is stated, the object
  and anything derived cell by cell from it (re-annotated objects, portal exports)
  should not be redistributed on Zenodo or the portal without written permission
  from the authors. Summary statistics (per cell type means, marker tables) are a
  lower risk but should still be cleared with them. The object will change when the
  paper is published, so record its sha256 and download date.
- The Zhao/Boyer DRG multiome cannot be used before 2027-04-29 unless the owners
  grant access. It is not needed for the main line; leave it out of the frozen plan
  and note it as a possible addition.
- The sympathetic ganglia multiome is open under CC BY and is the only public
  human sympathetic ganglion dataset found. It is relevant to the autonomic genes in
  HSAN and fits the ganglion side of the reference.
- The Yu Xenium data are CC0 and can be used and redistributed freely. With only
  100 genes and 4 sections it is a small spatial layer; the Heming sural nerve
  Xenium (GSE285984) and the GEO Visium series remain the main spatial inputs.

Commands for every check above are in `docs/runbook/01_datasets.md`; what was found at each step is in `docs/analysis_log/01_datasets.md`.

## Still to check

- Zenodo 10.5281/zenodo.14226218 (Heming 2025) and sensoryomics.com (Tavares-Ferreira Visium).
- Manual look at GSE189501 and GSE165569.
- Read the Dryad `README.md` after downloading it in a browser.
