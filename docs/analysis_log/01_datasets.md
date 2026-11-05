# Analysis log, stage 1: data inventory and access

What was found and decided at each step. Commands are in
`docs/runbook/01_datasets.md` under the same numbers.

## 1.2 GEO series verification

All 28 candidate GEO series from the novelty check exist, and their titles and
sample counts match what the novelty check recorded. Two series, GSE189501
(human glia from DRG and spinal cord, 6 samples) and GSE165569 (macaque DRG),
have no PubMed ID linked in GEO and need a manual look to confirm the source
paper. On a rerun GEO listed the taxa of GSE142541 in a different order, a
reminder that GEO records can change and that the stored table is the version
used here.

## 1.3 Access terms of the non-GEO sources

Checked the sources that are not on GEO: two SPARC datasets, one Dryad dataset
and the two painseq portals.

- SPARC 465, human sympathetic ganglia snMultiome (Yang 2025 preprint,
  PMID 41000833): CC BY 4.0, public, downloadable without login from the SPARC
  S3 bucket. 19.9 GB in 22 files: Cell Ranger ARC filtered matrix (2.2 GB),
  `meta.data.csv` (12 MB), ATAC fragments (17.8 GB), SDS subject and sample
  sheets. This is the only public human sympathetic ganglion dataset found, and
  it matters for the autonomic genes in HSAN.
- SPARC 464, human DRG snMultiome (Zhao/Boyer 2025 preprint, PMID 41040282):
  embargoed, release 2027-04-29, 237 GB in 1,038 files. The DOI quoted in the
  novelty check (10.26275/sojr-qhvm) still resolves to dataset 464, but
  Pennsieve now reports 10.26275/z6bl-slmq for the same version. No data use
  agreement is attached, so earlier access would need a request to the owners.
- Dryad 10.5061/dryad.gf1vhhmxq, human DRG Xenium (Yu 2024, PMID 39496796):
  CC0, version 3, 12.4 GB in five files (four Xenium output bundles and a
  README). The novelty check marked this DOI as unresolved; it resolves. These
  are the same four sections as GEO GSE273557, with a 100-gene custom panel
  and 1,340 manually segmented neurons.
- painseq portals: the humanDRG portal (Bhuiyan 2025 preprint, PMID 41279342)
  offers a 1.92 GB Seurat object on Dropbox; the harmonized portal (Bhuiyan
  2024, PMID 38905344) offers a 19.9 GB Dropbox folder with the full atlas
  objects. Neither page states a licence or terms of use. The preprint is
  CC BY-NC-ND, which covers the article, not the data object. There is no
  journal version of the preprint yet; raw data are promised for dbGaP/SPARC.

Conclusions:

- The Bhuiyan 2025 object can be the DRG backbone of the reference, with
  citation, but with no stated licence the object and anything derived from it
  cell by cell (re-annotated objects, portal exports) should not be
  redistributed on Zenodo or the portal without the authors' written
  permission. Summary statistics are lower risk but should also be cleared.
  The object may change when the paper is published, so its sha256 and the
  version used are recorded at download.
- SPARC 464 cannot be used before its embargo lifts and is not needed for the
  main line; it stays out of the frozen analysis plan and is noted as a
  possible later addition.
- The Yu Xenium data are free to use and redistribute, but with 100 genes and
  four sections they are a small supplementary spatial layer. The Heming sural
  nerve Xenium (GSE285984) and the GEO Visium series remain the main spatial
  inputs.

## 1.4 Data use agreement on SPARC 464

Pennsieve returns no data use agreement for dataset 464 (HTTP 204), so the
embargo is the only barrier and the terms of an early-access request were not
checked.

## 1.5 Donors and files in SPARC 465

The SDS `subjects.xlsx` gives 7 donors, 5 male and 2 female, aged 23 to 59.
The downloaded file's sha256 matches the value Pennsieve publishes, so the
anonymous S3 download delivers the published file. The donor counting was first
done with ad hoc code and then moved into `workflow/scripts/sds_subjects_summary.py`
so it can be rerun.

## 1.6 painseq portals and preprint licence

Both portals redirect to shinyapps.io applications. The download links are
public Dropbox links; no licence text appears anywhere on either page. The
bioRxiv API confirms the preprint licence is CC BY-NC-ND and that there is no
published version.

## 1.7 Dryad download

Scripted downloads from Dryad fail: the API returns 401 without an API token
and the web download link returns 403 to a script. The data are open (CC0), so
this is only a mechanics problem; the files will be fetched in a browser or
with a token when stage 3 needs them.

## Open items at the end of this pass

- Zenodo 10.5281/zenodo.14226218 (Heming 2025 processed objects) and
  sensoryomics.com (Tavares-Ferreira Visium) still to be checked.
- Manual look at GSE189501 and GSE165569.
- Dryad README to be read after a browser download.
- `config/datasets.yaml` is still empty and is to be filled from these results.

## 1.8 Heming 2025 on Zenodo

The Zenodo record (concept DOI 10.5281/zenodo.14226218, latest version
15108216, CC BY 4.0) holds 36 `.qs` files of figure data, not the processed
nuclei. The merged object has to be rebuilt from the GEO raw matrices
(GSE285983_RAW.tar plus the two metadata CSVs) or taken from the earlier
re-analysis in `pns-cmt-gene-map-natcomm2025`. The Zenodo files are still
useful for the pseudobulk DE results and the propeller abundance tests, which
can be compared with ours.

## 1.9 sensoryomics.com

The UT Dallas site is a shell around 11 Shiny applications. Three matter here:
`human-peripheral-nerves` (processed data of the DPN tibial and sural nerve
study, PMID 40828619, raw in dbGaP phs001158), `Human-SCG` (the sympathetic
chain spatial preprint, PMID 40667173) and `SC_Visium_Dashboard`
(Tavares-Ferreira DRG Visium). None of the apps can be read by a script; they
have to be opened in a browser to see what can be downloaded. The DPN nerve
paper's data statement confirms that the only open route to those data is the
app, so a dbGaP request stays the way to get nuclei-level data for diseased
nerve.

## 1.10 The two GEO series without a linked paper

GSE165569 is Kupari 2021 (Nat Commun, PMID 33686078); the GEO title is the
paper title and Europe PMC links them. Resolved.

GSE189501 is from Zhang Donghang's group at West China Hospital: snRNA-seq of
glia from DRG and spinal cord of 3 donors, raw UMI matrices on GEO. No paper
could be found; the group's eLife spinal cord atlas cites other series, and the
only citation of the accession is a 2026 PLoS One paper that reused it. It is
usable with a GEO citation, but 3 donors of glia only adds little to the DRG
reference; keep it as optional.

## 1.11 Sympathetic chain spatial study

PMID 40667173 in the novelty check is PubMed's record of the bioRxiv preprint,
not a journal article. Still unpublished.

## 1.12 Dataset table

`config/datasets.yaml` now has 22 entries across the reference, spatial and
bulk layers plus the macaque anchor, each with its accession, access class and
PMID. Mouse anchors and the genetics layer (GWAS, eQTL) are not in it yet.

## 1.13 What the sensoryomics apps actually offer

Opened the apps in a browser. Two of them lead to public Box folders with
processed data that the data availability statements did not mention:

- The DPN nerve app's "Download spatial data" button opens a Box folder with a
  merged, Harmony-integrated and annotated Seurat object of the DPN sural nerve
  Visium data (675 MB), a sample metadata file and five Loupe files, 3.6 GB in
  all. So the processed spatial data for diseased human nerve are open after
  all; only the raw reads need dbGaP. The bulk RNA-seq (6 control vs 6 DPN
  sural; 5 paired sural and tibial DPN) is only viewable gene by gene as qnTPM.
- The human DRG app ("RNA-Data") links a Box folder with the Tavares-Ferreira
  2022 DRG Visium data: 17 Space Ranger H5 files, 17 Loupe files, H and E
  images, the 31 supplementary tables and a combined Seurat object (186 MB),
  34.7 GB in all. This replaces "processed via portal only" in the earlier
  table: the Visium layer from this study is open.
- The sympathetic chain app shows CPM tables and plots per gene and has no
  download; those data stay portal-only until the preprint is published.
- `SC_Visium_Dashboard` is a mouse spinal cord atlas and is not relevant.

Neither Box folder states a licence. As with the painseq objects, they can be
used with citation but should not be redistributed without asking.

Consequence for the plan: the human spatial layer now has open processed data
for healthy DRG (Tavares-Ferreira Visium, 8 donors; Yu Xenium, 4 sections),
diseased sural nerve (Heming Xenium; DPN Visium) and DRG Visium from GEO
(Zhang, Vroman). The dbGaP request is for raw reads only and is not on the
critical path.

## 1.14 Genetics layer

Verified the sources for stage 5 against the live APIs. Five GWAS with full
summary statistics are usable: MVP hereditary and idiopathic peripheral
neuropathy (GCST90475844), MVP type 2 diabetes with neurological complications
(GCST90475676), diabetic neuropathy (GCST90129436), neuropathic pain
(GCST90859373) and carpal tunnel syndrome (GCST90018813), plus multisite
chronic pain (GCST008512) as a secondary trait. The DPN in T2D (ACCORD) and
CMT1A modifier studies have top associations only and are not used. There is
no GWAS for CMT itself, as expected for a monogenic disease. GTEx v10 tibial
nerve cis-eQTL and sQTL tarballs are in the `adult-gtex` bucket (2.56 GB and
1.96 GB); 667 donors with genotypes, 20,101 eGenes, 6,590 sGenes. The eQTL
Catalogue has the GTEx tibial nerve datasets QTD000286 to QTD000290 (532
samples) with SuSiE credible sets. Gene classes come from PanelApp GE panels
846 and 85, AU panels 3070 and 3439, and the PNPrank truth set (242 genes: CMT
59, dHMN 19, HSAN 16, syndromic 124, HSP 24). Two API quirks: the eQTL
Catalogue API v2 is retired (HTTP 410), so the dataset table on GitHub is used;
the GTEx open-access file metadata endpoint returns nothing for v10, so the
bucket listing is used. The per-tissue file names inside the GTEx tarballs
were not checked.

## 1.15 Analysis plan frozen

The plan now fixes, before any stage 2 code is run: which datasets go in each
layer (embargoed and controlled raw data out until granted); donor as the unit
of inference; scVI/scANVI vs Harmony chosen once by a fixed scIB rule;
cell2location for Visium and scANVI label transfer for targeted panels, with
pass/fail checks per section; BayesPrism primary and MuSiC secondary for
GTEx, validated on pseudo-bulk from held-out Heming donors; EWCE-style
enrichment and coloc with SuSiE for the genetics; and the stage 6 comparison
(PNPrank published features vs PNPrank plus pnatlas features, same truth set,
same splits and seeds, mean tier-1 AUROC as the primary metric, improvement
only if it beats the full across-seed range of the baseline with a paired
permutation test over cutoff years). Later changes go in the Amendments
section. One dependency is still open: PNPrank has no submission tag yet, and
stage 6 cannot start until it exists.
