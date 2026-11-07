# Decisions

Each non-obvious choice, with its reason. Newest last.

## Stage 0 (scaffold, 2026-11-02)

- **Scope: peripheral nervous system first, MASLD/MASH on hold.** A literature
  and data check on 2026-11-02 found each MASH layer already covered by
  2025-2026 papers, the largest staged single-cell cohorts under controlled
  access and only about 33 open spatial MASH biopsies. For ganglia and nerve,
  harmonised DRG atlases exist (Bhuiyan 2024, PMID 38905344; Bhuiyan 2025
  preprint, PMID 41279342), but no joint ganglia + nerve reference, no GTEx
  tibial nerve deconvolution with a nerve reference and no cell-state nerve
  eQTL colocalisation were found. Details in `docs/novelty_check.md`.
- **Separate from PNPrank.** PNPrank is submitted on its own; pnatlas is a
  follow-up that uses PNPrank as an application (stage 6). Merging them would
  hold back a finished manuscript.
- **PNPrank as a dependency, not copied code.** Added at stage 6, pinned to
  the PNPrank release tag used for submission, so the back-test code and
  splits are identical between versions.
- **Reuse processed atlases** from the PNPrank data root (`pnprank_root` in
  `config/paths.yaml`) instead of downloading them again.
- **Nerve layer starts from Heming 2025** (Nat Commun, PMID 40849297; GEO
  GSE285983/GSE285984), building on an earlier local re-analysis.
- **Workflow engine: Snakemake**, same layout as PNPrank (one rules file per
  stage). Packaging with `pyproject.toml`, `src/` layout, dependencies in
  `environment.yml`, `pip install -e . --no-deps`.
- **Environment:** spatial and deconvolution tools (cell2location, squidpy,
  scvi-tools) are added at the stage that needs them, possibly as a separate
  environment because of their torch dependency.
- **Gadi:** heavy integration runs on NCI Gadi. A Snakemake PBS profile is
  added with the first Gadi job, using project qn73 in committed files.
- **Random seed:** `20261102` (scaffold date).
- **Disk budget:** 300 GB under `data_root` (same as PNPrank); ask before
  exceeding.
- **Private notes:** `private_notes/` is kept locally only and is
  gitignored here.
- **Runbook and analysis log:** every command and manual step is written in
  `docs/runbook/` and what was found and decided at each step is written, in
  English, in `docs/analysis_log/`, so the analysis can be redone by hand and
  its reasoning followed. Command output is not kept as a record. Both are part
  of the public record and go on the public release whitelist. They carry no
  dates or times of the work.
## Stage 1 (data inventory and access)

- **DRG backbone from the Bhuiyan 2025 preprint object, not rebuilt.** It has
  126 donors; rebuilding would duplicate it. Its licence is not stated, so
  nothing cell-level derived from it is redistributed without author
  permission.
- **Integration by scVI/scANVI with Harmony as the comparator, chosen by a
  fixed scIB rule.** Written into the plan before any embedding is computed so
  the method is not picked on downstream results.
- **cell2location for Visium, scANVI label transfer for targeted panels.**
  Targeted panels (100 to 472 genes) are too sparse for a deconvolution model;
  label transfer on shared genes is the standard route.
- **BayesPrism primary, MuSiC secondary, validated on pseudo-bulk mixtures and
  GTEx pathology notes before any GTEx result is read.** Both are reference
  based; the validation thresholds are fixed in the plan.
- **Mendelian gene classes from the PNPrank truth set and PanelApp (GE 846
  v8.38, GE 85 v1.513, AU 3070 v3.1, AU 3439 v2.0); OMIM not used.** OMIM
  needs a licence key; PanelApp has an open API and version numbers.
- **GWAS limited to studies with full summary statistics in the GWAS Catalog.**
  Colocalisation needs the full statistics; studies with top hits only
  (GCST90104786, GCST008223, GCST008530) are listed but not used.
- **coloc priors p12 5e-6, PP.H4 above 0.8.** The default p12 of 1e-5 is known
  to be permissive; 5e-6 is the more conservative value in common use.
- **Cell-state interaction eQTL only as an amendment.** It needs GTEx donor
  genotypes under dbGaP phs000424; the main line does not wait for that.
- **Stage 6 improvement rule: beat the across-seed range of the published
  features and a paired permutation test over cutoff years.** A single
  AUROC difference inside the seed-to-seed spread is noise; the rule is fixed
  before the run.

## Stage 2 (joint ganglia and nerve reference)

- **Compute placement: data download and integration on the server, Gadi for
  the CPU-parallel stages, the laptop for development only.** The laptop
  (16 GB, no CUDA) and the ANZAC desktop (4 cores, 15 GB, no GPU, no conda)
  cannot run scVI on the full reference or cell2location; the server has 16
  cores, 60 GB and a 32 GB GPU, so stages 2 and 3 run there. Its disk is the
  limit (about 180 GB free), so only the files a stage needs are downloaded
  and raw downloads can be staged on the desktop's large disk on the same
  network. Gadi is kept for stages 4 and 5 (BayesPrism over 667 GTEx samples,
  colocalisation over many loci) and as the fallback when the server's GPU is
  in use by the sequencer, because its compute nodes have no internet, its
  scratch is purged after 100 days and its allocation is shared. The laptop
  keeps the same environment for code development on subsets, documentation
  and pull requests. In the repository the server is named only as "the
  server"; its address is passed at run time through `REMOTE`.
- **The painseq object is a label source, not the donor backbone.** Its
  metadata has no donor for most cells and it already contains the neurons of
  three datasets read separately here, so counting it as 126 donors would
  double count cells and inflate donor numbers. Jung 2023 becomes the DRG
  non-neuronal backbone; neuron subtypes follow the atlas labels. Recorded as
  plan amendment 1.
- **QC: author-filtered datasets keep all their cells, flags recorded only;
  Yu 2024 somata added to that group after the first QC run.** Refiltering
  annotated data would drop cells the labels and the label-transfer
  comparison depend on. The whole-cell mitochondrial threshold of 20
  percent is wrong for hand-picked neuronal somata (median 24 percent) and
  would have removed 71 percent of the authors' final set; the change was
  made after seeing the data and is stated in the analysis log (2.4).
- **Merge on the genes shared by every dataset (19,420 symbols), not the
  union.** The eight datasets were annotated with different gene sets
  (32,214 to 61,817 symbols; union 88,803) and the extra names are version
  specific lncRNA and clone-based identifiers. Filling them with zeros would
  make missing annotation look like absent expression in some studies and
  bias the batch models; the shared set holds the protein-coding genes and
  the common non-coding ones, which is what highly variable gene selection
  draws from anyway. The per-dataset QC files keep the full gene sets.
- **Batch unit for integration: donor where known, library otherwise.** The
  plan names the donor as the batch key. The painseq atlas has no donor and
  for most of its cells no library name either (orig.ident "SeuratProject"),
  so those cells take their atlas source study (U19_UTD, U19_HMS and so on)
  as the unit. Study is a separate categorical covariate in scVI and a
  second Harmony variable, so the batch unit never has to encode the
  dataset.
- **scANVI seed labels are the mapped states, with "unassigned" states and
  unlabelled datasets treated as unknown.** The Jung 2023 neurons carry only
  "DRG Neuron", and the trigeminal vascular and immune cells are not split
  by the authors; seeding scANVI with those would teach it a label that
  means "not annotated". Nguyen, Zhang and Yu have no author labels. These
  cells get predictions instead and are left out of the scIB metrics, which
  need a label per cell.
- **scIB metrics on 50,000 labelled cells per subset.** The silhouette
  metrics are quadratic in cell number; the nerve subset alone is 365,708
  nuclei. A fixed-seed subsample of 50,000 cells per subset is used for all
  three embeddings alike, so the comparison is fair even if the absolute
  values move a little with the sample. Batch ASW is the bras metric of
  scib-metrics 0.5, the library's batch silhouette.
- **Yang 2025 sympathetic ganglion author labels are not used.** The labels
  in the portal's metadata table do not follow the expression of the cells
  they are attached to (adjusted Rand index 0.005 against a clustering of
  the dataset, chance-level neighbour agreement, marker-positive cells
  labelled at the dataset's base rates), although the barcodes and per-cell
  counts match the matrix. Using them would seed scANVI with wrong labels
  and score the ganglion metrics against noise. The dataset enters
  unlabelled and is annotated by label transfer; the authors' labels stay
  in `author_label` for the record. Recorded in analysis log 2.5.
- **Label transfer is a nearest-neighbour vote in the scANVI space, within
  the cell's compartment and, where the authors gave a class, within that
  class.** scANVI's own predictions cross compartments (Jung 2023 neurons
  predicted as mast cells, Zhang glia as nerve fibroblast states), because
  the classifier sees every state at once. Restricting the candidate labels
  to the ganglion or nerve states, and to the author's class for the
  "unassigned" cells, keeps the transfer inside what the anatomy allows.
  Fifteen neighbours and a 0.5 vote threshold follow the plan's rule for
  the spatial mapping (label probability below 0.5 is unassigned).
- **Neurons of the sympathetic ganglion dataset are "sympathetic neuron".**
  The reference has no sympathetic seed cells (the authors' labels are
  unusable, see above), so these nuclei can only land among DRG and
  trigeminal neuron states; a sympathetic ganglion has no sensory neurons,
  so any nucleus that sits among neurons is a sympathetic one. TH, DBH and
  PRPH detection in those nuclei is reported as the check and is not used
  to assign, because nuclei with few counts miss any single marker.
- **The state rule counts known donors only.** The painseq atlas has no
  donor, so a DRG neuron family passes the rule only through the Jung,
  Nguyen and Yu neurons that the transfer puts into it. The pooled Yang
  2025 library counts as one donor unit.
