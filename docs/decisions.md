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
- **One cell2location signature per state name across compartments.** The
  plan says the signatures come from the joint reference, so a state that
  exists in both compartments (macrophage, endothelial, mural, adipocyte,
  myelinating and non-myelinating Schwann) gets one signature from all its
  cells; the compartment of the cells behind each signature is recorded.
  Splitting them would give a Visium DRG section two macrophage states to
  divide the same spots between.
- **Reference markers for check (a) are a one-vs-rest Wilcoxon test on a
  stratified reference subsample, top 20 by log fold change among genes
  detected in a quarter of the state's cells.** Stage 2 produced no marker
  table and the plan does not define one. For a targeted panel the test is
  run within the panel's genes, since the genome-wide top 20 of a state
  rarely overlaps a 100-gene panel; the plan's thresholds are kept.
- **Targeted panels are transferred with scVI and scANVI trained on the
  shared genes of the compartment's reference cells, then scArches for the
  query.** The stage 2 models use 3,000 highly variable genes, of which a
  panel shares few, so a model on the panel's genes is trained for each
  panel dataset; for the Yu 2024 neurons only the reference neurons are
  used, because the query holds only hand-segmented neurons.
- **Targeted panels are classified by the multinomial likelihood of each
  cell under the reference signatures over the shared genes (replaces the
  entry above; plan amendment 2).** On the Heming 2025 sural nerve Xenium
  the scANVI transfer put 93 percent of the cells into "mast" with
  probability 1.0 and no cell into a Schwann state: 83 percent of the
  cells have fewer than 50 transcripts over the panel, so the query batch
  differs from the reference nuclei mostly in depth, and the batch term
  absorbs the biology. The signature classifier has no batch term to
  misuse, its parameters are the signatures already estimated for the
  Visium mapping, and a count floor (20 transcripts over the shared
  genes) states openly which cells are too sparse to call instead of
  forcing a label on them. The candidate states stay within the
  compartment and the authors' class, as in the stage 2 transfer. Check
  (b) rotates over ten folds of the shared genes because a tenth of a
  96 gene panel gives a rank correlation over 10 genes, which moved by
  more than 0.3 between seeds; with every gene held out once the
  correlation runs over all shared genes. The Yu 2024 neurons, which
  passed under scANVI, pass under the classifier as well.
- **In a Visium section a state is checked on the spots where it is the
  dominant state by abundance share, and only when it dominates at least
  20 spots.** The plan's checks (a) and (b) are written for "spots or
  cells assigned to that state", and cell2location assigns a spot a
  mixture, not a label. The dominant state is the nearest thing to an
  assignment that does not need a threshold on the abundance scale, and
  20 spots is the floor below which an AUROC or a rank correlation says
  little. The cost is that states that are present but never the most
  numerous in a 30 cell spot, the DRG neuron families above all, go
  unchecked in Visium and are checked in the Yu 2024 Xenium neurons
  instead; the abundance estimates themselves are kept for every state.
- **Yu 2024: every segmented cell of the Dryad bundles is a second dataset
  (`yu2024_hdrg_xenium_cells`); the hand-segmented neuron table stays as it
  is.** The GEO table holds the authors' curated 1,415 neurons, whose checks
  are recorded; the Dryad bundles hold the 27,000 to 46,000 segmented cells
  per region that the table leaves out, so the non-neuronal ganglion states
  can be checked in this dataset and the anatomy panels of check (c) drawn
  on the post-Xenium H and E image. Keeping both ids leaves the recorded
  neuron checks untouched and keeps the two questions apart.
- **For a whole DRG section the three nerve sheath states (perineurial,
  epineurial fibroblast, endoneurial fibroblast) are candidates besides the
  ganglion states.** The compartment rule alone would offer none of the
  states check (c) asks about (the ganglion datasets label their sheath
  cells as fibroblast), so the Yu panels would be empty by construction,
  while a DRG section carries its capsule and the root entering it, which
  are nerve sheath tissue. The addition is per dataset (`extra_states` in
  the config); the cost, a ganglion fibroblast called a sheath state where
  the signatures are close over the panel's genes, is what checks (a) and
  (b) and the panels are there to show.
- **The sympathetic neuron state is not a candidate on a DRG section.** A
  DRG holds no sympathetic neurons, yet the state took 15,593 of the
  152,000 segmented cells of the Yu 2024 sections (analysis log 3.9), most
  likely satellite glia and other cells with thin, neuron-like panel
  counts, and it fell inside the plan's mapped rule only because it
  failed (a) in two sections of four. The removal is per dataset
  (`exclude_states` in the config) so that the sympathetic ganglion
  sections of a later dataset keep the state; the cells it held are
  reassigned among the remaining candidates or become unassigned, and
  the checks of those states are recorded again.
- **Check (c) is read by the analyst, on the panels, with a written
  rule.** The plan asks for one reader and a yes or no per section; the
  rule (runbook 3.11: where each sheath state must lie in a nerve and in
  a ganglion, "no" for a state spread without regard to the structure,
  "not present" below about 20 cells) is written down so that a second
  reader can repeat the reading, and every row carries a note of what
  was seen. The result is reported next to the mapped verdict and does
  not enter it, because the plan's rule names (a) and (b) only; changing
  the rule after seeing the result would be an amendment made in
  hindsight.
- **The UT Dallas spot coordinates are read as row,column.** The Seurat 5
  `GetTissueCoordinates` convention (x is the image row) put the spots
  beside the tissue; the orientation was chosen by scoring all eight
  against the tissue mask, not by eye, and is recorded in the config so
  the export itself is unchanged.
- **MuSiC is the primary deconvolution method by the rule fixed before
  the run, and BayesPrism's estimates are kept beside it for the rare
  states.** Only MuSiC meets the pseudobulk thresholds (analysis log
  4.3), so the rule names it; but it returns exact zeros for adipocyte,
  granulocyte, mast, lymphatic endothelial and repair Schwann in most
  GTEx samples, while BayesPrism's fractions for those states follow the
  marker genes in the bulk (Spearman 0.6 to 0.9). Changing the rule
  after seeing the result would be an amendment in hindsight, so the
  rule stands; the rare states are read from BayesPrism, and both tables
  are released.
- **Validation (b) reads the fat percentage in the GTEx pathology notes
  and is reported as uninformative.** A mention of fat is in 579 of 670
  notes, which leaves no contrast; the percentage gives one, but the
  adipocyte marker genes themselves follow it at Spearman 0.15 only,
  because the notes describe the piece before the fat was trimmed. The
  ceiling is kept in the table and the marker concordance is the
  evidence that the estimates carry signal.
- **The damage Schwann state is not interpreted in the bulk.** Neither
  method recovers it in the mixtures (Pearson about -0.08); it is kept
  in the tables and left out of the reading.
- **The back-test's ten seeds are PNPrank's global seed and the next
  nine, applied to its cutoff scoring; the published numbers are the
  first seed.** PNPrank's own back-test runs one seed per cutoff and
  the ten seeds of its config serve its repeated cross-validation; the
  plan names ten fold seeds and a noise baseline over them, so the same
  scoring function is run with those ten seeds for both feature sets.
  Nothing in the scoring, transform or evaluation code is changed.
- **The gene by state correlation feature uses the BayesPrism
  fractions.** MuSiC is the primary deconvolution method by the
  pre-set rule, but its exact zeros for the rare states in most samples
  leave nothing to correlate with; BayesPrism's fractions vary in every
  state and follow the marker genes (analysis log 4.3).
- **The gene class enrichment is not a PNPrank feature.** It is
  computed from the truth set and would leak the labels into the
  back-test; the colocalisation, specificity, compartment and
  deconvolution columns do not use the truth set.
- **The portal is a static site built from the result tables and is not
  committed.** Every value on it is read from `results/` at build time,
  so the page needs no server and can stay on GitHub Pages without
  maintenance; the pages workflow rebuilds it on pushes to `main`, and
  committing the built JSON (20 MB per build) would only duplicate the
  tables in the history.
- **The Bhuiyan 2025 cells are left out of the archived reference
  object.** The humanDRG portal states that cell-level derivatives are
  not to be redistributed; the state-level tables, which aggregate over
  many cells, keep every dataset. A user who needs the full object
  rebuilds it from the runbook with the portal download.
- **The public repository and the Zenodo upload wait for the
  manuscript.** Opening the repository before submission would amount to
  a preprint, and the archive version cited in the paper should be the
  one the final tables come from; the build and staging scripts are
  rehearsed now so that the release is a single command then.
