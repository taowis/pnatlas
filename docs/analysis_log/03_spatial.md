# Analysis log, stage 3: spatial mapping

What was found and decided at each step. Commands are in
`docs/runbook/03_spatial.md`.

## 3.1 Environment

cell2location 0.1.5 installs cleanly on top of the stage 2 environment but
fails to import against pyro 1.9.2 (a helper it imports from pyro's
autoguide utilities is gone); pinning pyro to 1.9.0 fixes it and
scvi-tools 1.4.2 still imports. openpyxl is added for the authors' xlsx
sample sheets. Both machines were updated the same way and the lock file
refreshed.

## 3.2 Spatial layer download

What the spatial sources hold, from the GEO file lists and the Box folder
pages, before downloading:

- Heming 2025 sural nerve Xenium (GSE285984): eight sections in one RAW
  tar (9.9 GB) with, per section, the cell feature matrix (H5), the cell
  and nucleus boundaries, the cell table and the transcript table
  (parquet) and the morphology image (OME-TIFF). The sample sheet gives
  the diagnosis per section: two CIDP, two CIAP, two vasculitic
  neuropathy, two controls. This is the main nerve spatial input.
- Yu 2024 DRG Xenium (GSE273557): the GEO copy has the per-region
  transcript tables (four regions, three donors: one lumbar section of
  donor 1, lumbar and two thoracic sections of donor 2 by the sample
  sheet) and a small matrix of 1,415 hand-segmented neurons by 87 panel
  genes with donor, level and position per neuron. The full Xenium
  bundles with all cells are on Dryad and need a browser download (stage
  1, 1.7); the neuron matrix is enough for a neuron-family transfer and
  is used first.
- Zhang 2024 DRG Visium (GSE243073): six donors (DRG 2 to 7), a filtered
  matrix archive and an image per section (0.3 GB). Whether the archives
  hold the spot positions is checked after download.
- Vroman 2023 DRG Visium (GSE215994): two sections of one 96-year-old
  female donor, Space Ranger H5, Loupe file and hires image per section
  (0.8 GB); no positions file is listed, so the positions may have to
  come from the Loupe files.
- Shiers 2025 DRG Visium with diabetes history (GSE295206): 16 sections
  from seven donors in three batches, filtered and raw H5, tissue
  positions, Loupe files and the authors' neuron selections per section
  (9.0 GB); the diabetes status per donor is in a docx on GEO.
- Tavares-Ferreira 2022 DRG Visium (Box): 16 Space Ranger H5 files (two
  sections each of four female and four male donors, ages 24 to 65, two
  sequencing runs), a sample sheet, and a combined Seurat object; the H
  and E images (26.6 GB) and Loupe files (10.4 GB) are left on Box for
  now. The H5 files carry no spot positions; the Seurat object is
  expected to.
- DPN sural nerve Visium (Box): one annotated, Harmony-integrated Seurat
  object of five samples (four type II diabetic amputation cases with
  moderate to severe axonal loss, one control with normal axonal density)
  and a sample sheet; the Loupe files are left.

The Box folders have no API without a token, but the shared-link pages
embed the file identifiers and the shared-file download endpoint serves
the files to a plain HTTP client, so the download is scripted like the
GEO one and recorded in the manifest with the exact URL per file.

All 28 files came down in one run and verify against the recorded md5
sums: 10.6 GB for the Heming Xenium tar, 9.7 GB for the Shiers tar, 0.9
GB for the Yu transcripts, 0.8 GB Vroman, 0.3 GB Zhang, 0.3 GB for the
Tavares-Ferreira H5 files and Seurat object and 0.7 GB for the DPN nerve
object. The server disk has 115 GB left after the download, which is
enough for extracting the tars one dataset at a time but not for the
Box images; those stay on Box until the anatomy checks need them.


## 3.3 Section inputs

What the archives hold, after extraction:

- Heming 2025 Xenium: per section a cell feature matrix with 99 gene
  probes (the other 442 features are control probes and codewords,
  dropped) and a cells table with the centroid in microns, cell and
  nucleus area and the transcript count. 4,784 to 20,286 cells per
  section, and few transcripts per cell (a median well under 50), so
  many cells carry almost no signal; the transfer will have to say how
  many cells it can call at all. The 99 genes are more than the 60 the
  plan requires for a targeted panel. Morphology images (1.4 GB to 2.8
  GB per section) stay in the tar.
- Yu 2024 Xenium: the GEO copy has only per-region transcript tables
  and the neuron matrix (87 genes, 1,415 neurons, positions in the cell
  id); the four regions are one lumbar section of donor 1 and one
  lumbar and two thoracic sections of donor 2. The neuron matrix is
  the stage 3 input; the transcript tables are left in the tar.
- Zhang 2024 Visium: six sections (DRG 2 to 7, ages 47 to 59), Space
  Ranger 1.3.1 filtered matrices with 2,000 to 3,000 spots and the
  36,601-gene 2020-A reference. No positions anywhere (the archives
  hold matrix, features and barcodes only; the paper gives GEO as the
  only deposit and Loupe for the image alignment), so these sections
  can take the deconvolution and the marker and held-out gene checks
  but not the anatomy check.
- Vroman 2023 Visium: two serial sections of one 96-year-old female
  donor, filtered H5 (2,005 and 1,775 spots, a 37,481-gene reference)
  and the high-resolution image; the only positions are inside the
  Loupe files, which no open reader parses. Same status as Zhang: no
  anatomy check.
- Shiers 2025 Visium: the docx lists 16 sections from six donors (five
  men aged 34 to 55, one woman of 48) in three batches run with Space
  Ranger 1.1.0, 1.3.0 and 2.0.0; it carries no diabetes status, which
  has to come from the paper. The batch 1 sections M2-1 and M2-2 are
  byte-identical to Tavares-Ferreira VM2-1 and VM2-2 (same md5; 58,288
  features against 36,601 for the other 14), so donor DN0003 is the
  Tavares-Ferreira donor M2 and is used once. Each section has a
  positions file (two header styles) and the authors' Loupe selections
  of Nageotte nodules and the neurons around them, kept as
  author_selection. Images and raw matrices are on the SPARC portal
  rather than GEO and are not needed yet.
- Tavares-Ferreira 2022 and the DPN nerve object: Seurat objects, read
  through the R export in the next step (counts, metadata, coordinates
  per image).

Section identifiers are the authors' own where they are short (S01,
DRG2, DRG2a, Region_1) and a donor plus slide code for Shiers. The
batch is the sequencing centre for Heming (three centres), the Space
Ranger batch for Shiers, and the donor otherwise.

The Shiers 2025 paper (Nature Communications, 2025) states that the 16
Visium sections come from six donors, five with diabetic peripheral
neuropathy and one diabetic donor taking analgesics, with no non-diabetic
donor in the Visium set; the donor numbers in the docx match those in the
paper, so the conditions are written per donor in the config. The
analgesic-taking diabetic donor (14 in the paper, DN0003 in the docx) is
the Tavares-Ferreira donor M2, which means two of the Tavares-Ferreira
"control" sections are from a diabetic donor; that is noted for the
stage 3 interpretation and the sections are read once, from the
Shiers copy, where the positions and the Nageotte annotations are.

The two Seurat exports:

- The Tavares-Ferreira combined object (Seurat 3.2.2) holds 3,952 spots
  from eight samples, one per donor, with 26 neuron-subtype scores and
  no images or coordinates: it is the authors' set of neuron-rich
  barcodes, not the sections. The sections are therefore read from the
  16 Space Ranger H5 files, 1,973 to 3,001 spots each, without
  positions, and the two M2 sections from the Shiers copy. The Box
  folder has the Loupe files, which would give positions, but no open
  tool reads them; the positions may also sit on the SPARC portal, to
  be checked when the anatomy check needs them.
- The DPN nerve object (Seurat 5.0.2) has six images and 12,591 spots
  with coordinates and scale factors for every spot: four diabetic
  amputation cases (moderate or severe axonal loss) and one control
  with two sections, plus the authors' eight-class cell type labels
  (Schwann, perineurial, VSMC/pericyte, endothelial, ECM fibroblast,
  immune, plasma, adipocyte), kept as author_CellType. This is the
  second nerve spatial input and the only Visium nerve data.

The first build of the five GEO datasets ran before the Shiers
conditions were written, so that dataset is rebuilt with them; the
sections of all seven datasets are then summarised in
`results/spatial/section_summary.tsv`.


## 3.4 Reference signatures

The signatures are one negative binomial regression over the joint
reference (642,605 cells with a named state; the 2,076 "unassigned"
cells are left out), state_final as the label, the donor key (dataset
and donor) as the batch and the study as a further categorical
covariate, with cell2location's default gene filter and training
settings (250 epochs, batches of 2,500, learning rate 0.002). States
with the same name in both compartments (macrophage, endothelial,
mural, adipocyte, the two Schwann states) get one signature each,
which is what "the joint reference" in the plan means; the compartment
of the cells behind each signature is recorded in signature_cells.tsv.

The probe on 20,000 cells and 5 epochs ran in under a minute and
peaked at 25 GB of memory, almost all of it the reference counts
(float64, 1.1e9 non-zeros, read whole by anndata before the subset).
The reader was changed to take the arrays with h5py as float32 and
subset before building the object, so the full run stays well under
the server's memory with the guard running.

The full run loads the reference in about a minute (11.5 GB resident,
17,661 genes pass the cell2location filter) and trains at about 17 s
per epoch on the GPU, 5.6 GB of GPU memory.

The plan's check (a) needs "the state's top 20 reference markers" and
stage 2 did not produce a marker table, so one is made here from the
reference itself: a one-vs-rest Wilcoxon test on a stratified sample of
up to 2,000 cells per state, markers being genes with adjusted p below
0.01 detected in at least a quarter of the state's cells, ranked by log
fold change, top 20. The markers come from the same cells as the
signatures but through a different statistic, so the check asks whether
the spots that the deconvolution gives to a state are the spots that
express that state's discriminating genes, which is what the plan
intends. The table is committed with the results.

How the Visium mapping is run per section (workflow/scripts/c2l_map.py):
the genes shared between the section and the signatures are split with
the project seed, a tenth held out; the rest is fitted with
Cell2location at the plan's 30 cells per spot and detection alpha 20,
one section per model (no batch term, since every section is its own
model), 30,000 epochs on the whole section; the posterior 5th percentile
of the abundance is the mapped abundance and the state with the largest
abundance is the spot's dominant state. Check (a) is the AUROC of the
marker score of each state in its dominant spots against the other
spots; check (b) is the Spearman correlation over the held-out genes
between their mean normalised expression in those spots and the state's
signature. A state is checked in a section only with at least 20
dominant spots, and the per-section rows are collected across sections
for the "not mapped" rule of the plan.

The full signature run took 72 minutes on the GPU (17 s per epoch, 250
epochs, 5.6 GB of GPU memory) and peaked at 35 GB of system memory
during the posterior export, which brought free memory down to 23
percent for a few minutes (guard warnings, no pause). The training ELBO
is flat from about epoch 100 on (0.11 percent change over the last 50
epochs), so the default 250 epochs are enough. The signatures cover
17,661 genes and 31 states, and their most specific genes are the
expected ones: SLC2A1, SLC22A3 and XACT for the perineurial state, MPZ
and PRX for myelinating Schwann cells, NPY, TH, DBH, SLC6A2 and TFAP2B
for the sympathetic neurons, MRGPRX1, SST and NPPB for C-NP, MYH11 and
ACTA2 for mural cells, PAMR1 and FGFR2 for endoneurial fibroblasts.
Six state names are shared by the two compartments and carry one
signature each (signature_cells.tsv gives the cell counts per side;
the nerve side dominates mural, non-myelinating Schwann and
endothelial, the ganglion side macrophage). The 20,000-cell probe
directory was deleted on the server after the full run.

## 3.5 Mapping the sections

The targeted-panel transfer was first run as planned, with scVI and
scANVI trained on the shared genes of the compartment's reference cells
and scArches for the query. For the Yu 2024 neurons (85 shared genes,
1,415 cells) that worked: 90 to 95 percent of the neurons assigned with
a median probability of 0.85 to 0.95, every checked family passing both
checks. For the Heming 2025 sural nerve Xenium (96 shared genes, 84,591
cells) it collapsed: 93 percent of the cells were called "mast" with
probability 1.0, no cell was called a Schwann cell, and the cells with
Schwann markers (SCN7A, S100B) sat in the "mast" bucket. The reason is
the data: 83 percent of the cells have fewer than 50 transcripts over
the panel (median 8 to 34 per section), so the query batch differs from
the reference nuclei mostly in depth, and the batch term of the query
model absorbs the biology while the classifier head defaults to the
state with the lowest panel expression. That attempt is kept in the log
and in the history of workflow/scripts/panel_transfer.py; the panel
transfer was replaced by a classifier on the reference signatures
themselves (plan amendment 2): the multinomial likelihood of each cell's
counts under every candidate state's signature over the shared genes,
flat prior, the posterior as the label probability with the plan's 0.5
threshold, and a minimum of 20 transcripts over the shared genes below
which a cell is unassigned for low counts. Because 10 percent of 96
genes is too few for a stable rank correlation, check (b) for panels
rotates over ten folds so that every shared gene is held out once and
the correlation runs over all shared genes; the threshold is unchanged.
With this the nerve sections show the expected picture (Schwann states,
perineurial, endoneurial and epineurial fibroblasts, mural and
endothelial cells), 16 to 63 percent of the cells assigned per section
and the rest almost all below the count floor. The checks of this run
are reported with the collected results.

Panel transfer results. The panel markers for check (a) are now the
reference markers computed within the panel's genes by the same script
as the genome-wide markers (one-vs-rest against all 31 states, so a
state's panel markers are the panel genes that single it out in the
whole reference, not only within its compartment, as the first attempt
had done); 2 to 20 markers per state for the Heming panel, 1 to 20 for
the Yu panel. The outputs of the scANVI attempt on the server (the
per-cell tables under the transfer directory and the three result
tables per dataset) were deleted and replaced by this run's. Heming
2025 sural nerve, 84,591 cells in 8 sections, 96 shared genes: 18 to
69 percent of the cells assigned per section (lowest in the two control
sections S22 and S30, which have the sparsest cells), 99 percent or
more of the unassigned cells are below the count floor and the rest
ambiguous; median posterior 0.90 to 1.00 per section. The
assigned cells are mural (9,637), damage Schwann (5,505), epineurial
fibroblast (5,213), non-myelinating Schwann (5,104), myelinating Schwann
(4,394), perineurial (3,763), endothelial (3,480), endoneurial
fibroblast (3,425), adipocyte (2,041), mast (700), macrophage (289),
repair Schwann (185) and lymphatic endothelial (149). Twelve of the
thirteen states with 20 or more cells in a section pass check (a) in
every section (median AUROC 0.76 for damage Schwann to 1.00 for
endothelial and lymphatic endothelial) and pass check (b) in at least
five of their sections (median Spearman 0.51 for lymphatic endothelial
to 0.77 for myelinating Schwann; endoneurial fibroblast and macrophage
each fail (b) in two sections, adipocyte and lymphatic endothelial in
one). Mast fails (b) in five of its seven sections (median 0.49) and
(a) in two (median AUROC 0.74 on two panel markers), so by the plan's
rule mast is not mapped in the sural nerve Xenium: the panel carries
two genes that mark mast cells, and 700 cells over eight sections is a
thin signal. Yu 2024 DRG neurons, 1,415 hand-segmented neurons in 4
sections, 85 shared genes: every neuron assigned (none below the count
floor, median posterior 1.0), seven neuron families with 20 or more
cells in a section, all passing both checks in every section (AUROC
0.81 to 1.00, Spearman 0.80 to 0.89). The run takes seconds on the CPU;
the same command on the laptop gives identical tables for Yu and tables
that differ from the server's in the sixth decimal of the AUROC for
Heming, with the same pass and fail flags.

Running the Visium sections: the single loop took about 17 minutes per
section with the GPU two thirds busy, so the remaining sections run in
two loops side by side, with a guard that pauses our jobs when another
person's compute process appears on the GPU or memory runs low (details
in the runbook). The guard's first version counted a graphics process
(someone's browser) as GPU work; it now counts compute processes only.
