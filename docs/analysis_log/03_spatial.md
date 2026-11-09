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

Outcome of the parallel run: all 44 sections finished with exit 0, the
guard never had to pause anything (no other compute process appeared on
the GPU and free memory stayed above 70 percent), and the memory guard
logged no warning during the mapping. The two loops did not shorten
the total time: one model alone trained in 13 to 19 minutes per
section, two side by side took 33 to 44 minutes each and three took 53
minutes, so the GPU was the bottleneck at about 96 percent utilisation
and the throughput stayed at one section per 17 to 20 minutes. The
loops and the guard are kept because they make a relaunch safe, not
because they are faster here; a second GPU would change that. Nothing
of ours was left running or paused on the server afterwards; the
section outputs (40 MB) stay there and are copied to the laptop's data
root.

Visium mapping results. The ELBO of every section model is flat over
the last tenth of the epochs (relative change below 0.0003). The
estimated number of cells per spot, the sum of the state abundances,
has a median of 42 (Shiers), 44 (Tavares-Ferreira), 35 (Zhang) and 66
(Vroman) in the DRG sections and 24 in the DPN sural nerve sections;
the dominant state holds a median of 29 to 40 percent of a spot's
abundance, so no spot is one cell type and the dominant state is a
label of convenience for the checks. In the DRG sections the dominant
states are satellite glia (31 to 48 percent of spots), endothelial (1
to 29 percent), non-myelinating Schwann (10 to 57 percent), myelinating
Schwann, macrophage and fibroblast; a neuron family is the dominant
state in under one percent of the spots (439 of 48,214 Shiers spots,
mostly A-Propr; 19 Tavares-Ferreira; 7 Zhang; none in Vroman), because
a 55 micron spot holds a few large neurons among dozens of glia and the
dominant state is the one with the most cells. The consequence for the
checks, which run on the spots a state dominates, is that the neuron
families are checked in the Visium sections only where A-Propr reaches
20 dominant spots (five sections, passing both checks) and otherwise
rest on the Yu 2024 Xenium neurons. In the DPN sural nerve sections the
dominant states are non-myelinating Schwann (37 percent of spots),
perineurial (22), endothelial (13), mural (11), myelinating Schwann (9)
and epineurial fibroblast (8), the picture a nerve cross-section should
give.

Mapping status (results/spatial/state_mapping_status.tsv, 28 state and
compartment rows). Not mapped by the plan's rule: non-myelinating
Schwann in the ganglion sections, which is the dominant state in a
third of the DRG spots but fails check (a) in 31 of 33 sections
(median AUROC 0.46; (b) passes in 27 of 33), so the spots the model
gives to this state do not carry its markers above the rest of the
section, and the abundance it receives in the DRG most likely belongs
to satellite glia and the nerve root glia whose signatures overlap it
on a 30 cell spot; mast in the sural nerve Xenium (above); epineurial
fibroblast in the ganglion sections, which reaches 20 dominant spots
in one section only and fails (b) there. Every other state is mapped.
The weaker passes are in the ganglion: endothelial fails (a) in 15 of
36 sections (median AUROC 0.74), satellite glia in 12 of 36 (0.75),
fibroblast in 4 of 9, and myelinating Schwann fails (b) in 10 of 30
sections (median Spearman 0.55); in the nerve every state passes (a) in
every section (median AUROC 0.89 to 0.99 in the Visium nerve sections)
and (b) in all but one or two sections. The held-out gene correlation
sits between 0.5 and 0.8 for every passing state, which is the range
the plan's threshold was set for. Check (c), anatomy, is open: it needs
the H&E or DAPI image next to the state map, judged by a reader, and is
possible where coordinates exist (the Heming and Yu Xenium cells, the
Shiers spots, the DPN nerve spots); Zhang, Vroman and Tavares-Ferreira
have no spot positions in the public files, so (c) cannot be done for
them.

## 3.6 Anatomy check inputs

The images for check (c) are fetched after the mapping, once it was
decided that the check would be done. What each source holds:

- Tavares-Ferreira 2022 (Box): the H and E images are 16 TIFF files of
  1.65 to 1.70 GB (24.7 GiB), one per section, in a folder next to the
  H5 and Loupe folders; the listing script now has an `--images` mode
  for them. The spot positions of these sections are still only in the
  Loupe files, so the images cannot be overlaid with the mapping; they
  are kept so that the sections can be looked at, and a panel is drawn
  only if positions turn up.
- Shiers 2025: the SPARC portal holds a copy of the dataset (discover
  dataset 490, version 1, DOI 10.26275/vz0c-uet7, CC BY) with the full
  Space Ranger spatial folder per section: high-resolution and
  low-resolution tissue images, scale factors, positions, the aligned
  fiducials and the Loupe files. The high-resolution image (4 to 7 MB)
  and the scale factors of each section are taken from the dataset's
  public S3 bucket, renamed with the config section name (the SPARC
  sample folders are the GEO sample names; donor DN0003's two sections
  are M2-1 and M2-2 there as on GEO). The positions already in use
  come from GEO and are the same files. A CFDE record describing
  "16 dorsal root ganglia from 8 patients (left and right)" looked like
  the Tavares-Ferreira data but is SPARC dataset 476, cervical DRG from
  fusion surgeries, a different study; Tavares-Ferreira is not on SPARC.
- Heming 2025: the DAPI morphology OME-TIFF of each section is in the
  GEO tar (0.6 to 2.8 GB gzipped per section, 10.5 GB in all); it is
  pyramidal and the cell centroids are in microns, so a low pyramid
  level with the pixel size from the OME metadata is enough for the
  panels.
- DPN sural nerve: the Seurat object carries the low-resolution tissue
  image of each of its six sections, which the R export now writes as
  PNG next to the coordinates (whose scale factors it already had); the
  Box folder has no other image.

The `tar -x ... -C` form used first put the image files next to the
tar instead of in the extracted directory; they were decompressed there
and moved. Nothing of the mapping results depends on the images.

All 49 image files (16 TIFF, 16 PNG, 17 JSON; DN0134-S1-1 has two scale
factor files from two Space Ranger runs, and the one whose scaled spot
positions fit inside the image is used) verify against the manifest.
The Xenium morphology file is a DAPI z-stack of 11 planes 3 microns
apart, JPEG 2000 compressed, whose first plane is nearly empty; the
panels use the maximum projection at pyramid level 3 (1.7 microns per
pixel). The server environment needed `imagecodecs` for the JPEG 2000
tiles, now in `environment.yml`.

Panels. For the Visium sections a dominant-state map would be empty
for the three nerve states in a ganglion (they are never the most
numerous cells of a 30 cell spot there, and in the sural nerve the
perineurial state dominates only part of the perineurium), so each
Visium panel shows the estimated number of cells of the state in every
spot instead, on the H and E; for the Xenium sections the cells
assigned to the three states are drawn on the DAPI image. The panels
are `results/spatial/anatomy/<dataset>/<section>.jpg` (30 sections: 8
Heming, 16 Shiers, 6 DPN nerve) and the reader's sheet is
`results/spatial/anatomy/anatomy_check.tsv`, one row per section and
state, to be filled with yes, no or not present. The judgement is the
reader's and is not recorded here until it is made. What the panels
show at first sight, as a note for the reader: in the Heming sural
nerve sections the perineurial cells form rings around the fascicles
(S01 most clearly, with 1,240 perineurial cells), the endoneurial
fibroblasts sit inside the fascicles and the epineurial fibroblasts
between them; in the DPN sural nerve sections the perineurial abundance
is highest in bands across the fascicles rather than only at their
rims, which the low-resolution image (600 pixels across) makes hard to
place; in the DRG sections the three states together hold 1 to 4
percent of the estimated cells, spread thinly over the ganglion with no
concentration at the capsule, so the expected answer there is "not
present" rather than a yes or a no.

## 3.7 Yu 2024 full Xenium bundles

It was decided to fetch the Dryad copy of the Yu 2024 Xenium runs (doi
10.5061/dryad.gf1vhhmxq, CC0, 12.4 GB in four `.rar` bundles and a
README). The GEO series gives only the per-region transcript tables and
the 1,415 hand-segmented neurons; the Dryad README says each bundle is
the complete Xenium onboard-analysis output of one section (cells with
segmentation, cell by gene matrix, transcripts, DAPI morphology images)
plus the post-Xenium H and E image of the section (`DRG_area_<n>.ome`)
and its alignment file to the morphology image. That gives, for this
dataset, both the non-neuronal cells that the hand-segmented table
leaves out and an H and E image for check (c). The README also states
which donor and level each region is: Region 1 donor N2 lumbar 2,
Region 2 donor N4 lumbar 3, Regions 3 and 4 donor N4 thoracic 12.

Dryad still refuses scripted downloads (API 401 without a token, web
link 403), so the files were saved from the dataset page in the
laptop's browser and verified afterwards against the sha-256 digests
the Dryad API publishes, with a small script that also writes the
manifest rows; the fetch is recorded as manual with the click sequence
in the runbook. The laptop system disk had only 14 GB free, so the
bundles were taken one at a time and moved to the repo data directory
on the external drive between downloads. The server was not reachable
over ssh while the downloads ran (connection refused), so the copy to
the server and anything done with the bundles there wait for it to come
back. Nothing in the mapping results so far uses these files.

The Region 1 bundle (4.9 GB, 6.2 GB unpacked, 153 entries) holds the
Xenium output folder as the README describes it: `cells.csv.gz` and the
cell by gene matrix (`cell_feature_matrix.h5` and the MEX folder), the
cell and nucleus boundaries, the transcripts in four formats
(`transcripts.csv` alone is 531 MB), the morphology images
(`morphology.ome.tif` 2.8 GB, the focus and maximum-projection images
0.3 GB each, and an unpacked `morphology_mip_extracted.tif` of 1.9 GB
that the README does not mention), the onboard clustering, and the
post-Xenium H and E image `DRG_area_1.ome.tif` (279 MB) with its
alignment to the morphology image (`keypoints.csv` and a 3 by 3
`matrix.csv`). `bsdtar` on the laptop lists the archive, so the bundles
can be unpacked without a separate rar tool.
All five files verified (sizes and sha-256 as published by Dryad; md5
and the web download link of each file in the manifest). The copy to
the server is still to be made.
The copy to the server went over a VPN at 4 to 6 MB/s and the
connection dropped once near the end of the last bundle; `rsync
--partial` resumed it and all five files verified on the server.

## 3.8 Where the project data lives on the server

The server's system disk holds every home directory and was 91 percent
full, 85 GB of it the project's `data/raw` and `data/interim` in the
checkout. The server has two arrays besides: the sequencer's output
array (7.3 TB, 1.7 TB free), which holds the lab's long-read runs and
is left alone, and a 22 TB data array with 8.3 TB free where the lab
keeps its analysis data and references. It was decided to move the two
large project directories to a project directory on the data array
and leave symbolic links in the checkout, so that no path in the code,
the config or this runbook changes. The array is a stripe without
redundancy, which is acceptable because every raw file is listed with
its source and checksum in the manifest and every intermediate file is
rebuilt by the workflow. The move script copies, compares the copy by
checksum, and removes the source only after the comparison is clean.
The move left both comparisons clean (83 GB raw, 13 GB interim on the
array), the system disk went from 85 to 169 GB free, and the Dryad
files verify through the links. The checkout's manifest and the
small metadata stay on the system disk with the code.

## 3.9 Yu 2024: every segmented cell, and H and E panels

It was decided to use the Dryad bundles for two things: the transfer of
every segmented cell of the four Yu sections (the GEO table holds only
the 1,415 hand-segmented neurons, so the non-neuronal states of the
ganglion were unchecked in this dataset) and the anatomy panels of check
(c), which the post-Xenium H and E image and its alignment make
possible. The all-cell data are a second dataset id next to the neuron
table rather than a replacement: the table is the authors' curated
neuron set, its checks are recorded, and the two answer different
questions (decisions).

Nine members are taken out of each bundle (cells table, cell by gene
matrix, the H and E image, its alignment files, run metadata), 0.1 to
0.3 GB per region, with libarchive through its Python binding, since the
server has no rar tool; on the laptop the same nine members taken out
with bsdtar have identical checksums.

Two things were learnt from the bundles about data already in use. The
GEO neuron table's cell ids (such as "3312-2030") had been read as x-y
image pixels. Against the Dryad cells table they are row-column of the
morphology image, as ImageJ names its ROIs: swapped and scaled by 0.2125
microns per pixel, the neurons fall a median 5.4 microns (90 percent
within 12) from a segmented cell centroid, and the cells they fall on
have a median area of 187 square microns against 126 for all cells, as
neuron somata should; in the original order the nearest cell is over a
millimetre away. The config now carries the axis order and the pixel
size, and the neuron sections have coordinates in microns. Nothing in
the transfer depends on coordinates, so the recorded neuron checks are
unchanged. Second, the alignment matrix exported by Xenium Explorer
(3 by 3) reproduces the keypoints when applied to the H and E side (for
Region 2, median residual 7.6 morphology pixels, 1.6 microns, largest
17), so it maps H and E pixels to morphology pixels, and its inverse
places the cells (microns over 0.2125) on the H and E, whose pixel is
0.5567 microns (the matrix scale is 2.62). The panels read the H and E
at pyramid level 1 (1.1 microns per pixel); the grey cells cover the
tissue outline in all four sections, so the placement is right.

Candidate states. The compartment rule alone would offer only the
ganglion states, among which none of the three states check (c) asks
about exists (the ganglion datasets label their sheath cells as
fibroblast), so the Yu panels would be empty by construction. A DRG
section carries its capsule and the root entering it, which are nerve
sheath tissue, so for this dataset the three nerve sheath states are
added to the ganglion candidates (24 states in all, including the
sympathetic and injured neuron states the ganglion compartment
contributes). The panel is the authors' 100 gene DRG panel, 98 of whose
genes are in the reference signatures; it was designed for neuron
subtypes, so the neuron families have their full 20 markers within the
panel but the non-neuronal states only 4 (macrophage, both Schwann
states) to 13 (satellite glia), and check (a) for them reads the
panel's coverage as much as the mapping.

Outcome. 151,974 cells in the four sections, 98 shared genes, 24
candidate states. Between 10 and 19 percent of the cells per section
are below the count floor (2,495 to 7,947 cells with fewer than 20
transcripts over the shared genes) and 106 to 232 are ambiguous; the
rest are assigned with a median posterior of 1.0. The assigned cells
are, in order, fibroblast (22,063 over the four sections), sympathetic
neuron (15,593), non-myelinating Schwann (14,008), myelinating Schwann
(13,075), endothelial (12,714), satellite glia (12,440), macrophage
(11,254), the ten DRG neuron families (20,167 together, A-PEP 9,086
and C-NP 4,660 the largest), injured neuron (1,889), mural (1,317),
epineurial fibroblast (779), adipocyte (250), endoneurial fibroblast
(232) and perineurial (230). The neuron count is far above the 1,415
hand-segmented somata, as expected from a nucleus-based segmentation
that cuts a large soma into several cells and lets the satellite glia
around it pick up neuron transcripts. Checks: every neuron family and
the injured neuron state pass (a) and (b) in every section (AUROC 0.85
to 0.99, Spearman 0.64 to 0.84); endothelial, fibroblast, macrophage,
both Schwann states, satellite glia, mural, adipocyte and endoneurial
fibroblast pass both in every section where they reach 20 cells
(median AUROC 0.86 to 0.96, Spearman 0.58 to 0.72). Epineurial
fibroblast fails (a) in three of four sections (median AUROC 0.67 on
six markers) and with the Visium section already counted is now not
mapped in the ganglion on five sections; perineurial fails (a) once and
(b) once in four (median AUROC 0.83, Spearman 0.53) and is mapped.
Sympathetic neuron takes 15,593 cells, fails (a) in two of four
sections and (b) in one (median AUROC 0.71, Spearman 0.51), and by the
plan's rule (more than half must fail) counts as mapped; a DRG holds no
sympathetic neurons, so this state is absorbing cells whose panel
counts are thin and neuron-like, most likely satellite glia and glial
cells with neuron transcripts, and whether to drop it from the
candidates of a DRG section is a decision for the user, not made here.
The status table gains the ganglion rows of A-LTMR, Ab-LTMR, C-LTMR,
adipocyte, endoneurial fibroblast, injured neuron, perineurial and
sympathetic neuron, the other ganglion rows count four sections more,
and no earlier mapped or not mapped verdict changes. The transfer's
method label in the collected checks now reads "signature transfer"
instead of the stale "scANVI transfer".

The first run of the script rebuilt the section summary from the
per-dataset section files, of which the server only held the two Yu
ones, so the summary and the collected checks were briefly reduced to
the Yu datasets; the summary was restored from the repository with the
four all-cell rows added, the collection was rerun, and the script now
replaces the Yu rows in the existing summary instead.

Panels. Four JPEGs (3,200 pixels wide) with the H and E on the left
and the assigned cells on the right, the three sheath states in colour
and the rest grey, and twelve rows in the reader's sheet (perineurial
97, 32, 51 and 50 cells in Regions 1 to 4; epineurial fibroblast 154,
96, 321 and 208; endoneurial fibroblast 34, 11, 124 and 74). The
coloured cells are a few hundred among 25,000 to 37,000 assigned, so
they are picked out only on zooming; the H and E shows the ganglion
with its capsule and the root fascicles attached, and the Region 2
perineurial cells lie along the capsule edge at first sight. The
judgement is the reader's and is not recorded here until made.
