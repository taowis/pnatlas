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
