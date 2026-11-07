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
