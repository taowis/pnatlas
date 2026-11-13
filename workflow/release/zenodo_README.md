# pnatlas data archive

Processed data and result tables of pnatlas, a multiscale reference of human
peripheral nervous system cell states built from public data: single-nucleus
and single-cell ganglia and nerve datasets joined into one reference, placed
in tissue with Xenium and Visium sections, estimated in GTEx tibial nerve, and
linked to inherited and acquired neuropathy genetics. The code, configuration,
runbook and analysis log are in the repository named in `CITATION.cff` there;
`data_README.md` and `MANIFEST.tsv` here record every public file the analysis
downloaded, with URL, version, access date and checksum, so that the raw
layers can be fetched again.

## Contents

Zenodo takes at most 100 files per record, so the archive is uploaded as
the reference object on its own, one `tar.gz` per component
(`reference_tables`, `spatial`, `deconvolution`, `results`, each unpacking
into the directory of that name) and the four record files (`README.md`,
`data_README.md`, `MANIFEST.tsv`, `MD5SUMS`). `MD5SUMS` lists every file
inside the bundles by its unpacked path.

- `pnatlas_reference_v0.1.0.h5ad`: the reference object, raw counts
  after quality control for 597,455 cells by 19,420 genes, with per cell the
  dataset, sample, donor, batch, tissue, compartment, protocol, the authors'
  label, the final class and state, how the state was assigned
  (`label_source`) and its transfer confidence, and the scANVI and scVI
  embeddings. The cells of the humanDRG portal release (Bhuiyan 2025) are not
  included, because that portal does not allow cell-level derivatives to be
  redistributed; they are 47,226 neurons of the full reference, and the
  runbook rebuilds the full object from the portal download.
- `reference/pnatlas_reference_cells.tsv.gz` (in `reference_tables.tar.gz`): the cell table of the object.
- `reference/state_specificity.tsv.gz`, `state_expression_cpm.tsv.gz`,
  `state_pseudobulk_cells.tsv`: per gene, the specificity and log CPM in each
  of the 31 states, computed over every dataset, and the cells behind each
  pseudobulk.
- `spatial/cell2location/<dataset>/<section>/`: cell2location abundances per
  spot (`abundance.tsv.gz`), the spot table, the training history and the
  checks of each Visium section; `spatial/panel_transfer/`: the state
  assigned to each cell of the Xenium sections; `spatial/signatures/`: the
  cells and the reference signatures used by cell2location.
- `deconvolution/`: the single-cell reference (counts, cells, genes, donors)
  used for GTEx deconvolution and the held-out pseudobulk mixtures with their
  true fractions.
- `results/`: every result table of the repository at the snapshot commit
  (reference, spatial, deconvolution, genetics, back-test).
- `MD5SUMS`: checksums of every file in the bundles, by unpacked path.

## Terms

Code is MIT licensed. The data derive from public datasets that remain under
the terms of their providers (GEO, SPARC CC BY 4.0, Dryad CC0, GTEx, GWAS
Catalog, eQTL Catalogue); the source and terms of each are listed in
`data_README.md`.
