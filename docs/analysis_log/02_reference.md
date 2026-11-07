# Analysis log, stage 2: joint ganglia and nerve reference

What was found and decided at each step. Commands are in
`docs/runbook/02_reference.md` and, for the environment, in
`docs/runbook/README.md`.

## 2.1 Environment

The conda environment was built for the first time here. `environment.yml`
was extended with the stage 2 tools before the build: scvi-tools (with the
CPU build of PyTorch from conda-forge), jax and jaxlib for scib-metrics,
leidenalg and igraph for clustering, h5py, matplotlib and seaborn. Three
packages are not on conda-forge for this platform and are pinned through pip
inside the environment: harmonypy, scib-metrics and scrublet.

The solve picked scvi-tools 1.4.2 rather than the latest 1.5.2 because the
file pins Python 3.11; the difference does not matter for scVI and scANVI as
used in the plan, so the pin stays. The environment takes 2.3 GB on the
laptop's internal disk. PyTorch sees the Apple GPU (MPS), which is enough for
model development on subsets; full-reference scVI runs go to the server. The
package tests pass in the new environment and the `pnatlas` entry point
resolves. The solved versions are in `docs/runbook/environment_lock.yml`,
exported without the machine-specific prefix line.

A machine config `config/paths.yaml` (not committed) was created with the
data root on the external drive, where 485 GB are free, and the disk budget of
300 GB from the example file.

## 2.2 Reference layer download

Where to compute was settled first (decision in `docs/decisions.md`, stage 2):
integration on the server, Gadi for the CPU-parallel stages later, the laptop
for development. The server had an idle GPU, 183 GB free and only the
sequencer's basecall service running when the download was started.

The download list was built by a script from the live GEO directory listings
and the Pennsieve file browser rather than typed, so the file names and sizes
are what the sources serve. For stage 2 only the human count matrices and the
authors' metadata are taken: 19 files, 11.3 GB. The SPARC 465 ATAC fragments
(17.8 GB) and the mouse files in the trigeminal series are left out; the
cross-species Jung 2023 parent series GSE201654 is not fetched because the
human sub-series GSE201586 holds the same human samples plus the annotated
object. The macaque anchor GSE165569 has no supplementary files on GEO; its
counts are per sample and are fetched when the anchor is used. All downloads
go through `download_files.py`, which records md5 and size in
`data/MANIFEST.tsv` as each file completes and can be rerun after an
interruption.

## 2.3 Reading the datasets

Each dataset is read by one script from a per-dataset block in
`config/reference_inputs.yaml` (format, files, which author columns give
sample, donor and labels), so the decisions about donors and labels are in a
config file rather than spread over notebooks. Three sources are R objects
(the painseq Seurat object, the Jung 2023 annotated object, the trigeminal
raw count RDS); a small R script exports them as matrix market plus CSV and
the Python side never loads R objects directly.

What the raw files turned out to contain:

- Nguyen 2021: six CSV count tables, one per preparation and donor, 1,837
  nuclei in all, no labels in the GEO files.
- Jung 2023 human sub-series: 18 Cell Ranger h5 files (filtered and raw per
  sample) from 7 donors; donor and labels come from the authors' object.
- Yu 2024: 1,136 single somata from 3 donors; the cell names encode donor and
  ganglion level.
- Zhang (GSE189501): 30,701 nuclei from 3 DRG runs.
- SPARC 465: the Cell Ranger ARC filtered matrix holds 36,601 genes and
  139,026 ATAC peaks; the peaks are dropped. The SDS samples sheet shows that
  library lib-1 pools three donors, so 8 libraries give 6 donor-level units
  (one of them a pool) rather than 7 donors; the pool is kept as its own unit
  and is not counted as three donors in any donor-level test.
- Heming 2025: 37 cellbender-filtered h5 files, one per donor, 365,708
  annotated nuclei after keeping the cells in the authors' metadata; two
  label levels, the Milbrandt sciatic nerve label and the finer cluster.

Two problems in the first pass. The server sync script only copied named
config files, so the new inputs file was missing on the server; it now
copies the whole config directory except the machine file. The Heming
build took 33 GB because all 37 samples were concatenated before the
unannotated cells were removed and the memory guard logged a warning at 19
percent free; the builder now subsets each sample to the annotated cells as
it is read. A metadata column with mixed text and numbers also failed to
write to h5ad; object columns are now written as text.

The R export failed at first for two reasons worth recording. The two GEO
`.rds.gz` and `.RDS.gz` files are gzip archives of files that are themselves
gzip-compressed RDS, so `readRDS` on the path stops with "unknown input
format"; the script now opens a second decompression layer when the first
read fails. The painseq object is a Seurat v4-style assay, and SeuratObject
5.5 has removed the old `slot` argument, so the counts are taken with
`layer = "counts"`. The R environment itself (R 4.4.3, Seurat 5.5.1) could
only be built once the Python environment build had finished, because two
conda solves at once lock the shared repodata cache.
