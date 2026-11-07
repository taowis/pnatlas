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
