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
