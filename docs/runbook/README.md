# Runbook

Every command and manual step used in the analysis, in run order, so the whole
analysis can be redone by hand. Each stage has its own file. Each step gives
the command, what it reads and writes, and a check.

## Setup

```bash
git clone https://github.com/taowis/pnatlas-dev.git
cd pnatlas-dev
export REPO="$PWD"

# Machine paths (large data never goes in git)
cp config/paths.example.yaml config/paths.yaml
# edit config/paths.yaml: data_root, pnprank_root
export DATA_ROOT=/path/to/pnatlas-data    # same value as data_root
```

## Environment

Built once at the start of stage 2 (stage 1 used only the system `python3`).
`environment.yml` lists the packages: scanpy and scvi-tools (PyTorch, CPU
build) for integration, jax and scib-metrics for the integration benchmark,
harmonypy for the comparison method, scrublet for doublets, snakemake for
the workflow. Three packages are not on conda-forge and are installed by pip
inside the environment, as listed in the `pip:` section of the file.

```bash
cd "$REPO"
conda env create -f environment.yml -n pnatlas -y
conda activate pnatlas
pip install -e .                      # the pnatlas package in src/, editable
python -c "import scanpy, scvi, harmonypy, scib_metrics, scrublet; print(scanpy.__version__, scvi.__version__)"
pytest -q
conda env export --no-builds > docs/runbook/environment_lock.yml
```

The solved versions are in `docs/runbook/environment_lock.yml`, written by the
last command; to rebuild the exact environment use that file instead of
`environment.yml`. The build was done on an Apple silicon laptop (osx-arm64,
16 GB); on Linux or a GPU server the same `environment.yml` is used and the
lock file will differ in platform-specific builds. Memory-heavy steps (scVI on
the full reference, cell2location) are run on the server, not on the laptop.

## Stages

| File | Stage |
|---|---|
| [01_datasets.md](01_datasets.md) | 1. Data inventory and access |
| [02_reference.md](02_reference.md) | 2. Joint ganglia and nerve reference |
