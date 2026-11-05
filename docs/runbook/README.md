# Runbook

Every command and manual step used in the analysis, in run order, so the whole
analysis can be redone by hand. Each stage has its own file. Blocks give the
date, the repo commit before the step, the command, the output and a check.

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

The conda environment (`conda env create -f environment.yml`) is not needed for
stage 1; the stage 1 scripts use the Python standard library only and are run
with the system `python3`. The environment build will be recorded here when
stage 2 first needs it.

## Stages

| File | Stage |
|---|---|
| [01_datasets.md](01_datasets.md) | 1. Data inventory and access (logs in `logs/01_datasets/`) |
