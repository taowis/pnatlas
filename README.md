# pnatlas

**Status:** Stage 3 (spatial mapping) done; stage 4 (GTEx tibial nerve deconvolution) next.

| Stage | Description | Status |
|---|---|---|
| 0 | Scaffold | done |
| 1 | Data inventory and access; analysis plan frozen | done (dbGaP raw-read request pending, not on the critical path) |
| 2 | Joint ganglia + nerve single-cell reference | done: reference integrated and annotated, state rule checked |
| 3 | Spatial mapping (Xenium, Visium) | done: 44 Visium sections and 3 Xenium panels mapped, checks (a), (b) and (c) collected; status per state in `results/spatial/state_mapping_status.tsv` |
| 4 | GTEx tibial nerve deconvolution | not started |
| 5 | Genetics: neuropathy gene classes, nerve eQTL colocalisation | not started |
| 6 | PNPrank v2: temporal back-test with the new features | not started |
| 7 | Static portal, Zenodo release, public repository | not started |

## What this is

A reference of human peripheral nervous system cell states that joins sensory
and autonomic ganglia with peripheral nerve, places those states in tissue
with spatial data, estimates them in bulk GTEx tibial nerve, and links them to
inherited and acquired neuropathy genetics. Only public data are used.

Harmonised DRG atlases already exist, so the reference on its own is not the
contribution. The questions are:

1. Which ganglion and nerve cell states carry the genes for each class of
   inherited neuropathy (CMT, dHMN, HSAN)?
2. How does the cell composition of tibial nerve vary across GTEx donors, and
   which nerve eQTLs act in specific cell states?
3. Do nerve and spatial features improve the PNPrank temporal back-test over
   the published PNPrank features? This comparison is fixed in
   `docs/analysis_plan.md` before it is run, and a null result will be
   reported as such.

Prior work and the remaining gap are in [`docs/novelty_check.md`](docs/novelty_check.md).

## Setup

```bash
conda env create -f environment.yml
conda activate pnatlas
pip install -e . --no-deps
cp config/paths.example.yaml config/paths.yaml   # then edit for this machine
pytest
```

## Layout

```
config/      parameters, dataset list, cell-type map; paths.yaml is local only
workflow/    Snakemake rules, one file per stage
src/pnatlas/ Python package and CLI
tests/
data/        see data/README.md; only README and small metadata are committed
results/     tables and figures; snapshot commit in results/SNAPSHOT.txt
docs/        novelty check, analysis plan, decisions, stage reports
portal/      static gene lookup site
manuscript/
```

## Related

- PNPrank (`taowis/pnprank-dev`): gene ranking for inherited peripheral
  neuropathy from single-cell atlases. Its processed atlases and truth set are
  reused here.

## Licence

MIT for code. Data remain under the terms of their original providers.
