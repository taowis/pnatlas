# Runbook, stage 2: joint ganglia and nerve reference

Commands in run order. `$REPO` and `$DATA_ROOT` are set as in
[README.md](README.md). Every step is run with the `pnatlas` environment
active unless stated otherwise.

### 2.1 Environment and machine config

The environment build is in the Environment section of [README.md](README.md).
Machine config for this run:

```bash
cd "$REPO"
cp config/paths.example.yaml config/paths.yaml
# data_root: a directory on the external drive with at least 300 GB free
# pnprank_root: the pnprank-dev clone (tag submitted-human-genetics)
mkdir -p "$DATA_ROOT"
conda activate pnatlas
python -c "import torch; print('mps', torch.backends.mps.is_available())"
```

Prints `mps True` on the laptop.
