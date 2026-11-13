#!/usr/bin/env bash
# Stage 6 (runbook 6.2): the pnatlas gene features, then PNPrank's temporal
# back-test with and without them. Prints "FEATURES <exit>" and
# "BACKTEST <exit>" and ends with BACKTESTDONE. Runs in the project checkout;
# the PNPrank checkout at its submission tag is the first argument. Extra
# arguments go to pnprank_backtest.py (probe: --probe).
# Usage: workflow/scripts/run_backtest.sh <pnprank checkout> <data_root> [results dir] [n_jobs] [backtest options...]
set -uo pipefail
repo=${1:?pnprank checkout}
root=${2:?data root}
res=${3:-results/backtest}
jobs=${4:-1}
shift 4 2>/dev/null || shift $#
mkdir -p logs "$res" "$root/processed"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pnprank-bt
nice -n 15 python workflow/scripts/pnatlas_features.py --specificity "$root/interim/genetics/specificity.tsv.gz" \
  --expression "$root/interim/genetics/expression_cpm.tsv.gz" --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz \
  --status results/spatial/state_mapping_status.tsv --gtex-counts "$root/interim/deconv/gtex_counts.tsv.gz" \
  --fractions results/deconvolution/gtex_bayesprism.tsv --coloc results/genetics/coloc_all_pairs.tsv \
  --universe "$root/external/pnprank/features_static.tsv.gz" --out "$root/processed/pnatlas_gene_features.tsv.gz" \
  > logs/pnatlas_features.log 2>&1
echo "FEATURES $?"
nice -n 15 python workflow/scripts/pnprank_backtest.py --pnprank-repo "$repo" \
  --features-static "$root/external/pnprank/features_static.tsv.gz" \
  --pnatlas-features "$root/processed/pnatlas_gene_features.tsv.gz" --out "$res" --n-jobs "$jobs" "$@" \
  > logs/pnprank_backtest.log 2>&1
echo "BACKTEST $?"
echo BACKTESTDONE
