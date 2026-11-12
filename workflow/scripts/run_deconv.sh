#!/usr/bin/env bash
# Stage 4 (runbook 4.2): prepare the nerve reference, the validation
# mixtures and the GTEx tibial nerve counts, deconvolve both bulk sets with
# BayesPrism and MuSiC, then evaluate (validation (a) and (b), covariate
# models). Each step prints a "<STEP> <name> <exit code>" line; the run ends
# with DECONVDONE. Runs in the project checkout; the Python steps use the
# pnatlas environment, the R steps pnatlas-r (both activated here).
# Usage: workflow/scripts/run_deconv.sh <data_root> [results dir] [prepare options...]
set -uo pipefail
root=${1:?data root}
res=${2:-results/deconvolution}
shift; shift || true
prep=$root/interim/deconv
gtex=$root/raw/gtex_v10_tibial
mkdir -p logs "$res"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pnatlas
nice -n 15 python workflow/scripts/deconv_prepare.py --counts "$root/interim/reference_merged.h5ad" \
  --annotation "$root/interim/reference_annotation.tsv.gz" --gtex-counts "$gtex/gene_reads_v10_nerve_tibial.gct.gz" \
  --out-dir "$prep" "$@" > logs/deconv_prepare.log 2>&1
echo "PREPARE $?"
cores=$(python -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['deconvolution']['bayesprism_cores'])")
conda activate pnatlas-r
for m in bayesprism music; do
  for b in mixtures gtex; do
    nice -n 15 Rscript workflow/scripts/deconv_run.R --prep-dir "$prep" --bulk "$prep/${b}_counts.tsv.gz" \
      --method $m --out "$res/$([ $b = mixtures ] && echo validation || echo gtex)_$m.tsv" --cores "$cores" \
      > "logs/deconv_${m}_${b}.log" 2>&1
    echo "RUN $m $b $?"
  done
done
conda activate pnatlas
python workflow/scripts/deconv_evaluate.py --prep-dir "$prep" --results-dir "$res" \
  --sample-attributes "$gtex/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt" \
  --subject-phenotypes "$gtex/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt" > logs/deconv_evaluate.log 2>&1
echo "EVALUATE $?"
echo DECONVDONE
