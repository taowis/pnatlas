#!/usr/bin/env bash
# Stage 5 (runbook 5.2 and 5.3): gene class enrichment on the reference
# specificity, then the colocalisation of every GWAS of the plan with the
# tibial nerve eQTL, and the collection of the results. Each step prints a
# "<STEP> <name> <exit code>" line; the run ends with GENETICSDONE. Runs in
# the project checkout; Python steps use the pnatlas environment, R steps
# pnatlas-r. Extra arguments go to coloc_prepare.py (probe: --max-regions N).
# Usage: workflow/scripts/run_genetics.sh <data_root> [results dir] [prepare options...]
set -uo pipefail
root=${1:?data root}
res=${2:-results/genetics}
shift; shift || true
gen=$root/interim/genetics
eq=$root/raw/eqtl_catalogue_gtex_tibial
mkdir -p logs "$res" "$gen/coloc"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pnatlas
nice -n 15 python workflow/scripts/gene_class_enrichment.py --specificity "$gen/specificity.tsv.gz" \
  --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz --truth data/metadata/pnprank_truth_set.tsv \
  --panels-dir data/metadata --results-dir "$res" > logs/gene_class_enrichment.log 2>&1
echo "ENRICHMENT $?"
p1=$(python -c "import yaml; c=yaml.safe_load(open('config/config.yaml'))['genetics']['coloc']; print(c['p1'], c['p2'], c['p12'])")
read -r P1 P2 P12 <<<"$p1"
for g in gwas_diabetic_neuropathy:GCST90129436 gwas_mvp_t2d_neurological:GCST90475676 gwas_mvp_hereditary_idiopathic_pn:GCST90475844 \
         gwas_neuropathic_pain:GCST90859373 gwas_carpal_tunnel:GCST90018813; do
  id=${g%%:*}; acc=${g##*:}
  conda activate pnatlas
  nice -n 15 python workflow/scripts/coloc_prepare.py --gwas-id "$id" --gwas "$root/raw/$id/$acc.h.tsv.gz" \
    --meta "$root/raw/$id/$acc.h.tsv.gz-meta.yaml" --eqtl "$eq/QTD000286.all.tsv.gz" \
    --credible-sets "$eq/QTD000286.credible_sets.tsv.gz" --lbf "$eq/QTD000286.lbf_variable.txt.gz" \
    --out-dir "$gen/coloc" "$@" > "logs/coloc_prepare_$id.log" 2>&1
  echo "PREPARE $id $?"
  conda activate pnatlas-r
  nice -n 15 Rscript workflow/scripts/coloc_run.R --dir "$gen/coloc/$id" --out "$gen/coloc/$id/coloc.tsv" \
    --p1 "$P1" --p2 "$P2" --p12 "$P12" > "logs/coloc_run_$id.log" 2>&1
  echo "COLOC $id $?"
done
conda activate pnatlas
python workflow/scripts/coloc_collect.py --coloc-dir "$gen/coloc" --specificity "$gen/specificity.tsv.gz" \
  --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz --status results/spatial/state_mapping_status.tsv \
  --results-dir "$res" > logs/coloc_collect.log 2>&1
echo "COLLECT $?"
echo GENETICSDONE
