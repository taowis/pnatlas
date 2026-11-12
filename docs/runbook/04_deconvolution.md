# Stage 4: GTEx tibial nerve deconvolution

Plan section 5. The nerve part of the stage 2 reference (the Heming 2025
sural nerve nuclei, 37 donors, 16 states) is used to estimate the cell
state composition of the GTEx v10 tibial nerve samples with BayesPrism
(primary) and MuSiC (secondary), after a validation on pseudo-bulk
mixtures of held-out donors and against the GTEx pathology notes. All
steps run on the server; the laptop holds the code and the result
tables. `$REMOTE` is `user@host` of the server.

## 4.1 Inputs and environment

The five GTEx files (per-tissue read counts and TPM of "Nerve - Tibial",
sample attributes with their data dictionary, subject phenotypes) are
listed in `data/metadata/stage4_files.tsv` with their sizes from HEAD
requests to the bucket and fetched with the stage download script, which
records them in the manifest:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && setsid -f nohup python3 -I workflow/scripts/download_files.py \
  data/metadata/stage4_files.tsv --data-root data --manifest data/MANIFEST.tsv \
  --fetch "workflow/scripts/download_files.py data/metadata/stage4_files.tsv" > logs/stage4_download.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && grep -c ^done logs/stage4_download.log'   # 5
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/MANIFEST.tsv" data/MANIFEST.tsv
```

The R environment gains the dependencies of the two methods
(`environment-r.yml`; on the existing server environment the same
packages are added with conda) and the methods themselves from GitHub at
fixed commits, as neither repository has release tags
(`workflow/scripts/deconv_install.R`: BayesPrism 2.2.3, MuSiC 1.0.0):

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && cd ~/Projects/pnatlas-dev && \
  nice -n 15 conda install -y -n pnatlas-r -c conda-forge -c bioconda r-remotes r-nmf r-snowfall r-gplots r-nnls \
    r-mcmcpack r-ggplot2 bioconductor-biobase bioconductor-singlecellexperiment bioconductor-scran \
    bioconductor-biocparallel bioconductor-toast > logs/stage4_r_install.log 2>&1 && \
  conda activate pnatlas-r && nice -n 15 Rscript workflow/scripts/deconv_install.R > logs/deconv_install.log 2>&1; \
  tail -2 logs/deconv_install.log'
```

Check: the last two lines of `logs/deconv_install.log` give the versions
of BayesPrism and MuSiC.

## 4.2 Reference, mixtures, deconvolution and evaluation

Parameters are the `deconvolution` block of `config/config.yaml`: nerve
compartment, states with at least 200 nerve cells, 7 of the 37 donors
held out (chosen by the project seed), up to 1,000 reference cells per
state, 200 mixtures of 500 cells with flat Dirichlet proportions, the
pass thresholds of the plan (mean Pearson above 0.8, RMSE below 0.1),
the GTEx tissue label and the cores for BayesPrism. A probe with ten
mixtures, 100 cells per state and 20 GTEx samples checks the chain
first; then the full run, both under the memory guard:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader && who && free -g | sed -n 2p && \
  ~/bin/memguard.sh start logs/memguard.log "deconv_prepare|deconv_run|deconv_evaluate" && \
  setsid -f nohup workflow/scripts/run_deconv.sh data results/deconvolution_probe --n-mixtures 10 --max-cells-per-state 100 --gtex-samples 20 \
    > logs/deconv_probe.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/deconv_probe.log'    # PREPARE 0, four RUN ... 0, EVALUATE 0, DECONVDONE
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && rm -rf data/interim/deconv results/deconvolution_probe && \
  setsid -f nohup workflow/scripts/run_deconv.sh data results/deconvolution > logs/deconv.log 2>&1 < /dev/null'
```

The probe writes its interim files to the same `data/interim/deconv`
directory, which is why it is removed before the full run. When
`logs/deconv.log` ends with DECONVDONE, stop the guard and bring the
results back:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/deconv.log && ~/bin/memguard.sh stop logs/memguard.log; grep -E "WARN|PAUSE" logs/memguard.log | tail -3; \
  grep -E "^reference|^nerve|^genes|^mixtures" logs/deconv_prepare.log'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/deconvolution/" results/deconvolution/
rsync -a --include "*.tsv" --exclude "*" "$REMOTE:~/Projects/pnatlas-dev/data/interim/deconv/" "$DATA_ROOT/interim/deconv/"
cp "$DATA_ROOT/interim/deconv/prepare_summary.tsv" "$DATA_ROOT/interim/deconv/donors.tsv" results/deconvolution/
```

Check: `results/deconvolution/validation_pseudobulk.tsv` has one row per
method with the pass column; `method_choice.tsv` names the primary
method; `gtex_proportions.tsv` has one row per tibial nerve sample with
the covariates and the state proportions; `validation_histology.tsv`
and `gtex_covariate_models.tsv` hold validation (b) and the models.

## 4.3 Evaluation rerun with the marker concordance

The evaluation step was rerun on the laptop after the full run, with
three additions to `deconv_evaluate.py`: the Benjamini-Hochberg step
skips states without a p value (a state estimated at zero in every
sample has no variance; before the fix every q value of the term came
out missing), validation (b) reads the fat percentage written in the
pathology notes instead of a yes or no mention and runs for both
methods, and a marker concordance table compares each method's
estimates with the signal in the bulk RNA. The marker step needs the
stage 5 specificity (runbook 5.1); `run_deconv.sh` now passes those
files when they exist, so a replay of 4.2 after 5.1 gives the same
tables. Inputs brought from the server first:

```bash
rsync -a --include "GTEx_Analysis_v10_Annotations_S*DS.txt" --exclude "*" \
  "$REMOTE:~/Projects/pnatlas-dev/data/raw/gtex_v10_tibial/" "$DATA_ROOT/raw/gtex_v10_tibial/"
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/interim/deconv/gtex_counts.tsv.gz" "$DATA_ROOT/interim/deconv/"
conda activate pnatlas
mkdir -p logs
python workflow/scripts/deconv_evaluate.py --prep-dir "$DATA_ROOT/interim/deconv" --results-dir results/deconvolution \
  --sample-attributes "$DATA_ROOT/raw/gtex_v10_tibial/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt" \
  --subject-phenotypes "$DATA_ROOT/raw/gtex_v10_tibial/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt" \
  --gtex-counts "$DATA_ROOT/interim/deconv/gtex_counts.tsv.gz" --specificity "$DATA_ROOT/interim/genetics/specificity.tsv.gz" \
  --expression "$DATA_ROOT/interim/genetics/expression_cpm.tsv.gz" --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz \
  > logs/deconv_evaluate_local.log 2>&1
cp "$DATA_ROOT/interim/deconv/prepare_summary.tsv" "$DATA_ROOT/interim/deconv/donors.tsv" results/deconvolution/
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && rm -rf results/deconvolution_probe'    # probe outputs, superseded by the full run
```

Check: `results/deconvolution/marker_concordance.tsv` has one row per
method and state with the ten marker genes and the Spearman
correlation; `validation_histology.tsv` has two rows per method and the
"markers" ceiling row; `gtex_covariate_models.tsv` has a method column
and no missing q values.
