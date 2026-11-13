# Stage 6: PNPrank v2 temporal back-test

Plan section 7. Two feature sets through PNPrank's own back-test code at
its submission tag: (A) the published features, (B) the published
features plus the pnatlas gene features. Local variables as in the other
runbooks (`$DATA_ROOT`, `$REMOTE`); the PNPrank checkout is
`$PNPRANK` (locally `../pnprank-dev`, on the server
`~/Projects/pnprank-dev`).

## 6.1 Inputs

PNPrank's static feature matrix is not committed in its repository; it
is rebuilt there by `snakemake stage4_model` from the atlas downloads of
its `config/atlases.yaml`. The copy used here is recorded in the
manifest (`data/README.md`, "Stage 6"); the checkout is detached at the
tag so that the code, configs and truth set are the submitted ones
(`pnprank_backtest.py` checks the tag's commit and that nothing under
`src`, `workflow`, `config` and `data/truth` differs from it):

```bash
PNPRANK=../pnprank-dev
git -C "$PNPRANK" rev-parse submitted-human-genetics^{commit}      # 8cba969...
mkdir -p "$DATA_ROOT/external/pnprank" "$DATA_ROOT/processed"
cp "$PNPRANK/data/processed/features_static.tsv.gz" "$DATA_ROOT/external/pnprank/"
md5 -q "$DATA_ROOT/external/pnprank/features_static.tsv.gz"          # 6fa97af70aabe56790fa0dea20fccc63
diff <(cut -f1-4 "$PNPRANK/data/truth/neuropathy_genes.tsv") <(cut -f1-4 data/metadata/pnprank_truth_set.tsv)   # empty
```

The back-test environment pins the library versions of PNPrank's lock
file (`conda-lock.osx-arm64.txt` at the tag: numpy 2.4.6, pandas 2.3.3,
scikit-learn 1.9.1, scipy 1.17.1, joblib 1.6.0, python 3.11), which the
local `pnatlas` environment happens to match; on the server a separate
environment holds them:

```bash
rsync -a --exclude manuscript --exclude presentations --exclude results --exclude "data/raw" --exclude "data/interim" \
  --exclude "data/processed" --exclude "data/external" --exclude notebooks --exclude __pycache__ --exclude .snakemake \
  "$PNPRANK/" "$REMOTE:~/Projects/pnprank-dev/"
ssh "$REMOTE" 'cd ~/Projects/pnprank-dev && git checkout -q --detach submitted-human-genetics && git rev-parse HEAD'
rsync -a "$DATA_ROOT/external/pnprank/features_static.tsv.gz" "$REMOTE:~/Projects/pnatlas-dev/data/external/pnprank/"
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda create -y -n pnprank-bt -c conda-forge \
  python=3.11 numpy=2.4.6 pandas=2.3.3 scikit-learn=1.9.1 scipy=1.17.1 joblib=1.6.0 pyyaml'
```

## 6.2 Features and back-test

`pnatlas_features.py` builds the pnatlas gene features (docstring: the
state specificity and log CPM, the compartment sums over the mapped
states, the gene by state correlations in the GTEx samples with the
BayesPrism fractions, the colocalisation flags), keyed by the Ensembl
ids of PNPrank's universe. `pnprank_backtest.py` runs
`oof_scores_at_cutoff` of PNPrank for both feature sets at every
positive set, cutoff and seed (PNPrank's config: cutoffs 2005 to 2020
every three years, five scoring folds, three inner folds, seeds
20261006 to 20261015), and writes the per-seed table, the means per
cutoff and the summary with the plan's rule. A probe (one seed, the
first cutoff, tier 1) first, locally, then the full run on the server:

```bash
conda activate pnatlas
python workflow/scripts/pnatlas_features.py --specificity "$DATA_ROOT/interim/genetics/specificity.tsv.gz" \
  --expression "$DATA_ROOT/interim/genetics/expression_cpm.tsv.gz" --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz \
  --status results/spatial/state_mapping_status.tsv --gtex-counts "$DATA_ROOT/interim/deconv/gtex_counts.tsv.gz" \
  --fractions results/deconvolution/gtex_bayesprism.tsv --coloc results/genetics/coloc_all_pairs.tsv \
  --universe "$DATA_ROOT/external/pnprank/features_static.tsv.gz" --out "$DATA_ROOT/processed/pnatlas_gene_features.tsv.gz"
python workflow/scripts/pnprank_backtest.py --pnprank-repo "$PNPRANK" --features-static "$DATA_ROOT/external/pnprank/features_static.tsv.gz" \
  --pnatlas-features "$DATA_ROOT/processed/pnatlas_gene_features.tsv.gz" --out /tmp/backtest_probe --n-jobs 2 --probe
# A at 2005, tier 1: AUROC 0.764, the published value
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && free -g | sed -n 2p && who | wc -l && \
  ~/bin/memguard.sh start logs/memguard.log "pnprank_backtest|pnatlas_features" && \
  setsid -f nohup workflow/scripts/run_backtest.sh ~/Projects/pnprank-dev data results/backtest 24 > logs/backtest.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/backtest.log && ~/bin/memguard.sh stop logs/memguard.log'   # FEATURES 0, BACKTEST 0, BACKTESTDONE
rsync -a --exclude scores "$REMOTE:~/Projects/pnatlas-dev/results/backtest/" results/backtest/
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/processed/pnatlas_gene_features.tsv.gz" results/backtest/
```

Check: `results/backtest/backtest_summary.tsv` has the row with
`primary` true (tier 1, AUROC) and its `improvement` verdict;
`backtest_by_cutoff.tsv` has twelve rows per positive set (six cutoffs
by two feature sets) with A's means equal to PNPrank's published
`backtest_performance.tsv` (fold_rank, all_after, full) at the seed
20261006 row of `backtest_seeds.tsv`.

## 6.3 Checks against PNPrank and clean-up

The published back-test table is read from the PNPrank checkout at the
tag; the check script compares it with feature set A at the global
seed. The per-gene score cache is kept on the laptop data root and
removed from the server once the checksums agree.

```bash
git -C "$PNPRANK" show 'submitted-human-genetics^{commit}:results/backtest/backtest_performance.tsv' > "$DATA_ROOT/external/pnprank/backtest_performance.tsv"
python workflow/scripts/backtest_check_published.py --published "$DATA_ROOT/external/pnprank/backtest_performance.tsv" \
  --seeds results/backtest/backtest_seeds.tsv          # 12 of 12 rows equal the published values
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/backtest/scores/" "$DATA_ROOT/interim/backtest/scores/"
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev/results/backtest/scores && md5sum *' | sort > /tmp/scores_server.md5
(cd "$DATA_ROOT/interim/backtest/scores" && md5sum *) | sort | diff - /tmp/scores_server.md5 && \
  ssh "$REMOTE" 'rm -rf ~/Projects/pnatlas-dev/results/backtest/scores'
```

Check: the script prints `12 of 12 rows equal the published values`
and exits 0; `diff` prints nothing before the cache is removed.
