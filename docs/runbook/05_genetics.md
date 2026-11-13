# Stage 5: genetics

Plan section 6. Gene classes of inherited neuropathy are scored for
enrichment in the cell states of the reference, and the GWAS signals of
acquired neuropathies are colocalised with tibial nerve eQTL. `$REMOTE`
is `user@host` of the server; downloads run there (the data root) and
the small snapshots are fetched on the laptop into `data/metadata/`.

## 5.1 Inputs

Summary statistics and QTL. The five GWAS Catalog studies with full
summary statistics (harmonised files with their tabix index and
metadata), the eQTL Catalogue GTEx tibial nerve gene expression dataset
QTD000286 (all nominal associations with index, permutation results,
SuSiE credible sets and the per-variant log Bayes factors) and the GTEx
v10 eQTL tar are listed in `data/metadata/stage5_files.tsv` (21 files,
14 GB, sizes from HEAD requests) and fetched on the server with the stage
download script, which records them in the manifest; the Ensembl 105
GTF and cDNA (GENCODE v39, the annotation of the GTEx v10 expression
files) are a second list, `data/metadata/stage5_annotation.tsv`,
fetched after the first has finished because both append to the
manifest:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && setsid -f nohup nice -n 15 python3 -I workflow/scripts/download_files.py \
  data/metadata/stage5_files.tsv --data-root data --manifest data/MANIFEST.tsv \
  --fetch "workflow/scripts/download_files.py data/metadata/stage5_files.tsv" > logs/stage5_download.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && grep -c ^done logs/stage5_download.log'   # 21
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && setsid -f nohup nice -n 15 python3 -I workflow/scripts/download_files.py \
  data/metadata/stage5_annotation.tsv --data-root data --manifest data/MANIFEST.tsv \
  --fetch "workflow/scripts/download_files.py data/metadata/stage5_annotation.tsv" > logs/stage5_annotation_download.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && grep -c ^done logs/stage5_annotation_download.log'   # 2
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/MANIFEST.tsv" data/MANIFEST.tsv
```

A file that fails part way (one GWAS file did, the connection dropped)
is finished by rerunning the same list: `curl -C -` resumes the partial
file and the files already in the manifest are skipped. The rerun and
the annotation list were chained in one detached shell so that only one
run appends to the manifest at a time:

```bash
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && setsid -f nohup bash -c "python3 -I workflow/scripts/download_files.py data/metadata/stage5_files.tsv \
  --data-root data --manifest data/MANIFEST.tsv --fetch \"workflow/scripts/download_files.py data/metadata/stage5_files.tsv\" && \
  python3 -I workflow/scripts/download_files.py data/metadata/stage5_annotation.tsv --data-root data --manifest data/MANIFEST.tsv \
  --fetch \"workflow/scripts/download_files.py data/metadata/stage5_annotation.tsv\"; echo STAGE5DOWNLOADDONE" \
  > logs/stage5_download_resume.log 2>&1 < /dev/null'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/MANIFEST.tsv" data/MANIFEST.tsv
```

The md5 of every GWAS file equals the `data_file_md5sum` of its
metadata yaml. The `data/README.md` section "Stage 5: genetics" lists
all 23 files with their md5.

Gene lists. The four PanelApp panels are fetched at the versions the
plan fixes (Genomics England 846 v8.38 and 85 v1.513, PanelApp
Australia 3070 v3.1 and 3439 v2.0; the API serves a version on request)
into `data/metadata/panelapp_<site>_<id>_v<version>.tsv`, committed as
snapshots. The PNPrank truth set is read from the sibling repository at
the tag recorded in stage 1 (runbook 1.16):

```bash
python3 -I workflow/scripts/panelapp_fetch.py --out-dir data/metadata
git -C ../pnprank-dev show submitted-human-genetics:data/truth/neuropathy_genes.tsv > data/metadata/pnprank_truth_set.tsv
```

Check: the four panel files report 376, 288, 307 and 19 genes (231,
159, 249 and 18 green); the truth set has 242 rows.

Gene length and GC. Per gene, the union length of its exons from the
GTF and the GC fraction of its transcripts from the cDNA file, for the
matched background of the enrichment, written as a snapshot of release
105 (on the server once the two files are there, then brought back):

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  nice -n 15 python workflow/scripts/gene_length_gc.py --gtf data/raw/ensembl105/Homo_sapiens.GRCh38.105.gtf.gz \
    --cdna data/raw/ensembl105/Homo_sapiens.GRCh38.cdna.all.fa.gz --out data/metadata/ensembl105_gene_length_gc.tsv.gz'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/metadata/ensembl105_gene_length_gc.tsv.gz" data/metadata/
```

## 5.2 Specificity and gene class enrichment

Specificity of every gene per state from donor pseudobulk of the whole
reference (`workflow/scripts/reference_specificity.py`: counts summed per
donor and state, pseudobulks under `specificity_min_cells` dropped, CPM,
mean over donors, divided by the sum over states), on the server; the
two tables come back to the laptop's data root:

```bash
ssh "$REMOTE" 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas && cd ~/Projects/pnatlas-dev && \
  setsid -f nohup nice -n 15 python workflow/scripts/reference_specificity.py --counts data/interim/reference_merged.h5ad \
    --annotation data/interim/reference_annotation.tsv.gz --out-dir data/interim/genetics > logs/reference_specificity.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && tail -1 logs/reference_specificity.log'   # pseudobulks kept, 31 states, 19420 genes
rsync -a "$REMOTE:~/Projects/pnatlas-dev/data/interim/genetics/" "$DATA_ROOT/interim/genetics/"
```

Enrichment (`workflow/scripts/gene_class_enrichment.py`): for each class
(the five disease categories of the truth set and the whole set, the
green genes of each PanelApp panel) and each state, the mean
specificity of the class genes against 10,000 draws of background genes
matched on exon length and GC (ten quantile bins of each), p with the
+1 correction, Benjamini-Hochberg across states within a class. Runs on
the laptop in a minute:

```bash
python workflow/scripts/gene_class_enrichment.py --specificity "$DATA_ROOT/interim/genetics/specificity.tsv.gz" \
  --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz --truth data/metadata/pnprank_truth_set.tsv \
  --panels-dir data/metadata --results-dir results/genetics
```

Check: `results/genetics/gene_class_enrichment.tsv` has one row per
class and state (10 classes by 31 states) with the `enriched` column;
`gene_classes.tsv` lists each class gene with its state of highest
specificity.

## 5.3 Colocalisation

For each of the five GWAS (`coloc_prepare.py`): lead variants below
`gwas_p` clumped at `window_bp`, GWAS and eQTL rows of each region read
through tabix, variants matched on position and allele pair, genes with
at least 50 shared variants kept as pairs, genes with a SuSiE credible
set in the region flagged and their log Bayes factors pulled from the
lbf file. Then (`coloc_run.R`) coloc.abf on every pair and coloc.bf_bf
between the GWAS single-signal log Bayes factors and the eQTL SuSiE ones
where a credible set exists, priors from the config; `coloc_collect.py`
gathers the pairs of every GWAS, calls a pair colocalised above
`pp_h4`, and gives each colocalised gene its state of highest
specificity and the compartments where that state is mapped (stage 3
status table). All in one script on the server, enrichment included, a
probe with two regions per GWAS first:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && free -g | sed -n 2p && \
  setsid -f nohup workflow/scripts/run_genetics.sh data results/genetics_probe --max-regions 2 > logs/genetics_probe.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/genetics_probe.log'   # ENRICHMENT 0, PREPARE and COLOC 0 for each GWAS, COLLECT 0, GENETICSDONE
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && rm -rf data/interim/genetics/coloc results/genetics_probe && \
  setsid -f nohup workflow/scripts/run_genetics.sh data results/genetics > logs/genetics.log 2>&1 < /dev/null'
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && cat logs/genetics.log'
rsync -a "$REMOTE:~/Projects/pnatlas-dev/results/genetics/" results/genetics/
rsync -a --include "*/" --include "regions.tsv" --include "pairs.tsv" --include "study.tsv" --include "coloc.tsv" --exclude "*" \
  "$REMOTE:~/Projects/pnatlas-dev/data/interim/genetics/coloc/" "$DATA_ROOT/interim/genetics/coloc/"
```

After the probe, `coloc_run.R` was corrected to read the case-control
flag (pandas writes "True") and `coloc_collect.py` to label genes
absent from the reference; the probe's run step was repeated by hand on
one GWAS directory before the full run:

```bash
REMOTE=user@host workflow/remote/sync.sh
ssh "$REMOTE" 'cd ~/Projects/pnatlas-dev && source ~/miniconda3/etc/profile.d/conda.sh && conda activate pnatlas-r &&   Rscript workflow/scripts/coloc_run.R --dir data/interim/genetics/coloc/gwas_mvp_hereditary_idiopathic_pn --out /tmp/cc_probe.tsv   --p1 1e-4 --p2 1e-4 --p12 5e-6 && cut -f1-5,11,12 /tmp/cc_probe.tsv | head'   # gwas_type cc
```

Check: `results/genetics/coloc_summary.tsv` has one row per GWAS with
the regions, pairs tested and colocalised genes; `coloc_colocalised.tsv`
the colocalised pairs with their state assignment.
