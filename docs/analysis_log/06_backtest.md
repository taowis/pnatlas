# Stage 6: PNPrank v2 temporal back-test

## 6.1 Design

The question of the plan (section 7) is whether the pnatlas features
improve PNPrank's temporal back-test, with PNPrank's own code at its
submission tag so that nothing in the evaluation changes. Reading that
code before writing anything showed three points where the plan's
wording and the code differ, and each is resolved inside the plan's rule
rather than by changing it.

PNPrank's back-test scores every cutoff with one seed (its global seed,
20261006); the ten seeds the plan names are the ones PNPrank's config
defines for its repeated cross-validation (the global seed and the next
nine). The back-test here runs the same scoring function at each cutoff
with each of those ten seeds, for both feature sets, so that the noise
baseline of the plan (the range of A's mean AUROC over seeds) exists at
every cutoff and is computed exactly as written. The published values
are the seed 20261006 rows. PNPrank's back-test does not compute the
recall of post-cutoff genes in the top one percent; it is computed here
from the same within-fold percentiles that PNPrank's primary transform
produces (a gene is in the top one percent when its percentile is at
least 0.99), which is the plan's secondary metric. PNPrank caches its
per-cutoff scores and reads them back whatever features produced them,
so the scores here are cached by positive set, cutoff, seed and feature
set under a directory of this project.

The pnatlas features, one row per Ensembl gene of PNPrank's universe
(17,154 of its 19,229 genes have them; the rest are absent from the
reference and are imputed by PNPrank's own median imputer with a missing
indicator, as PNPrank does for genes absent from an atlas): the
specificity and log CPM of the gene in each of the 31 reference states;
the sum of its specificity over the states mapped in the ganglion (20
states) and in the nerve (12) in stage 3; the Spearman correlation,
across the 670 GTEx tibial nerve samples, between its log CPM and the
BayesPrism fraction of each of the 16 nerve states (the BayesPrism
fractions because the feature needs variation across samples for every
state, and MuSiC's exact zeros for the rare states would make the
correlation meaningless; decisions.md); and three colocalisation
columns, whether the gene was in any tested pair, its largest PP.H4 (0
when untested) and whether that passes 0.8. The gene class enrichment
of stage 5 is not a feature: it is computed from the truth set and
would leak the labels. The two reference atlases PNPrank already uses
(the Heming sural nerve and the Jung DRG) are among the eleven in its
own features, at PNPrank's harmonised cell types; the pnatlas columns
add the joint reference's finer states and the spatial, bulk and
genetic layers, which is what the comparison is about.

Feature set A is PNPrank's 239 static columns plus its 11 co-expression
label features (recomputed inside every fold from the pre-cutoff genes,
as in PNPrank); B is A plus the 83 pnatlas columns. Both go through the
same logistic model, the same C grid, the same folds (the folds depend
on the labels and the seed only, which the script asserts), the same
percentile transform and the same evaluation (post-cutoff positives of
the set against the 18,982 non-candidate genes). The primary comparison
is tier-1 AUROC averaged over the six cutoffs: for each seed the mean
over cutoffs is taken for A and for B, the noise baseline is the range
of A's means over the ten seeds, and an improvement is claimed only when
B's mean exceeds A's by more than that range and an exact one-sided
sign-flip test over the six cutoff differences gives p below 0.05 (with
six cutoffs the smallest attainable p is 1/64, so the test needs every
cutoff in B's favour, or five with the sixth nearly tied). The same
statistics are reported for the top one percent recall and for tiers 1
and 2 together, as the plan's secondary metrics. A null result is
reported as such.

The probe (seed 20261006, cutoff 2005, tier 1) reproduced PNPrank's
published AUROC for A exactly (0.764), which is the check that the
evaluation is the submitted one.

## 6.2 Results

The full run on the server scored 240 jobs (two positive sets, six
cutoffs, ten seeds, two feature sets) with the feature matrix built
there (17,154 genes by 83 columns, identical in content to the one
from the local probe; 5,026 of the universe genes were in a tested
colocalisation pair and 16 of them pass PP.H4 0.8). The memory guard
logged no warning. Feature set A at seed 20261006 equals PNPrank's
published back-test at every one of the twelve positive set by cutoff
rows: AUROC to three decimals, average precision to four and the
positive counts (`workflow/scripts/backtest_check_published.py`), so
the comparison below is against the submitted evaluation, not a
re-implementation of it.

Primary comparison, tier 1 AUROC averaged over the six cutoffs. A's
mean is 0.7714, and over the ten seeds it ranges from 0.7708 to 0.7721
(a range of 0.0013); B's mean is 0.7784, a difference of 0.0070, more
than five times that range, and B's per-seed means (0.7776 to 0.7790)
do not overlap A's. B is ahead at five of the six cutoffs (2005 by
0.007, 2008 by 0.003, 2014 by 0.004, 2017 by 0.018, 2020 by 0.013) and
behind at 2011 (by 0.003); the exact sign-flip test over the six
differences gives p 0.031. Both conditions of the plan's rule hold, so
the pnatlas features improve the temporal back-test on its primary
metric. The gains are concentrated at the late cutoffs, where the
models train on the most pre-cutoff positives and are evaluated on the
fewest post-cutoff ones (19 at 2017 and 13 at 2020), so the late
differences rest on small evaluation sets even though the ten seeds
agree on their direction.

Secondary metrics. The recall of post-cutoff tier 1 genes in the top
one percent is lower with the pnatlas features: 0.127 against 0.146
averaged over cutoffs, behind at five cutoffs and tied at 2020, and
the loss (0.019) is larger than A's across-seed range (0.008). Average
precision is unchanged (0.0137 against 0.0142, within the range). For
tiers 1 and 2 together the AUROC difference is 0.0022 with B ahead at
three cutoffs only (p 0.20), which the rule reads as null, and the top
one percent recall and average precision are null as well. The
improvement is therefore in the ordering of the whole ranking for the
strict positives, not in how many of the later discoveries reach the
top of the list; the top of the list is where PNPrank's recall is
measured for prioritisation, and the new features do not help there.

Two details worth recording. The inner cross-validation chose C 0.001,
the smallest value of PNPrank's grid, in 599 of 600 fits for A and in
all 600 for B, so both models are heavily regularised and the 83 added
columns enter at small weights; the grid is PNPrank's and was not
extended, because changing it would change the evaluation. And the
feature matrix the back-test used is committed as
`results/backtest/pnatlas_gene_features.tsv.gz`, with the per-seed,
per-cutoff and summary tables beside it.

Clean-up on the server. The per-gene score cache of the run (240
files, 58 MB, one per positive set, cutoff, seed and feature set) was
copied to the laptop data root (`interim/backtest/scores/`, checksums
identical) and deleted on the server; a README in the server's
`results/backtest/` says so. The PNPrank checkout at the tag (132 MB)
and its static feature matrix stay on the server, since any rerun of
this stage needs both; the `pnprank-bt` environment stays with them.
