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
