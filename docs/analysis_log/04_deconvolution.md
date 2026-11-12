# Stage 4: GTEx tibial nerve deconvolution

Plan section 5; commands in `docs/runbook/04_deconvolution.md`.

## 4.1 Design

The nerve part of the reference is one dataset, the Heming 2025 sural
nerve nuclei (365,708 cells with a named state, 37 donors, 16 states
from B cells and granulocytes with a few hundred cells to perineurial
cells with 84,000). The plan asks for donor pseudobulk signatures as
the reference and for BayesPrism and MuSiC as the methods. Both methods
take single cells and build their own reference from them, BayesPrism a
state profile with a per-sample correction, MuSiC a profile with
cross-donor variance weights that need the donor of every cell, so the
cells are handed over as cells with their state and donor, and the
"donor pseudobulk" of the plan is what MuSiC computes internally. The
reference given to the methods is a random subsample of up to 1,000
cells per state from the reference donors (16,000 cells at most), which
keeps BayesPrism's dense matrices to a few gigabytes; the states with
fewer cells contribute all of them.

Validation (a) holds out 7 of the 37 donors, chosen by the project seed
(S02, S08, S15, S22, S23, S26, S30), before the reference is built. Each
of the 200 pseudo-bulk mixtures takes one held-out donor in turn, draws
its state proportions from a flat Dirichlet, draws 500 cells of that
donor with those proportions (with replacement; a state the donor lacks
is drawn from the other held-out donors) and sums their counts, so the
mixtures carry donor effects the reference has not seen, as GTEx samples
will. The truth is the realised cell fractions. The plan's statistic is
the mean over mixtures of the Pearson correlation across states, with the
RMSE over all entries; the per-state correlation across mixtures is
recorded as well because a flat Dirichlet gives every state the same
range and the mean across states hides which ones are estimated.

The GTEx side uses the v10 per-tissue read counts of "Nerve - Tibial"
(GENCODE v39 gene ids mapped to symbols through the file's Description
column, duplicate symbols summed), the open sample attributes (ischaemic
time, pathology notes) and the open subject phenotypes (sex, age
bracket, Hardy scale). Validation (b) reads the pathology notes for
mentions of adipose tissue and of fibrosis and compares the adipocyte
and the fibroblast (endoneurial plus epineurial) fractions between the
samples with and without the mention, one-sided in the expected
direction. The covariate model is ordinary least squares of each state's
proportion on age (bracket midpoint), sex, ischaemic time and the Hardy
scale as a factor, all together, with Benjamini-Hochberg across states
within each term.

Both methods come from GitHub at fixed commits because neither
repository has release tags (BayesPrism 2.2.3, MuSiC 1.0.0); their
dependencies are in `environment-r.yml`. A probe with ten mixtures, 100
cells per state and 20 GTEx samples runs the chain before the full run.
The prepare step found 19,385 genes shared between the reference and
the GTEx symbols.

## 4.2 Probe

The probe (ten mixtures of the held-out donors, 100 cells per state,
20 GTEx samples) ran the whole chain. Two faults were fixed on the way:
the BayesPrism version installed (2.2.3) names the outlier gene
fraction argument differently from its documentation, and the covariate
model has to drop a Hardy level absent from a small sample set. On the
probe MuSiC reached a mean Pearson of 0.83 (RMSE 0.03) and BayesPrism
0.73 (RMSE 0.05), with 100 cells per state only; the full run decides.
BayesPrism took about 13 seconds per bulk sample on 12 cores at that
reference size.

## 4.3 Results

The nerve reference for the deconvolution has 16 states (every nerve
state with at least 200 cells; the DRG and ganglion states are not in
the nerve compartment) and 14,549 cells from 30 donors; the 7 held-out
donors give 200 mixtures of 500 cells. The GTEx table has 670 tibial
nerve samples, 660 with every covariate (ten have no Hardy scale).

Validation (a), on the mixtures: MuSiC passes (mean Pearson across
states within a mixture 0.84, median 0.89, RMSE 0.034) and BayesPrism
does not (0.73, median 0.77, RMSE 0.050, against the thresholds of 0.8
and 0.1). Per state, BayesPrism's failure comes from adipocyte
(estimated mean 0.136 against a true 0.061, RMSE 0.107) and from the
granulocyte and non-myelinating Schwann fractions, which it halves.
Neither method recovers the damage Schwann state (Pearson across
mixtures about -0.08 for both; MuSiC estimates it at 0.009 against a
true 0.056): in a pseudobulk its profile is not separable from the
other Schwann states, and it is left out of the interpretation. The
other 15 states are recovered by MuSiC at Pearson 0.84 to 0.95 across
mixtures.

By the rule written before the run (BayesPrism unless it fails and
MuSiC passes) MuSiC is the primary method, and `gtex_proportions.tsv`
carries its estimates. Both methods' GTEx estimates are kept, because
they disagree in a way the mixtures do not show. MuSiC is a weighted
non-negative least squares and returns exact zeros: in the GTEx samples
it gives zero adipocyte in 77 percent, zero granulocyte in 97, zero
mast in 86, zero lymphatic endothelial in 84 and zero repair Schwann in
90 percent. The marker concordance (each method's fraction against the
mean log CPM of the state's ten most specific expressed protein coding
genes, across the 670 samples) shows BayesPrism tracking the bulk
signal better for 13 of the 16 states: adipocyte 0.72 against 0.49,
mast 0.87 against 0.19, granulocyte 0.70 against 0.20, T/NK 0.61
against -0.08, B 0.58 against 0.51; the two methods agree on the large
states (myelinating Schwann 0.96 and 0.93, epineurial fibroblast 0.86
and 0.85, mural 0.93 and 0.83, endothelial 0.70 and 0.69), and MuSiC is
better for non-myelinating Schwann (0.87 against 0.66). The mixtures
are built from nuclei of the same technology as the reference, which is
the setting MuSiC fits well; the bulk is a different platform, which is
what BayesPrism's reference update is built for. The rule is kept
because it was fixed in advance, but the rare states should be read
from BayesPrism, and the MuSiC zeros for them are not estimates.

Validation (b), the histology notes, fails for both methods, and the
reason is in the notes. The pathologists wrote a fat percentage in 393
of the 670 notes (median 10 percent, range 0 to 90; "no fat" and
"trace fat" read as 0), and the adipocyte fraction does not follow it
(Spearman 0.05 for BayesPrism, 0.04 for MuSiC). But the adipocyte
marker genes themselves follow it only weakly (the ten-gene score at
0.15, ADIPOQ, PLIN1, LEP and GPAM singly at 0.12 to 0.16): the notes
describe the piece as received, and many say the fat was "delineated",
marked for removal before extraction, so the percentage is not the fat
in the RNA. The ceiling row is kept in the table. Fibrosis is mentioned
in 27 notes only and the fibroblast fraction does not differ. The
histology check is therefore uninformative rather than failed, and the
marker concordance is the evidence that the estimates carry signal.

The covariate models (660 samples, age, sex, ischaemic time and the
Hardy scale together) give the same leading findings with both methods:
the non-myelinating Schwann fraction rises with age (MuSiC 0.0017 per
year, BayesPrism 0.0021, q below 1e-17 for both) and the myelinating
Schwann fraction falls (-0.0011 and -0.0015 per year, q below 1e-8),
which is the known loss of myelinated fibres with age in sural nerve
morphometry; the myelinating Schwann fraction is higher in women (0.016
and 0.023, q below 0.003); the endothelial fraction is lower in Hardy 2
deaths (fast natural deaths, -0.026 and -0.031, q below 1e-5) and falls
with age in MuSiC; the endoneurial fibroblast and the macrophage
fractions rise with age in both. MuSiC gives 23 and BayesPrism 26 state
by term pairs at q below 0.05; the full tables are
`results/deconvolution/gtex_covariate_models.tsv`.

The probe outputs on the server (`results/deconvolution_probe`) were
deleted after the full run; the full results and the small interim
tables are in the repository and the local data directory. The memory
guard recorded no warning during the full run.
