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
