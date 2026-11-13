# Stage 5: genetics

Plan section 6; commands in `docs/runbook/05_genetics.md`.

## 5.1 Inputs and design

Gene lists. The PNPrank truth set is read at the tag of the submitted
version (242 genes; CMT 59, dHMN 19, HSAN 16, SYND 124, HSP 24). The
four PanelApp panels were fetched at the versions the plan fixes; the
Genomics England panels and the autonomic panel are still at those
versions, while PanelApp Australia's hereditary neuropathy panel has
moved on to version 3.3, so version 3.1 was requested from the API
explicitly, and the four files are committed as snapshots (376, 288,
307 and 19 genes; 231, 159, 249 and 18 green).

Background for the enrichment. The plan asks for length and GC matched
background genes. The Ensembl BioMart archive for release 105 (GENCODE
v39, the annotation of the GTEx v10 expression files) no longer answers
BioMart queries, so length and GC are computed from the release 105 GTF
and cDNA files: the union length of a gene's exons, and the GC fraction
over its transcripts. The background is the protein-coding genes with a
specificity, a length and a GC; it is cut into ten quantile bins of
each, and every bootstrap draw replaces each class gene by a background
gene from the same length by GC bin.

Specificity. Donor by state pseudobulks of the whole reference (689
with at least 20 cells, over 31 states), CPM, mean over donors, divided
by the sum over states, so a gene's specificities sum to one. Cells of
unknown donor form one pseudobulk per dataset and state.

Colocalisation. The GWAS Catalog harmonised files (GRCh38, tabix
indexed) of the five studies and the eQTL Catalogue tibial nerve gene
expression dataset (QTD000286: all nominal associations with index, the
SuSiE credible sets and the per-variant log Bayes factors) are the
inputs, 14 GB in all. Lead variants are clumped by distance (500 kb);
each region is tested against every gene with at least 50 shared
variants. coloc.abf runs on every pair; where the gene has a SuSiE
credible set in the region, coloc.bf_bf combines the GWAS single-signal
log Bayes factors (from coloc.abf's Wakefield approximation) with the
eQTL SuSiE log Bayes factors, and the largest PP.H4 over signal pairs is
the reported value; otherwise the abf value is. The GWAS is handled as
case-control with the study's case fraction when its metadata gives
one. The colocalised eGenes are assigned to the state of their highest
specificity, with the compartments where that state is mapped.

The harmonised files turned out to be in the GWAS-SSF layout
(chromosome, base_pair_location, effect_allele, other_allele, beta or
odds_ratio, standard_error, effect_allele_frequency, p_value, rsid),
not the older hm_ columns the preparation script was first written
for, and the script was rewritten before any run. Two of the five
studies (the MVP PheCode scans) give odds ratios with no standard
error, only the 95 percent interval, so the effect is the log odds
ratio and the standard error comes from the interval; they also give
the case and control counts per variant, and the case count for the
case-control model is the median over variants. The carpal tunnel
study's metadata lists each sample set twice, and identical entries
are counted once (480,201 European plus 178,726 East Asian). The
diabetic neuropathy study (a UK Biobank chronic pain scan run in
Regenie) and the carpal tunnel meta-analysis give no case count, so
they are treated as quantitative in coloc, which is the usual fallback
when the case fraction is unknown. The harmonisation code is not used
as a filter: alleles are matched as an unordered pair and coloc works
on squared effects, so strand and orientation do not matter for the
posterior.

## 5.2 Gene class enrichment

Background: 17,227 protein-coding genes with a specificity, a length and
a GC (of 19,951 protein-coding genes in release 105), in 100 length by
GC bins. Of the truth set 221 genes are in the background (CMT 53, dHMN
18, HSAN 14, SYND 113, HSP 23); of the PanelApp green genes 213, 146,
226 and 16.

The classes land where their biology says. CMT is enriched in
myelinating Schwann cells (3.2 fold, q 0.003) and in Ab-LTMR neurons
(1.6 fold). HSAN is enriched in the peptidergic neuron states (C-PEP
2.5 fold, PEP 2.7, C-NP 2.5, A-PEP 2.1), the nociceptors whose loss
defines the disease. dHMN is enriched in Ab-LTMR (2.0 fold) and C-LTMR
(2.4) neurons, the only neuron states of the reference, which stand in
for the motor neurons the reference does not hold. HSP is enriched in
Ab-LTMR (2.0 fold). SYND, the syndromic class, is enriched nowhere,
as expected for genes whose expression is broad. The whole truth set
and the three hereditary neuropathy panels agree with each other:
enriched in myelinating Schwann cells, the injured neuron state and the
Ab-LTMR, C-LTMR, PEP, C-PEP, A-PEP and C-NP families, at folds of 1.2
to 1.7, and in nothing else. The autonomic neuropathy panel (16 genes)
is enriched in the sympathetic neuron state alone, 3.8 fold. No class
is enriched in any nerve or ganglion non-neuronal state other than
myelinating Schwann cells, and none in the non-myelinating, repair or
damage Schwann states, fibroblasts, perineurial cells, immune or
vascular states. Full table: `results/genetics/gene_class_enrichment.tsv`;
each class gene with its state of highest specificity:
`results/genetics/gene_classes.tsv`.

## 5.3 Colocalisation probe

The probe (two lead regions per GWAS, the whole chain from enrichment
to the collector) ran through. The lead regions of the full scans are
few for the smaller studies: the diabetic neuropathy scan has one
region (the HLA, chromosome 6 at 32.6 Mb) and the neuropathic pain scan
one (chromosome 11 at 113.0 Mb, near DRD2), the hereditary and
idiopathic neuropathy scan 22, the type 2 diabetes with neurological
manifestations scan 241 (the diabetes loci, TCF7L2 first). In the probe
regions 38 to 140 genes shared at least 50 variants with the GWAS, and
3 to 24 of them had a SuSiE credible set, for which the log Bayes
factors were written and coloc.bf_bf ran.

Two faults were found and fixed before the full run. The case-control
flag written by pandas is the text "True", which R does not read as a
logical, so every study had been run as quantitative; with the flag
read, the two MVP scans run as case-control with their case fraction
(the probe's one colocalised pair, ILRUN at the chromosome 6 region of
the hereditary and idiopathic neuropathy scan, kept PP.H4 0.81 under
both). And a colocalised gene whose Ensembl 105 symbol differs from the
reference's (ILRUN was C6orf106 in the reference's annotation) had an
empty state; the collector now says "not in reference" so that the
case is visible, and such genes are looked up by hand in the reading of
the results.

## 5.4 Colocalisation results

Lead regions after clumping at 500 kb: 1 for the diabetic neuropathy
scan, 1 for neuropathic pain, 17 for carpal tunnel syndrome, 22 for the
hereditary and idiopathic neuropathy PheCode scan and 241 for the type
2 diabetes with neurological manifestations scan. They gave 9,858
region by gene pairs with at least 50 shared variants, 1,112 of them
with a SuSiE credible set for the eGene in the region (tested with
coloc.bf_bf), the rest with coloc.abf. 81 pairs reach PP.H4 above 0.5
and 27 pass the 0.8 threshold, in 21 regions, all from the two MVP
scans; the three smaller scans (one or 17 regions) give none. The full
table is `results/genetics/coloc_all_pairs.tsv`, the colocalised pairs
with their state assignment `coloc_colocalised.tsv`, and the
specificity of those genes across all states `coloc_gene_specificity.tsv`.

Hereditary and idiopathic neuropathy: one eGene, ILRUN, colocalises in
two neighbouring regions of chromosome 6p21 (PP.H4 0.81 with the SuSiE
signal at the 34.76 Mb lead, 0.98 with coloc.abf at the 35.42 Mb lead).
ILRUN is not in the reference specificity (it was C6orf106 in older
annotations and is below the expression filter of the pseudobulk
table under either name), so it has no state in this reference.

Type 2 diabetes with neurological manifestations: 25 pairs in 19
regions, and they are the diabetes loci rather than nerve loci: WFS1
(Wolfram syndrome gene, PP.H4 0.83, highest specificity in Ab-LTMR
neurons at 0.17), JAZF1 (0.96), TFAP2B (0.81, the one clearly
state-specific gene, sympathetic neuron at 0.78), SCN2A (0.98, highest
in the sympathetic neuron state at 0.22), CCND2-AS1 (1.00), CBFA2T3
(0.95, granulocyte 0.19), CEP68 (0.97), NLRC3 (0.94, T/NK 0.40), PGAP3
and PNMT at the ERBB2 locus, DVL3, CLUAP1, PLEKHH3 and CNTNAP1 at
17q21, RALY-AS1 in two regions, DNAH10OS, H1-7, RMC1, NIPSNAP2,
RCCD1-AS1 and four unnamed genes (ENSG00000286176 next to WFS1,
ENSG00000236283 next to SCN2A, ENSG00000234936, ENSG00000275228). Most
of these have a flat specificity (highest state below 0.2), which is
expected for broadly expressed genes whose eQTL is shared across
tissues; the state column is the one of highest specificity as the plan
says and is not a claim of specificity. CNTNAP1 (a paranodal protein,
the gene of a congenital hypomyelinating neuropathy) at 17q21 is the
one colocalised gene with a direct nerve role; its top state in the
reference is injured neuron at 0.10.

The reading: the nerve eQTL colocalisation of these GWAS is thin, which
follows from the scans themselves (the two non-MVP neuropathy scans
have one lead each), not from the eQTL dataset (1,112 pairs had a
credible set). The gene class enrichment (5.2) is the stronger genetic
link between the reference and neuropathy; the colocalised genes go to
PNPrank v2 as a feature in stage 6 as the plan says, with their PP.H4
and state.

The per-region GWAS and eQTL extracts on the server
(`data/interim/genetics/coloc`, 2 GB) were deleted after the small
tables of each GWAS (regions, pairs, study, coloc) had been copied to
the local data directory; `coloc_prepare.py` rebuilds them from the raw
files, which stay on the server and are in the manifest. The memory
guard recorded no warning during the run.
