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
