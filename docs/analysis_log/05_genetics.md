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
