# Analysis plan

Frozen at the end of stage 1. The text below is not edited after freezing;
any later change is added under "Amendments" with its number and reason, and
the stage reports say which amendment they follow.

## 1. Datasets per layer

The dataset list is `config/datasets.yaml`; verification records are in
`data/metadata/` and `docs/datasets.md`.

Inclusion: human tissue; ganglia (DRG, trigeminal, sympathetic) or peripheral
nerve; processed data (counts or a processed object) obtainable without an
access application, or through a public portal that states what it provides.
Non-human data (macaque GSE165569, mouse series) are used only as anchors for
cell-type labels, never as donors in the reference.

Exclusion: datasets under embargo at the time of the stage that would use
them (SPARC 464 until its release); raw data under controlled access unless
and until access is granted, in which case the addition is an amendment;
datasets with fewer than two donors in a layer where donors are the unit of
inference (they can still be used for marker checks); datasets with no stated
licence and no author permission are used for analysis with citation but
nothing derived from them at cell level is redistributed.

Layer 1 reference: the Bhuiyan 2025 preprint object is the DRG backbone;
Nguyen 2021, Jung 2023 and Yu 2024 add DRG donors where they are not already in
that object (checked by donor ID and study); Yang 2022 for trigeminal; SPARC
465 for sympathetic ganglia; Heming 2025 (GSE285983) for nerve. Layer 2
spatial: Heming 2025 Xenium (GSE285984), Zhang 2024 and Vroman 2023 Visium,
Yu 2024 Xenium, Shiers 2025 DPN DRG, and the Tavares-Ferreira 2022 DRG Visium
and DPN sural nerve Visium processed objects from the UT Dallas Box folders
(no licence stated: analysis with citation, no redistribution). Layer 3 bulk: GTEx v10 tibial nerve
(670 samples) as the primary cohort; Guo 2020 and Hur 2011 sural nerve DPN as
secondary cohorts. Layer 4 genetics: PanelApp panels and the PNPrank truth set
for Mendelian genes; GWAS Catalog studies with full summary statistics for
acquired neuropathy and nerve traits; GTEx v10 tibial nerve cis-eQTL and the
eQTL Catalogue GTEx tibial nerve datasets for colocalisation.

## 2. Unit of inference

The donor. Cells, nuclei, spots and sections are not replicates. Differential
expression and abundance comparisons use donor-level pseudobulk (sum of
counts per donor per cell state) with DESeq2 or edgeR, or a mixed model with
donor as a random effect when a donor contributes several sections or regions.
A cell state is reported only if it is found in at least three donors from at
least two studies, unless the tissue has fewer donors available in total, in
which case the limitation is stated.

## 3. Reference construction

Integration: scVI with donor as the batch key and study as a covariate,
followed by scANVI seeded with the author labels mapped through
`config/celltype_map.tsv`. Harmony on the same PCA space is run as the
comparison method.

Choice between methods is made once, by a fixed rule, before any downstream
analysis: the scIB metrics (batch ASW, iLISI, graph connectivity for batch
removal; cell-type ASW, NMI and ARI against the mapped author labels for
biological conservation) are computed for both embeddings on the ganglia and
the nerve subsets separately, and the method with the higher mean biological
conservation score is kept provided its batch removal score is within 0.05 of
the other method's; otherwise the method with the higher overall scIB score is
kept.

Annotation: DRG neuron and non-neuron names follow the Bhuiyan 2025
nomenclature; trigeminal and sympathetic labels are mapped to it where the
type exists and kept as their own names otherwise; nerve labels follow
Heming 2025. `config/celltype_map.tsv` records every author label and its
unified label and is the only place where that mapping is changed. Labels at
the level used for genetics and deconvolution (the "state" level) must each
contain at least 50 cells from at least three donors.

## 4. Spatial mapping

Visium sections: cell2location with reference signatures estimated from the
joint reference (negative binomial regression per cell state, donor and study
as batch), default priors, 30 cells per spot as the expected cell count
unless the section's nuclear count says otherwise.

Xenium and other targeted panels: scANVI label transfer on the genes shared
between the panel and the reference, after confirming that at least 60 shared
genes remain; cells with a label probability below 0.5 are "unassigned" and
reported as such.

Checks that a mapped state is real, applied per section and reported for
every section used: (a) marker agreement, the AUROC of the state's top 20
reference markers for spots or cells assigned to that state versus the rest,
required above 0.7; (b) held-out genes, 10 percent of the shared genes are
left out of the mapping and their observed expression per state is correlated
with the reference prediction, required Spearman above 0.5; (c) anatomy, the
perineurial, epineurial and endoneurial states must lie where the H&E or DAPI
image shows those compartments, judged by one reader and recorded as a yes or
no per section. A state that fails (a) or (b) in more than half of the
sections is reported as not mapped.

## 5. Deconvolution

Reference: the nerve part of the joint reference at the state level, donor
pseudobulk signatures. Method: BayesPrism as the primary method, MuSiC as the
secondary; the two are compared on the validation sets and the primary is
kept unless it fails validation.

Validation before any GTEx result is interpreted: (a) pseudo-bulk mixtures
built from held-out donors of GSE285983 with known proportions, 200 mixtures,
required mean Pearson correlation above 0.8 between estimated and true
proportions across states and RMSE below 0.1; (b) agreement with histology,
GTEx tibial nerve pathology notes (open, on the portal) that mention adipose
or fibrotic tissue are compared with the estimated adipocyte and fibroblast
fractions, required direction as expected with a Wilcoxon test. The GTEx
analysis then tests cell-state proportions against age, sex, ischaemic time
and death classification (Hardy scale) in a linear model with all four
together. Diabetes and BMI are tested only if dbGaP phs000424 access is
granted, as an amendment.

## 6. Genetics

Gene classes: CMT, dHMN, HSAN and syndromic neuropathy from the PNPrank truth
set (242 curated genes, tiers 1 and 2, categories CMT, dHMN, HSAN, SYND, HSP);
PanelApp Genomics England panel 846 v8.38 green genes, panel 85 v1.513, and
PanelApp Australia panel 3070 v3.1 and 3439 v2.0 (autonomic) as the broader
lists; version numbers are fixed at these values. OMIM is not used because it
needs a licence key.

Enrichment: for each gene class and each cell state, EWCE-style bootstrap
enrichment of specificity scores computed from donor pseudobulk, 10,000
bootstraps with length and GC matched background genes, Benjamini-Hochberg
across states within a class; a class is "enriched" in a state at adjusted
p below 0.05.

Colocalisation: GWAS Catalog studies with full summary statistics
(GCST90475844 hereditary and idiopathic peripheral neuropathy, GCST90475676
type 2 diabetes with neurological manifestations, GCST90129436 diabetic
neuropathy, GCST90859373 neuropathic pain, GCST90018813 carpal tunnel) against
GTEx v10 tibial nerve cis-eQTL and the eQTL Catalogue GTEx tibial nerve
dataset QTD000286 (SuSiE credible sets available). Method: coloc with SuSiE
where credible sets exist and coloc.abf otherwise, priors p1 1e-4, p2 1e-4,
p12 5e-6, a colocalised pair at PP.H4 above 0.8. Cell-state assignment of a
colocalised eGene uses its specificity in the reference and its mapped
spatial compartment. Cell-state interaction eQTL in GTEx needs the donor
genotypes (dbGaP phs000424) and is added only as an amendment if access is
granted.

## 7. Stage 6 primary comparison

Fixed before it is run. Two feature sets are compared with identical
evaluation code from the PNPrank submission tag (`submitted-human-genetics`,
commit 8cba969 of pnprank-dev): (A) the published PNPrank
features; (B) the published features plus the pnatlas features (nerve
cell-state specificity, spatial compartment scores, GTEx deconvolution
derived nerve expression, colocalisation flags). Same truth set (PNPrank
`data/truth/neuropathy_genes.tsv` at the submission tag), same cutoff years
(2005, 2008, 2011, 2014, 2017, 2020), same temporal splits and the same 10
random fold seeds.

Primary metric: tier-1 AUROC in the temporal back-test, averaged over cutoff
years. Noise baseline: the spread of (A) over the 10 fold seeds. An
improvement is claimed only if the mean AUROC of (B) exceeds that of (A) by
more than the full range of (A) across seeds, and a paired permutation test
over the six cutoff years gives p below 0.05. Secondary metrics: recall of
post-cutoff tier-1 genes in the top 1 percent, and the same two statistics on
tier 2. A null result (no improvement, or an improvement within the noise
baseline) is reported as such and is a finding of the paper.

## Amendments

1. (Stage 2, step 2.3) The painseq humanDRG object is the harmonized DRG
   neuron atlas (53,596 neurons, nine source datasets, 22 subtype labels,
   no recoverable donors), not a 126-donor whole-tissue dataset. It is used
   as the neuron label source and label-transfer reference; its cells that
   come from Jung 2023, Nguyen 2021 and Tavares-Ferreira 2022 are excluded
   before integration, and its remaining neurons enter with library as the
   unit and donor marked unknown. The DRG non-neuronal backbone is Jung 2023
   (7 donors) with Nguyen 2021 and Zhang. Section 1 is read with this
   change; the donor rule in section 2 is unchanged.
