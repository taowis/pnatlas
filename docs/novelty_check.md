# Human PNS (DRG + nerve) four-layer resource: novelty check and public data inventory

Checked 2026-11-02. Every paper below was retrieved live from PubMed (E-utilities), Europe PMC or the bioRxiv API, and the abstract was read. GEO accessions were checked with NCBI esummary (db=gds), GTEx counts with the GTEx portal API (v2), GWAS IDs with the GWAS Catalog REST API v2, and SPARC/Zenodo DOIs with DataCite. Numbers in square brackets are the record IDs in `novelty_refs.ris`. When I could not check something directly, it is marked UNVERIFIED.

## 1. Short verdict

Most of layer 1 (the harmonised single-cell reference) has already been done, and layers 2 and 4 are partly done for DRG. Nobody has yet put all four layers together for the DRG and the peripheral nerve as one system.

- A harmonised cross-species DRG + trigeminal atlas with a portal already exists: Bhuiyan 2024 Sci Adv [7]. It covers 6 species, 31 datasets, 18 neuronal and 11 non-neuronal types, and has a portal at harmonized.painseq.com, which is live.
- A much larger human DRG reference atlas came out as a preprint in Nov 2025: Bhuiyan 2025 bioRxiv [21]. It has 126 donors, covers cervical, thoracic and lumbar levels, resolves 22 neuronal and 10 non-neuronal types, and adds cross-species integration plus Xenium and MERSCOPE spatial data. It already has a public portal (humanDRG.painseq.com, "NIH PRECISION Human Pain Network DRG Atlas", live, 53,596 neuronal cells/nuclei), and a Seurat object can be downloaded from that page. The preprint says raw and processed data will go to dbGaP/SPARC. Its bioRxiv record does not yet show a journal publication.
- Mapping inherited neuropathy genes onto human ganglion cell types has also started. Yang 2025 bioRxiv [22] (Gereau lab) built a snMultiome atlas of human sympathetic ganglia + DRG and mapped disease genes, with an explicit HSAN section and a Shiny portal. The companion DRG multiome is Zhao/Boyer 2025 [23], on SPARC.
- Pain GWAS has already been mapped onto human DRG cell types. Toikumo 2025 JCI [13] took a chronic pain GWAS of more than 1.2 million people [60] onto hDRG scRNA-seq and found enrichment in hPEP.TRPV1/A1.2. Kupari 2021 [2] did the same with macaque neurons.
- There is a human sural nerve snRNA-seq + Xenium atlas across polyneuropathies: Heming 2025 Nat Commun [29]. It has 33 patients + 4 controls and 365,708 nuclei, open data on GEO, and a portal at pns-atlas.mzhlab.com.
- I found nothing that deconvolves GTEx tibial nerve with a human nerve single-cell reference. Donovan 2020 [47] deconvolved 28 GTEx tissues using mouse signatures, and the word "nerve" does not appear anywhere in its full text. Kim-Hellmuth 2020 [48] estimated 7 broad cell types across GTEx tissues. Ray 2019 [39] analysed sex differences in 248 GTEx tibial nerve samples, but without deconvolution.
- The sural nerve layer starts from an existing local re-analysis of Heming 2025 [29] (GEO GSE285983) with inherited neuropathy gene mapping. The closest published items are Heming 2025 [29] (nerve, acquired polyneuropathies) and Yang 2025 [22] (ganglia, HSAN genes).

The open space is therefore not "a DRG atlas". It is a joint DRG + nerve reference across all inherited neuropathy classes, used to deconvolve GTEx tibial nerve, with eQTL colocalisation. Details are in section 3.

## 2. Existing work, by layer

### 2a. Human sensory ganglia single-cell atlases

| # | First author, year | Title (short) | Journal | PMID | DOI | What it covers (from abstract) |
|---|---|---|---|---|---|---|
| 1 | Nguyen, 2021 | Single-nucleus transcriptomic analysis of human DRG neurons | eLife | 34825887 | 10.7554/eLife.71752 | Human DRG snRNA-seq; neuronal classes vs mouse; several human classes with no clear mouse match. NIH Ryba/Davidson group. GEO GSE168243. |
| 2 | Kupari, 2021 | Single cell transcriptomics of primate sensory neurons... chronic pain | Nat Commun | 33686078 | 10.1038/s41467-021-21725-z | Macaque DRG neurons (Ernfors lab); mapped chronic pain GWAS loci onto neuron types. GEO GSE165569. |
| 3 | Tavares-Ferreira, 2022 | Spatial transcriptomics of DRG identifies molecular signatures of human nociceptors | Sci Transl Med | 35171654 | 10.1126/scitranslmed.abj8186 | Visium-based near-single-neuron profiles from 8 organ donors; 12 neuronal clusters; sex differences. Raw data on dbGaP phs001158; processed data on sensoryomics.com. |
| 5 | Jung, 2023 | Cross-species transcriptomic atlas of DRG | Nat Commun | 36690629 | 10.1038/s41467-023-36014-0 | snRNA-seq of mouse, guinea pig, cynomolgus and human DRG (Genentech); human from 6 + 1 donors. GEO GSE201654 / GSE201586. |
| 6 | Yang, 2022 | Human and mouse trigeminal ganglia cell atlas... migraine | Neuron | 35349784 | 10.1016/j.neuron.2022.03.003 | Human + mouse TG snRNA-seq and snATAC-seq; migraine GWAS mapped onto cell types. GEO GSE197289. |
| 7 | Bhuiyan, 2024 | Harmonized cross-species cell atlases of TG and DRG | Sci Adv | 38905344 | 10.1126/sciadv.adj9173 | 31 datasets, 6 species, 18 neuronal + 11 non-neuronal types; new human (10 donors) and axolotl DRG data. Portal: harmonized.painseq.com. |
| 8 | Yu, 2024 | Deep single-soma RNA-seq... human somatosensation | Nat Neurosci | 39496796 | 10.1038/s41593-024-01794-1 | Single-soma Smart-seq of 1,136 hDRG neurons, 16 types; Xenium validation; microneurography. GEO GSE249746 (soma), GSE273557 (Xenium). |
| 9 | Lu, 2024 | Developing human sensory neurons and organoid modelling | Cell | 39536745 | 10.1016/j.cell.2024.10.023 | Embryonic human DRG single-cell + spatial atlas; organoids. |
| 10 | Arendt-Tranholm, 2026 | Transcriptome of human C2 DRG in C1-2 arthrodesis | Brain | 41032656 | 10.1093/brain/awaf370 | Bulk, snRNA-seq and spatial data from C2 DRG of 22 neck-pain surgical patients, compared with lumbar DRG and TG. |
| 14 | Avraham, 2022 | Satellite glia single-cell, rodent vs human | Pain | 35503034 | 10.1097/j.pain.0000000000002628 | SGC profiles in rodent and human DRG. GEO GSE169301. |
| 21 | Bhuiyan, 2025 (preprint) | A Reference Atlas of the Human DRG | bioRxiv | 41279342 | 10.1101/2025.11.05.686654 | 126 donors, cervical/thoracic/lumbar; 22 neuronal + 10 non-neuronal types; cross-species; Xenium (472 genes) + MERSCOPE (300 genes); microneurography. Portal humanDRG.painseq.com; data "will be" on dbGaP/SPARC. |
| 22 | Yang, 2025 (preprint, v2 Feb 2026) | Multi-omic atlas of human autonomic and sensory ganglia implicates cell types in peripheral neuropathies | bioRxiv | 41000833 | 10.1101/2025.09.18.677119 | snMultiome (RNA + ATAC) of human sympathetic ganglia + DRG; SG vs DRG differences; GRNs; disease-gene mapping incl. HSAN types I to IX. SPARC 10.26275/oqxx-si9h; Shiny portal. |
| 23 | Zhao, 2025 (preprint; v2 lists Boyer K first) | Single nucleus multiomic atlas of human DRG... non-neuronal cell types to pain | bioRxiv | 41040282 | 10.1101/2025.09.21.677655 | Postmortem hDRG snRNA + snATAC; 9 major types; peak-gene links; pain genes in non-neuronal cells. SPARC 10.26275/sojr-qhvm. |
| 27 | He, 2026 (preprint) | Multi-omic profiling of human and mouse DRG... gene delivery to nociceptors | bioRxiv | 41959226 | 10.64898/2026.03.05.709931 | Human + mouse DRG snMultiome; enhancer AAV screen; sequence model of nociceptor cis-regulation. |
| 28 | Jankelunas, 2026 (preprint) | Canine DRG cell types by cross-species comparison | bioRxiv | 42039410 | 10.64898/2026.04.14.718477 | Dog DRG mapped onto existing cross-species atlases. Listed to show the harmonised framework is still being extended. |
| 63 | Meriau, 2026 (preprint) | Integrated single-cell atlas of peripheral glia | bioRxiv | none | 10.64898/2026.01.29.702614 | More than 200,000 DRG + sciatic nerve cells integrated across studies and injury models (rodent, from the abstract); 28 types; SGC and Schwann sub-states; progenitor-like cells. Not in PubMed. |

### 2b. Human DRG spatial and disease (diabetes, pain)

| # | First author, year | Title (short) | Journal | PMID | DOI | Covers |
|---|---|---|---|---|---|---|
| 4 | Ray, 2023 | RNA profiling of human DRG... sex differences in neuropathic pain | Brain | 35867896 | 10.1093/brain/awac266 | Bulk RNA-seq of 70 DRG (50 passed QC) from thoracic vertebrectomy patients, pain vs no pain. Raw data location UNVERIFIED (likely UT Dallas dbGaP phs001158). |
| 11 | Shiers, 2025 | Nageotte nodules in human DRG... DPN | Nat Commun | 40325011 | 10.1038/s41467-025-59538-z | 90 organ donors screened; spatial sequencing of DPN DRG; nodules made of SGC and non-myelinating Schwann cells. GEO GSE295206 (processed), dbGaP phs001158 (raw). |
| 12 | Hall, 2022 | Transcriptomic analysis of human sensory neurons in painful DPN | Sci Rep | 35304484 | 10.1038/s41598-022-08100-8 | Bulk RNA-seq of L4/L5 DRG, 5 DPN vs 7 controls; macrophage inflammation, loss of neuronal genes. Accession not checked (UNVERIFIED). |
| 24 | Sankaranarayanan, 2026 (preprint) | Progressive neurodegeneration in human DRG from diabetes to painful neuropathy | bioRxiv | 41648154 | 10.64898/2026.01.16.700028 | Multi-omic DRG data from healthy, diabetic-no-neuropathy and DPN donors; loss of A-fibre subsets; non-neuronal remodelling. |
| 25 | Mazhar, 2026 (preprint) | Satellite microglia-like cells in human DRG... DPN | bioRxiv | 42182418 | 10.64898/2026.05.12.724479 | Large hDRG snRNA-seq + spatial data with and without DPN; TMEM119+/P2RY12+ cells in the satellite envelope; also sural nerve. |
| 26 | Mancilla Moreno, 2025 (preprint) | Spatial transcriptomics of human paravertebral sympathetic chain ganglia... diabetes | bioRxiv | 40667173 | 10.1101/2025.06.24.661373 | Lumbar sympathetic chain spatial data from diabetic vs non-diabetic donors. |
| 15 | Ahlgreen, 2026 | Mapping SGC heterogeneity... DRG | Adv Sci | 41134052 | 10.1002/advs.202511569 | Mostly mouse SGC subtypes mapped in situ; human DRG shows inner and outer SGC layers. |
| 16 | Choi, 2026 | Spatially localized immune metaprograms in human DRG | PLoS One | 42546053 | 10.1371/journal.pone.0354750 | Reanalysis only: GSE189501 immune nuclei, NMF programmes projected onto Xenium GSE273557, checked on GSE168243. |
| 17 | Schwabenland, 2026 | Spatial immune profiling... microglia-like cells in human DRG | Eur J Immunol | 42153599 | 10.1002/eji.70211 | Imaging mass cytometry with 41 markers on human DRG. |
| 18 | Zhang, 2024 | Spatial + snRNA-seq atlas of adult human spinal cord | eLife | 38289829 | 10.7554/eLife.92046 | Spinal cord atlas; includes a DRG Visium series (GEO GSE243073). |
| 19 | Vroman, 2023 | Matrisome expression in murine and human DRG | Front Mol Neurosci | 37664243 | 10.3389/fnmol.2023.1232447 | Integrates bulk, single-cell and spatial data; human DRG Visium in GEO GSE215994 (2 samples). |
| 13 | Toikumo, 2025 | Cell-type-specific genetic architecture of chronic pain in brain and DRG | J Clin Invest | 41055971 | 10.1172/JCI197583 | Chronic pain GWAS × hDRG scRNA-seq; hPEP.TRPV1/A1.2 enriched; full meta-analysis summary statistics only on request. |

### 2c. Human peripheral nerve

| # | First author, year | Title (short) | Journal | PMID | DOI | Covers |
|---|---|---|---|---|---|---|
| 29 | Heming, 2025 | Perineurial hyperplasia and lipid-associated nerve macrophages in human polyneuropathies | Nat Commun | 40849297 | 10.1038/s41467-025-62964-8 | Sural nerve snRNA-seq of 33 PNP patients + 4 controls (365,708 nuclei) + Xenium; CXCL14 perineurial marker; repair Schwann cells up. GEO GSE285983/GSE285984; Zenodo 10.5281/zenodo.14226218; portal pns-atlas.mzhlab.com. |
| 33 | Tavares-Ferreira, 2025 | Peripheral nerves shift toward inflammatory phenotypes in DPN | J Clin Invest | 40828619 | 10.1172/JCI184075 | Bulk transcriptomes of tibial + sural nerves from amputations (mostly diabetic), control sural nerves, and sural spatial data. dbGaP phs001158; portal sensoryomics.shinyapps.io/human-peripheral-nerves. |
| 32 | Morchio, 2026 | Spatially resolved single-cell analysis... human neuromas | Pain | 41553171 | 10.1097/j.pain.0000000000003907 | Trigeminal nerve and neuroma snRNA-seq + spatial data, painful vs non-painful. |
| 30 | Yim, 2022 | Glial diversity in peripheral nerves at single-nuclei resolution | Nat Neurosci | 35115729 | 10.1038/s41593-021-01005-1 | Mouse nerve snRNA-seq atlas; motor-associated Pmp2+ Schwann cells; human ALS nerve check. GEO GSE182098 (mouse). |
| 31 | Stefanovska, 2025 | FGF signalling... chondrocyte-like nerve fibroblast in aging | Nat Commun | 41238529 | 10.1038/s41467-025-65297-8 | Mouse sciatic nerve aging snRNA-seq; human perineurial fibroblasts in vitro. |
| 34 | Wolbert, 2020 | Heterogeneity of peripheral nerve cells in health and autoimmunity | PNAS | 32295886 | 10.1073/pnas.1912139117 | Rodent nerve scRNA-seq (GEO GSE142541). |
| 35 | Gerber, 2021 | Sciatic nerve ATlas (SNAT) | eLife | 33890853 | 10.7554/eLife.58591 | Mouse sciatic nerve single-cell profiles. |
| 36, 37 | Hur 2011; Guo 2020 | Sural nerve DPN transcriptomes | Brain; Clin Epigenetics | 21926103; 32787975 | 10.1093/brain/awr228; 10.1186/s13148-020-00913-6 | Hur: microarray of 50 sural nerve samples from a 52-week trial. Guo: sural nerve RNA-seq + methylation from T2D patients with DPN. |
| 38 | Svaren, 2019 | Schwann cell transcript biomarkers for hereditary neuropathy skin biopsies | Ann Neurol | 30945774 | 10.1002/ana.25480 | CMT1A skin biopsies with a Nanostring panel, PMP22 normalised to Schwann cell genes. Not genome-wide. |
| 40 | Gummadi, 2026 (preprint) | Small RNA-seq of human sural nerves in DPN | bioRxiv | 41648488 | 10.64898/2026.01.15.699768 | miRNA only (GEO GSE324215, 12 samples). |

### 2d. GTEx nerve, deconvolution, eQTL

| # | First author, year | Title (short) | Journal | PMID | DOI | Covers |
|---|---|---|---|---|---|---|
| 47 | Donovan, 2020 | Cellular deconvolution of GTEx tissues | Nat Commun | 32075962 | 10.1038/s41467-020-14561-0 | 28 GTEx tissues deconvolved with mouse signatures; cell-type-interaction eQTLs for liver and skin. Full text has no mention of nerve. |
| 48 | Kim-Hellmuth, 2020 | Cell type-specific genetic regulation across human tissues | Science | 32913075 | 10.1126/science.aaz8528 | Cell-type-interaction QTLs for 7 broad cell types (xCell-type estimates). |
| 39 | Ray, 2019 | Transcriptome analysis of human tibial nerve... sexually dimorphic | Front Mol Neurosci | 30890918 | 10.3389/fnmol.2019.00037 | 248 GTEx tibial nerve samples; 149 sex-DE genes; no deconvolution. |
| 49 | Parisien, 2017 | Human genetic variability on gene expression in DRG | Cell Rep | 28564610 | 10.1016/j.celrep.2017.05.018 | Human DRG eQTL map (arrays, GEO GSE78150/GSE78148). |
| 62 | Eraslan, 2022 | Single-nucleus cross-tissue reference maps | Science | 35549429 | 10.1126/science.abl4290 | GTEx snRNA-seq of 8 tissues, 16 donors. The abstract does not list the tissues, so whether tibial nerve was included is UNVERIFIED (I believe it was not). |

### 2e. Mouse references (layer 1 alignment)

Usoskin 2015 [41] (GSE59739), Sharma 2020 [42] (GSE139088), Zeisel 2018 [43] (mousebrain.org; CELLxGENE collection 3d79a8b3...), Qi 2024 [44] (DRG genetic toolkit), Renthal 2020 [45] (injury, GSE154659), Mapps 2022 [46] (SGC in sympathetic and sensory ganglia, GSE175421).

## 3. What is already done, and the gap

Already done, and should not be claimed as new:
- A harmonised cross-species DRG/TG atlas with a portal [7], now overtaken by a 126-donor human DRG reference with cross-species integration, spatial data and a portal [21].
- Human SG + DRG multiome with HSAN gene-to-cell-type mapping [22, 23].
- Chronic pain GWAS enrichment in hDRG neuron types [13, 2], and migraine in TG [6].
- A human sural nerve snRNA-seq + Xenium atlas across acquired polyneuropathies, with a portal [29].
- Spatial data on human DRG neurons and SGC neighbourhoods (Visium [3], Xenium [8, 21], MERSCOPE [21], DPN spatial [11, 25]).
- DPN effects in DRG, sympathetic ganglia and nerve [11, 12, 24, 25, 26, 33].

Gaps I could not find filled, as of Nov 2026:
1. A joint human ganglia + nerve reference. The DRG atlases include little nerve tissue, and the nerve atlas [29] includes no ganglia. Nobody has annotated SGCs, ganglion Schwann cells, nerve myelinating/Remak/repair Schwann cells, endo/peri/epineurial fibroblasts and nerve macrophages in one shared nomenclature. Meriau 2026 [63] does DRG + sciatic nerve glia, but in rodents according to its abstract.
2. Inherited neuropathy genes across all classes (CMT1/2/4/X, dHMN, HSAN, from PanelApp/OMIM) scored jointly against neuronal subtypes, SGCs and nerve Schwann/fibroblast/perineurial states. Yang 2025 [22] covers HSAN in ganglia only. This is the strongest angle for this group, given the existing gene-ranking tool.
3. GTEx tibial nerve (v10: 670 RNA-seq samples, 667 with eQTL) deconvolved with a human nerve reference, and cell-state proportions tested against age, sex and death circumstances. I found no published example. Diabetes and BMI live in the controlled GTEx phenotype tables (dbGaP phs000424), so the diabetes analysis needs a dbGaP application.
4. Cell-state-aware tibial nerve eQTL colocalisation (cell-type-interaction eQTL) with carpal tunnel [50], neuropathic pain [51], mono/polyneuropathy [52] and DPN GWAS. Kim-Hellmuth [48] and Donovan [47] did not resolve nerve-specific states.
5. A static gene lookup portal that covers ganglia + nerve + GTEx-derived proportions + genetics in one place. The existing portals are Shiny apps for single datasets.

Risk to novelty: Renthal/Bhuiyan (PRECISION network) and Gereau/Yang are clearly heading this way, and the 126-donor atlas [21] is not yet in a journal. Expect it, and possibly a nerve extension, within 12 months. Building the layer 1 DRG reference from scratch would duplicate [21]. The better plan is to use [21]/[7] as the DRG backbone and put the effort into the nerve side, the GTEx deconvolution, and the inherited neuropathy genetics.

## 4. Dataset inventory

Access: "open" means processed counts can be downloaded without an application. "controlled" means dbGaP or similar.

### Layer 1: single-cell, human

| Dataset | Accession | Technology | Donors / samples | Condition | Access | PMID |
|---|---|---|---|---|---|---|
| Human DRG reference atlas | humanDRG.painseq.com (Seurat object on Dropbox); dbGaP/SPARC promised | sc/snRNA-seq (3', FLEX), soma-seq, Visium | 126 donors (53,596 neuronal cells in portal) | organ donors, cervical/thoracic/lumbar | Seurat object open; raw data pending | 41279342 (preprint) |
| Harmonized DRG/TG atlas | harmonized.painseq.com; GEO GSE255436 (axolotl only) | integration of 31 datasets | 10 new human donors | healthy | portal open; human raw data location UNVERIFIED | 38905344 |
| Nguyen human DRG | GEO GSE168243 (6 samples); CELLxGENE collection 741abcf7 | snRNA-seq | 6 samples | healthy donors | open | 34825887 |
| Jung cross-species | GEO GSE201654 (52), GSE201586 (human, 18) | snRNA-seq | human 7 donors | healthy | open | 36690629 |
| Yu single-soma | GEO GSE249746 (1,136 cells) | Smart-seq2 single soma | n=3 donors in a key figure (full count UNVERIFIED) | healthy | open | 39496796 |
| Human glia DRG + spinal cord | GEO GSE189501 (6) | snRNA-seq | 6 samples | | open | none linked |
| Avraham SGC | GEO GSE169301 (8, rat + human) | scRNA-seq | | | open | 35503034 |
| Xu 2026 two human DRG | GEO GSE310000 (2) | snRNA-seq | 2 | | open | 41501451 |
| Human DRG (itch study) | GEO GSE296776 (9) | RNA-seq | 9 | | open | 40436117 [66] |
| SG + DRG multiome | SPARC 10.26275/oqxx-si9h (SG); 10.26275/sojr-qhvm (DRG) | snMultiome | multiple donors (count not checked) | organ donors | open via SPARC (terms UNVERIFIED) | 41000833; 41040282 |
| Human TG | GEO GSE197289 (36, human + mouse) | snRNA + snATAC | | | open | 35349784 |
| Human sural nerve PNP | GEO GSE285983 (37) | snRNA-seq | 33 PNP + 4 controls | CIDP, vasculitic, diabetic etc. | open | 40849297 |
| C2 DRG, neck pain | not checked | bulk + snRNA + spatial | 22 patients | acute vs chronic pain | UNVERIFIED (likely dbGaP phs001158) | 41032656 |
| UT Dallas human DRG / nerve programme | dbGaP phs001158 | various raw data | many | organ donors, DPN, pain surgery | controlled | 35171654, 40325011, 40828619 |

Macaque/mouse for alignment: macaque DRG GSE165569 (Kupari, PMID 33686078); macaque within GSE201654; mouse GSE59739 (Usoskin), GSE139088 (Sharma), GSE154659 (Renthal injury), GSE175421 (SG/DRG satellite glia), GSE182098 (Yim nerve), GSE142541 (Wolbert nerve).

### Layer 2: spatial, human

| Dataset | Accession | Technology | Samples | Condition | Access | PMID |
|---|---|---|---|---|---|---|
| Yu hDRG Xenium | GEO GSE273557 (4); Dryad 10.5061/dryad.gf1vhhmxq (not resolved) | Xenium | 4 | healthy | open | 39496796 |
| Shiers DPN DRG | GEO GSE295206 (16) | spatial sequencing (platform per paper) | 16 | diabetes / DPN donors | open processed; raw on dbGaP | 40325011 |
| Tavares-Ferreira DRG Visium | sensoryomics.com; dbGaP phs001158 | Visium | 8 donors | healthy | processed via portal; raw controlled | 35171654 |
| Bhuiyan 2025 Xenium + MERSCOPE | pending | Xenium 472 genes, MERSCOPE 300 genes | lumbar DRG | healthy | not yet released | 41279342 |
| Zhang DRG Visium | GEO GSE243073 (7) | Visium | 7 | | open | 38289829 |
| Vroman DRG Visium | GEO GSE215994 (2) | Visium | 2 | | open | 37664243 |
| Heming sural nerve Xenium | GEO GSE285984 (8) | Xenium | 8 | PNP + control | open | 40849297 |
| DPN sural nerve spatial | dbGaP phs001158 | spatial | | DPN | controlled | 40828619 |
| Sympathetic chain spatial | not checked | spatial | | diabetes vs not | UNVERIFIED | 40667173 |

### Layer 3: bulk

| Dataset | Accession | n | Condition | Access | PMID |
|---|---|---|---|---|---|
| GTEx v8 tibial nerve | GTEx portal / dbGaP phs000424 | 619 RNA-seq (200 F / 419 M), 532 eQTL | postmortem donors aged 20 to 70 | expression open; donor diabetes/BMI controlled | |
| GTEx v10 tibial nerve | GTEx portal | 670 RNA-seq (211 F / 459 M), 667 eQTL; 20,101 eGenes | same | same | |
| Ray 2023 DRG | UNVERIFIED (likely phs001158) | 50 DRG after QC | neuropathic pain vs not | controlled? | 35867896 |
| Hall 2022 DPN DRG | UNVERIFIED | 5 DPN + 7 control | DPN | | 35304484 |
| Parisien DRG eQTL | GEO GSE78150 (428 arrays), GSE78148, GSE77968 (214) | 214 DRG | pain phenotypes | expression open; genotypes not checked | 28564610 |
| Sural nerve DPN RNA-seq | GEO GSE148059 (77); methylation GSE148060 | 77 | T2D DPN by HbA1c | open | 32787975 |
| Sural nerve DPN microarray | GEO GSE24290 (35) | 35 | progressive vs non-progressive DPN | open | 21926103 |
| Tibial + sural nerve DPN | dbGaP phs001158 | | DPN amputations vs controls | controlled | 40828619 |
| ENCODE tibial nerve snATAC | GEO GSE333099, GSE333055, GSE333024 | 1 each | ENCODE donors | open | |
| CMT patient nerve/skin | none found genome-wide; only Nanostring skin panels (Svaren 2019) | | CMT1A | | 30945774 |

No public genome-wide CMT/dHMN patient nerve or skin RNA-seq turned up in GEO (searched "Charcot-Marie-Tooth OR hereditary neuropathy", human). Only iPSC-derived and fibroblast datasets exist (e.g. GSE69988, GSE97851, GSE144208). The disease arm of layer 3 will have to use acquired polyneuropathy (GSE285983) and DPN nerves.

### Layer 4: genetics

| Trait | GWAS Catalog ID | Sample | Full summary stats | PMID |
|---|---|---|---|---|
| Carpal tunnel syndrome | GCST007581 | 12,312 cases / 389,344 controls (UKB) | yes | 30833571 |
| Carpal tunnel syndrome | GCST90018813 | 18,855 cases / 461,346 controls + EAS | yes | 34594039 |
| Carpal tunnel (largest) | GCST90104786 | 48,843 / 1,190,837 | no | 35332129 [64] (Skuladottir 2022 Nat Commun) |
| Trigger finger (shared with CTS) | GCST90104907 | 2,908 / 436,579 | yes | 36043126 |
| Neuropathic pain, all-cause | GCST90859373 | 34,278 / 140,134 (UKB) | yes | 40899795 |
| Sciatica; central NP; post-surgical NP; fibromyalgia | GCST90859374 to GCST90859379 | various | yes | 40899795 |
| Diabetic neuropathy | GCST90129436 | 772 / 435,199 | yes | 37219871 |
| T2D with neurological manifestations (MVP, EUR) | GCST90475676 | 59,205 / 375,439 | yes | 39024449 |
| Hereditary and idiopathic peripheral neuropathy (MVP, EUR) | GCST90475844 | 36,894 / 386,060 | yes | 39024449 |
| Neuropathic pain in T2D | GCST003063 | ~961 cases | no | 26629533 |
| DPN in T2D (ACCORD) | GCST008223 | 4,384 / 784 | no | 31127053 |
| Multisite chronic pain | GCST008512 | 387,649 | yes | 31194737 |
| Chronic pain, more than 1.2M (Toikumo) | not in Catalog | 1,235,695 | on request | 40093235 / 41055971 |
| Mono- and polyneuropathies (FinnGen + UKB) | not in Catalog | 983,477 | UNVERIFIED (FinnGen R-releases open) | 41518141 |
| DPN, MVP + UKB + AoU | preprint | | UNVERIFIED | 42326776 |
| CMT1A modifier GWAS | GCST008530 | 471 | no | 30958311 [65] (Tao 2019 J Neuromuscul Dis) |

eQTL: GTEx v10 tibial nerve (667 samples, 20,101 eGenes, 6,590 sGenes); GTEx v8 (532); human DRG eQTL [49].
Mendelian genes: Genomics England PanelApp (hereditary neuropathy panels) and OMIM still need to be pulled. PanelApp has an open API; OMIM needs a licence key.

## 5. Main risks and blockers

1. Scoop risk: the 126-donor atlas [21] and the SG/DRG multiome [22] will probably be published soon, and both have portals. Layer 1 on its own is not publishable as new.
2. Access: much of the US human DRG/nerve data (UT Dallas phs001158, the Bhuiyan 2025 raw data, the Toikumo GWAS) is controlled or "on request". GTEx diabetes status also needs dbGaP phs000424.
3. Disease data for inherited neuropathy is thin. There is no public CMT patient nerve or skin transcriptome; the only disease references are acquired PNP [29] and DPN.
4. Batch effects: there are many chemistries (3', FLEX, soma Smart-seq, Visium spots, Xenium panels). Neuronal nuclei make up a small share of DRG snRNA-seq, so nerve-to-ganglion integration needs careful benchmarking.
5. Deconvolution of tibial nerve depends on having a nerve reference with adipocytes and epineurial cells. GTEx tibial samples include epineurial fat, so check that the reference [29] covers those cells.

## 6. Files

- References: `novelty_refs.ris` has 66 records. ID n in the RIS is the number in square brackets above, and AN holds the PMID (record 63 is a preprint with no PMID). Abstracts were stripped from the RIS.
- The searches were run against live PubMed, Europe PMC, bioRxiv, GEO, the GTEx API v2, the GWAS Catalog API v2, CELLxGENE and DataCite on 2026-11-02.
