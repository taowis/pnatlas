# Analysis log, stage 2: joint ganglia and nerve reference

What was found and decided at each step. Commands are in
`docs/runbook/02_reference.md` and, for the environment, in
`docs/runbook/README.md`.

## 2.1 Environment

The conda environment was built for the first time here. `environment.yml`
was extended with the stage 2 tools before the build: scvi-tools (with the
CPU build of PyTorch from conda-forge), jax and jaxlib for scib-metrics,
leidenalg and igraph for clustering, h5py, matplotlib and seaborn. Three
packages are not on conda-forge for this platform and are pinned through pip
inside the environment: harmonypy, scib-metrics and scrublet.

The solve picked scvi-tools 1.4.2 rather than the latest 1.5.2 because the
file pins Python 3.11; the difference does not matter for scVI and scANVI as
used in the plan, so the pin stays. The environment takes 2.3 GB on the
laptop's internal disk. PyTorch sees the Apple GPU (MPS), which is enough for
model development on subsets; full-reference scVI runs go to the server. The
package tests pass in the new environment and the `pnatlas` entry point
resolves. The solved versions are in `docs/runbook/environment_lock.yml`,
exported without the machine-specific prefix line.

A machine config `config/paths.yaml` (not committed) was created with the
data root on the external drive, where 485 GB are free, and the disk budget of
300 GB from the example file.

## 2.2 Reference layer download

Where to compute was settled first (decision in `docs/decisions.md`, stage 2):
integration on the server, Gadi for the CPU-parallel stages later, the laptop
for development. The server had an idle GPU, 183 GB free and only the
sequencer's basecall service running when the download was started.

The download list was built by a script from the live GEO directory listings
and the Pennsieve file browser rather than typed, so the file names and sizes
are what the sources serve. For stage 2 only the human count matrices and the
authors' metadata are taken: 19 files, 11.3 GB. The SPARC 465 ATAC fragments
(17.8 GB) and the mouse files in the trigeminal series are left out; the
cross-species Jung 2023 parent series GSE201654 is not fetched because the
human sub-series GSE201586 holds the same human samples plus the annotated
object. The macaque anchor GSE165569 has no supplementary files on GEO; its
counts are per sample and are fetched when the anchor is used. All downloads
go through `download_files.py`, which records md5 and size in
`data/MANIFEST.tsv` as each file completes and can be rerun after an
interruption.

## 2.3 Reading the datasets

Each dataset is read by one script from a per-dataset block in
`config/reference_inputs.yaml` (format, files, which author columns give
sample, donor and labels), so the decisions about donors and labels are in a
config file rather than spread over notebooks. Three sources are R objects
(the painseq Seurat object, the Jung 2023 annotated object, the trigeminal
raw count RDS); a small R script exports them as matrix market plus CSV and
the Python side never loads R objects directly.

What the raw files turned out to contain:

- Nguyen 2021: six CSV count tables, one per preparation and donor, 1,837
  nuclei in all, no labels in the GEO files.
- Jung 2023 human sub-series: 18 Cell Ranger h5 files (filtered and raw per
  sample) from 7 donors; donor and labels come from the authors' object.
- Yu 2024: 1,136 single somata from 3 donors; the cell names encode donor and
  ganglion level.
- Zhang (GSE189501): 30,701 nuclei from 3 DRG runs.
- SPARC 465: the Cell Ranger ARC filtered matrix holds 36,601 genes and
  139,026 ATAC peaks; the peaks are dropped. The SDS samples sheet shows that
  library lib-1 pools three donors, so 8 libraries give 6 donor-level units
  (one of them a pool) rather than 7 donors; the pool is kept as its own unit
  and is not counted as three donors in any donor-level test.
- Heming 2025: 37 cellbender-filtered h5 files, one per donor, 365,708
  annotated nuclei after keeping the cells in the authors' metadata; two
  label levels, the Milbrandt sciatic nerve label and the finer cluster.

Two problems in the first pass. The server sync script only copied named
config files, so the new inputs file was missing on the server; it now
copies the whole config directory except the machine file. The Heming
build took 33 GB because all 37 samples were concatenated before the
unannotated cells were removed and the memory guard logged a warning at 19
percent free; the builder now subsets each sample to the annotated cells as
it is read. A metadata column with mixed text and numbers also failed to
write to h5ad; object columns are now written as text.

The R export failed at first for two reasons worth recording. The two GEO
`.rds.gz` and `.RDS.gz` files are gzip archives of files that are themselves
gzip-compressed RDS, so `readRDS` on the path stops with "unknown input
format"; the script now opens a second decompression layer when the first
read fails. The painseq object is a Seurat v4-style assay, and SeuratObject
5.5 has removed the old `slot` argument, so the counts are taken with
`layer = "counts"`. The R environment itself (R 4.4.3, Seurat 5.5.1) could
only be built once the Python environment build had finished, because two
conda solves at once lock the shared repodata cache.

What the three R objects are:

- The trigeminal RDS is a bare count matrix of 38,028 nuclei whose cell names
  match the authors' barcode metadata one to one (donor, library, class,
  subtype).
- The Jung 2023 object (Seurat 4.0.2) holds 111,925 human nuclei from the 18
  libraries with one label level of nine cell types. It has no donor column;
  the GEO sample characteristics give the donor for each library (7 donors,
  Hom_1 to Hom_7, 1 to 4 libraries each, FACS or density gradient
  isolation), now in `data/metadata/geo_samples.tsv` and copied into the
  inputs file as a donor map.
- The painseq "humanDRG" object (Seurat 5.1.0) is not the whole-tissue
  dataset the plan assumed. It is the harmonized DRG neuron atlas: 53,596
  neurons with 22 subtype labels, drawn from nine datasets (a `Dataset`
  column: four U19 sites plus Jung 2023, Nguyen 2021, Tavares-Ferreira 2022
  and two Bhuiyan 2024 sets). Most of the new U19 cells have `orig.ident`
  "SeuratProject", so donors cannot be recovered from the object, and about
  6,400 of its cells are the neurons of datasets already read here.

Consequence: the painseq object is used as the neuron nomenclature and label
source (its 22 subtypes become the fine neuron level of
`config/celltype_map.tsv`) and as a label-transfer reference for neurons,
not as a 126-donor backbone. Cells whose `Dataset` is Jung 2023, Nguyen
2021 or Tavares-Ferreira 2022 are excluded before integration to avoid
duplicates, and the remaining U19 and Bhuiyan 2024 neurons enter the
reference with the library as their donor-level unit, flagged as donor
unknown so they are left out of donor-level tests. The DRG non-neuronal
backbone is Jung 2023 (7 donors) with Nguyen 2021 and Zhang. This changes
section 1 of the plan and is recorded as the first amendment.

Of the 111,925 annotated Jung 2023 cells, 110,773 are present in the
filtered Cell Ranger matrices; the remaining 1,152 are cells the authors
kept from the raw matrices and are dropped here, since the filtered
matrices are the input for every other 10x dataset.

## 2.4 QC and the cell type map

The QC thresholds were written into the config before the first run and are
deliberately plain: 200 genes, 5 percent mitochondrial reads for nuclei (20
for the whole-cell somata), scrublet per sample. Five datasets come with the
authors' own filtering and annotation (Heming already carries scDblFinder
calls, all singlets), so for those every annotated cell is kept and the
flags are only recorded; a second filtering pass would throw away cells the
labels depend on and make the label-transfer comparison harder to read.
The three unannotated sets (Nguyen, Yu, Zhang) are filtered.

The cell type map has three levels. `class` (neuron, glia, stromal,
vascular, immune) is what every dataset can reach. `state` is the level used
for genetics and deconvolution: DRG neuron families follow the painseq atlas
(C-NP, C-PEP, A-PEP, Ab-LTMR, A-LTMR, A-Propr, C-Thermo, C-LTMR, plus the
ATF3 injured state), nerve states follow Heming 2025 (myelinating,
non-myelinating, repair and damage Schwann; perineurial; endoneurial and
epineurial fibroblast; endothelial; lymphatic endothelial; mural; the
immune types), and ganglion non-neuronal states use the same names where the
cell type is the same (satellite glia, Schwann, fibroblast, endothelial,
mural, macrophage). `fine` keeps the authors' own subdivisions. Trigeminal
labels are mouse-style (NF, PEP, NP, SST, TRPM8, cLTMR): NP, SST, TRPM8 and
cLTMR map onto atlas families, PEP stays unsplit, and NF is kept as its own
name because the atlas splits A-LTMR differently. Sympathetic neurons have
no DRG counterpart and keep their own state. The coarse labels in the
Heming metadata are transferred from the mouse sciatic atlas and are wrong
for the endothelial clusters (they come out as "adip"), so the authors'
fine clusters are the source of the map and the transferred label is only
kept as a column.

QC results (`results/reference/qc_summary.tsv`). Nuclei datasets have
median mitochondrial fractions well under 1 percent and almost no cells
under 200 genes. The doublet calls are modest except in the Zhang glia set,
where scrublet flags 3,211 of 30,701 nuclei (10 percent) and they are
removed, leaving 27,325. Nguyen loses 300 nuclei to the mitochondrial
threshold and keeps 1,535. Among the kept-all datasets the painseq neurons
(10,731 flagged, median 1.1 percent) and Jung 2023 (3,946 flagged) carry
the most mitochondrial flags; neurons have more mitochondrial transcripts
than other nuclei, the flags stay in obs and are not applied.

One deviation from the thresholds as written. The Yu 2024 single somata
have a median mitochondrial fraction of 24 percent, and the 20 percent
whole-cell threshold would remove 808 of 1,136 cells. These are complete
neuronal somata picked by hand and sequenced with Smart-seq2, the file on
GEO is the authors' final set (the 1,136 cells of the paper), and a high
mitochondrial fraction is expected for whole large neurons rather than a
sign of damage. The dataset was moved to the keep-all group after this was
seen, which is a change made after looking at the data and is recorded as
such; the flags are kept in obs so the effect of the threshold can be
checked later. No other threshold was changed.

### 2.5 Merge and integration

The eight QC'd datasets were merged into one AnnData of raw counts. All of
them use gene symbols, but the annotation sets differ (32,214 symbols in
Yang 2022 to 61,817 in the painseq object) and the union of 88,803 names
is padded with version-specific lncRNA and clone identifiers, so the merge
keeps the 19,420 genes present in every dataset (recorded in
`docs/decisions.md`). The cell type map was applied without a missing
label: every author label in the h5ad files has a row in
`config/celltype_map.tsv`, which is now checked by code
(`src/pnatlas/labels.py`, with tests). Following amendment 1, the painseq
cells whose source is Jung 2023 (3,071), Tavares-Ferreira 2022 (1,956) or
Nguyen 2021 (1,343) were removed, leaving 47,226 atlas neurons; the 44,363
of those with no library name are grouped by their atlas source study for
the batch unit, and their donor is unknown.

The merged reference has 644,681 cells, 87 batches and 65 known donors
(`results/reference/merged_summary.tsv`). By compartment it is 278,973
ganglion cells and 365,708 nerve nuclei; by class, 58,372 neurons, 189,753
glia, 199,633 stromal, 121,517 vascular, 45,410 immune and 29,996
unlabelled cells (Nguyen, Zhang and Yu, which carry no author labels). The
state table (`merged_state_by_dataset.tsv`) shows where the labels come
from: the DRG neuron families come almost entirely from the atlas, with
small trigeminal contributions to C-NP, C-Thermo and C-LTMR; the 6,243 Jung
neurons are "neuron unassigned" until label transfer; nerve states come
only from Heming; satellite glia, fibroblast, macrophage, endothelial,
mural and the two Schwann states are shared by two to four ganglion
datasets, which is what the integration has to line up. The merge peaked
at about 30 GB of memory (one warning from the memory guard, no pause).

The first integration launch stopped at the highly variable gene step:
scanpy's seurat_v3 flavour fits a loess curve with scikit-misc, which was
not in the environment. It was added to `environment.yml`, installed on
both machines, the lock file was refreshed and the job relaunched; nothing
else changed.

First integration run. scVI ran its 30 epochs (about 10 seconds per epoch
on the GPU, 2.2 GB of GPU memory), scANVI its 20, and Harmony followed on
the CPU; the job peaked at 38 GB of host memory, one memory guard warning
at 13 percent free and no pause. The scIB metrics on 50,000 labelled cells
per subset (`results/reference/scib_metrics.tsv`) put scANVI first by the
rule: mean bio conservation 0.586 against 0.511 for Harmony, batch
correction 0.569 against 0.579, within the 0.05 margin. The nerve subset
integrates well (graph connectivity 0.95, NMI 0.82 for scANVI); the
ganglion subset is harder (NMI 0.51, iLISI 0.04), partly because its
labels are confounded with dataset.

Looking at the scANVI predictions showed a problem with one dataset. The
predictions agree with the seed labels for 94 percent of the atlas
neurons, 96 percent of the nerve nuclei, 88 percent of the trigeminal
cells and 75 percent of Jung 2023, but for only 7 percent of the
sympathetic ganglion cells (Yang 2025), where most cells of every seeded
state were predicted as "sympathetic neuron". A nearest-neighbour check
gave the same picture in every embedding and in the uncorrected PCA: the
neighbours of a sympathetic ganglion cell carry the same label mix
whatever its own label is, and within the dataset only 26 percent of a
cell's 15 neighbours share its label, which is what chance gives (the
other datasets are at 83 to 95 percent). The cause is not on our side.
The join of `meta.data.csv` to the matrix is by barcode (verified in the
code), and the authors' nCount_RNA agrees with our per-cell count sums
(correlation 0.9998), so each label sits on the barcode the authors gave
it. But the expression of those cells does not follow the labels: the
cells with TH or DBH counts, with PTPRC, with MPZ or with PECAM1 are
labelled in the proportions of the whole dataset, every author label has
the same marker means, and a k-means clustering of the dataset on its own
PCA (12 clusters with clean T cell, macrophage, endothelial, mural,
fibroblast, satellite glia and myelinating Schwann markers) has an
adjusted Rand index of 0.005 against the author labels. The matrix is
sound; the label column of the metadata file does not belong to the
barcodes it is listed with, most likely a cell-order mix-up when the
authors wrote the table. The same test on the other four labelled
datasets is now a script (`check_label_alignment.py`, results in
`results/reference/label_alignment.tsv`) so the problem would be caught
at build time in future.

Decision: the Yang 2025 author labels are treated as unusable
(`integration.labels_unusable` in `config/config.yaml`; the map rows stay
as a record of the authors' label set). The dataset enters the reference
unlabelled, like Nguyen, Zhang and Yu, and gets its states by label
transfer in the next step; "sympathetic neuron" will have to be identified
from the data (TH, DBH, PRPH), and a request to the authors or SPARC for a
corrected metadata table is worth making. Because the wrong labels had
seeded scANVI and entered the ganglion metrics for 22 percent of the
labelled ganglion cells, the merge, integration and evaluation were rerun
with the dataset unlabelled; the first run's outputs were set aside on
the server (`data/interim/run1`) and are not used. The rerun results
follow.

Label alignment check (`results/reference/label_alignment.tsv`). On their
own PCA, the author labels of the atlas, Jung 2023, Yang 2022 and Heming
2025 give k-means adjusted Rand indices of 0.28 to 0.62 and 73 to 92
percent same-label neighbours (chance 8 to 33 percent); Yang 2025 gives
0.006 and 27 percent against a chance level of 26 percent. The other four
label sets are sound.

Rerun with Yang 2025 unlabelled. The merged reference is unchanged in
size (644,681 cells, 87 batches, 65 known donors); 551,360 cells now seed
scANVI and 82,946 are unlabelled. scVI and scANVI ran the same epochs with
the same validation ELBO to three figures. The metrics
(`results/reference/scib_metrics.tsv`, 50,000 labelled cells per subset):

| subset | embedding | NMI | ARI | label ASW | batch ASW | iLISI | graph conn. | batch | bio | total |
|---|---|---|---|---|---|---|---|---|---|---|
| ganglion | scVI | 0.574 | 0.308 | 0.528 | 0.807 | 0.029 | 0.887 | 0.574 | 0.470 | 0.512 |
| ganglion | scANVI | 0.662 | 0.421 | 0.589 | 0.792 | 0.033 | 0.933 | 0.586 | 0.557 | 0.569 |
| ganglion | Harmony | 0.573 | 0.419 | 0.547 | 0.812 | 0.060 | 0.846 | 0.573 | 0.513 | 0.537 |
| nerve | scVI | 0.766 | 0.546 | 0.579 | 0.805 | 0.125 | 0.940 | 0.623 | 0.630 | 0.628 |
| nerve | scANVI | 0.815 | 0.667 | 0.699 | 0.763 | 0.136 | 0.958 | 0.619 | 0.727 | 0.684 |
| nerve | Harmony | 0.726 | 0.519 | 0.618 | 0.883 | 0.197 | 0.851 | 0.644 | 0.621 | 0.630 |

Applying the rule to the means over the two subsets: scANVI has the
higher bio conservation (0.642 against 0.567 for Harmony) and its batch
correction (0.603) is within 0.05 of Harmony's (0.608), so scANVI is the
reference embedding (`results/reference/integration_choice.tsv`). The
ganglion metrics improved once the mislabelled cells were out (NMI 0.51
to 0.66 for scANVI), which is the effect of the labels, not of the
embedding. Harmony mixes batches a little more (higher iLISI and batch
ASW) at the cost of label structure; scANVI keeps labels together, which
is what the downstream label transfer and the deconvolution signatures
need. iLISI is low for every method on the ganglion subset because most
states there come from one or two datasets, so a cell's neighbours are
from its own study by construction; graph connectivity, which asks only
that each state form a connected neighbourhood, is 0.93 to 0.96 for
scANVI.

scANVI's own predictions agree with the seed labels for 90 percent of
the seeded cells (93 percent atlas, 97 percent nerve, 84 percent
trigeminal, 66 percent Jung 2023, where the authors' coarse non-neuronal
classes split across nerve-style states). For the sympathetic ganglion
cells the predictions now come out as 22,338 satellite glia, 9,762
macrophages, 7,886 non-myelinating Schwann cells and 1,694 endothelial
cells, almost exactly the counts of the authors' label set (22,868;
9,685; 7,113; 1,871), which confirms that the authors' labels are the
right set attached to the wrong barcodes. The predictions are kept as a
column (`scanvi_pred`) but are not the annotation: they cross
compartments (767 Jung neurons predicted as mast cells, Zhang glia spread
over nerve fibroblast states), so the next step transfers labels by
nearest neighbours in the scANVI space within the cell's own
compartment, and assigns the sympathetic neurons and the states of the
unlabelled DRG datasets from that.

Outputs on the server: `data/interim/reference_merged.h5ad` (counts,
2.6 GB), `data/interim/reference_embeddings.h5ad` (obs and the four
embeddings, 651 MB), `data/interim/models/{scvi,scanvi}`; the first run is
kept in `data/interim/run1` until the next step is done and will then be
deleted. The result tables are in `results/reference/`.

### 2.6 Label transfer and the state rule

Every cell now has a state. Author-labelled cells keep their mapped state.
The rest (Nguyen, Zhang, Yu and the sympathetic ganglion dataset, which
carry no usable labels, plus the "neuron unassigned" cells of Jung 2023
and the "immune unassigned" and "vascular unassigned" cells of the
trigeminal set) take the majority state of their 15 nearest reference
cells in the scANVI space, looking only at reference cells of the same
compartment and, where the authors gave a class, of that class. A cell
whose top state holds fewer than half of its neighbours is "<class>
unassigned" (or "unassigned" if no class reaches half either). The rule
is in `config/config.yaml` (`transfer`) and `docs/decisions.md`.

The sympathetic ganglion needed its own rule, because the reference has
no sympathetic neuron seed. Looking at the catecholamine genes in those
52,950 nuclei first: TH is detected (one or more counts) in 1,503 nuclei,
DBH in 1,168, PRPH in 3,513, both TH and DBH in 244, and a marker score
(mean log-normalised TH, DBH, PRPH) is above 1 in 362 nuclei. The scANVI
predictions of 2.5 put about 450 of these nuclei into neuron states, and
those are where the marker-positive nuclei sit (among the 56 nuclei
predicted as injured neuron, 31 have a score above 1; among the 22,338
satellite glia, 37). The authors' own metadata counts 1,030 "Neuronal"
cells. Neurons are rare in this nuclei preparation, and low-count nuclei
miss any one marker, so the rule is anatomical: a sympathetic ganglion
nucleus whose neighbours are neurons is a sympathetic neuron, and the
marker detection in those nuclei is tabulated as the check
(`results/reference/sympathetic_neuron_markers.tsv`).

The reliability of the transfer is measured on the labelled cells
themselves: within each compartment, the batches are split into five
folds and each fold's cells are predicted from the other four, with all
classes as candidates, which is the harder case of the unlabelled
datasets. The held-out cells come from batches the reference has not
seen, so the numbers are what to expect for a new donor.

Results (`results/reference/label_transfer_summary.tsv`,
`final_state_by_dataset.tsv`). The run takes a few minutes on the CPU
with no memory warning. Of the 93,321 cells without a seed state, 91,245
received one and 2,076 (2 percent) stay unassigned. By dataset:

- Jung 2023 neurons (6,243): 5,095 go to the atlas families (C-NP 1,785;
  A-PEP 987; A-LTMR (NF) 754; Ab-LTMR 368; C-PEP 342; PEP 263; A-Propr
  253; A-LTMR 224; the rest small) and 1,148 stay "neuron unassigned",
  with a median vote of 0.73. These nuclei sit between families more
  often than the atlas cells do, and the two trigeminal neuron states
  (A-LTMR (NF), PEP) compete with the atlas families for them.
- Nguyen 2021 (1,535): 1,422 neurons across the families, a few glia and
  immune cells, 87 unassigned.
- Yu 2024 single somata (1,136): 1,123 placed in the neuron families
  (A-PEP 331, Ab-LTMR 284, C-NP 221, A-Propr 124, C-PEP 95, injured 29,
  A-LTMR 19, C-LTMR 9, C-Thermo 8), 12 "neuron unassigned". The
  Smart-seq2 somata land inside the nuclei-based atlas families, which
  is the integration doing its job.
- Zhang DRG glia (27,325): fibroblast 11,185; satellite glia 7,867;
  endothelial 2,211; myelinating Schwann 1,714; macrophage 1,648;
  adipocyte 706; mural 657; lymphocyte 505; 119 neurons; 709 unassigned
  at some level.
- Trigeminal unassigned cells (4,132): the vascular cells split into
  1,857 endothelial and 471 mural, the immune cells into 1,779
  macrophages and 25 lymphocytes, none unassigned.
- Sympathetic ganglion (52,950): satellite glia 30,670; macrophage 9,937;
  fibroblast 6,481; endothelial 1,804; myelinating Schwann 1,414; mural
  1,084; lymphocyte 1,065; sympathetic neuron 322; adipocyte 51; 119
  unassigned. The 322 sympathetic neurons come from all six donor units.
  Of them, 196 have at least one of TH, DBH and PRPH detected and 110 at
  least two; TH is detected in 38 percent of them against 2.3 percent of
  the satellite glia and 2.1 percent of the macrophages of the same
  nuclei (`sympathetic_neuron_markers.tsv`). The authors' metadata
  counts 1,030 "Neuronal" cells; 322 is what the embedding supports.
  Sixteen of the 65 fully unassigned nuclei carry two or three of the
  markers and may be neurons with mixed neighbourhoods.

Held-out check (`transfer_holdout.tsv`). Predicting held-out batches
from the rest, 98.4 percent of ganglion cells and 99.4 percent of nerve
nuclei reach the vote threshold, and of those 95.2 and 97.8 percent get
their own state back. The atlas neuron families are recovered at 89 to
94 percent of assigned cells, with the confusion between neighbouring
families (A-PEP with Ab-LTMR, C-PEP with A-PEP); the nerve states at 93
to 99 percent except the graded Schwann states (damage Schwann 66
percent, mostly called non-myelinating; repair Schwann 88 percent).
Two ganglion states do not survive the check and are flagged for the
downstream steps rather than relabelled: the Jung 2023 "Adipocyte"
nuclei (659) are called fibroblast when their batch is held out (0 of
659 recovered), so the ganglion adipocyte calls (Jung 659, Zhang 706,
sympathetic ganglion 51) are a fibroblast-like population whose name
rests on one dataset's label; and the trigeminal non-myelinating Schwann
nuclei (702) are called satellite glia (11 percent recovered), so within
ganglia this state is not separable from satellite glia in the embedding
and the two should be read together there. Neither affects the nerve
compartment, where the deconvolution reference comes from (adipocyte
99.6 percent, non-myelinating Schwann 98.4 percent recovered).

State rule (`state_rule.tsv`). All 31 final states have at least 50
cells from at least 3 known donors, so no state has to be dropped or
merged at this point. The DRG neuron families now have 12 to 22 known
donors each (Jung, Nguyen, Yu, Zhang and the trigeminal set together),
on top of the atlas cells whose donors are unknown; the sympathetic
neuron state has 322 cells from 6 donor units (one of them a pool of
three donors); the smallest state is granulocyte with 365 nuclei from 24
nerve donors.

Outputs: the per-cell annotation is `data/interim/reference_annotation.tsv.gz`
on the server and in the local data root (one row per cell of the
merged reference with the seed state, the transferred state and its vote,
the final state and class, and the label source); the summary tables are
in `results/reference/`. The first integration run of 2.5
(`data/interim/run1` on the server, 681 MB) was deleted after this step;
everything in it is superseded by the rerun and reproducible from the
runbook.
