# Stage 7: portal, data archive and public copy

## 7.1 Static gene lookup site

The portal is a static site, as decided at the start of the project: no
server, no database, every number precomputed from the committed result
tables, so it can be hosted on GitHub Pages indefinitely and the paper's
data availability rests on the archive DOI rather than on a service.
The builder reads the gene by state feature matrix of stage 6 (which
already holds, per gene, the specificity and log CPM in every reference
state, the compartment sums, the correlations with the nerve state
fractions in GTEx and the colocalisation flags), the gene classes, the
PNPrank truth set, the colocalisation pairs, the state tables of stages
2 to 4 and the enrichment of stage 5, and writes one JSON shard per
two-letter symbol prefix (557 shards for the 17,154 genes of the
universe), a one megabyte index for the search box, a state table and
the site metadata. The page offers three views: a gene (states ordered
by specificity with bars, the GTEx correlations, colocalisation tests,
class and truth-set membership, links to Ensembl and the GTEx portal),
the 31 states (cells, donors, datasets, spatial mapping per compartment,
GTEx mean fractions by both methods, enriched gene classes) and an
about page with the downloads. Deep links of the form `#gene/MPZ` work.
The built site is not committed; the pages workflow rebuilds it from
the tables on every push to `main` that touches them, and a test builds
it in the test suite and reads MPZ back (myelinating Schwann first,
CMT tier 1).

Checked in a browser: the MPZ page lists myelinating Schwann at
specificity 0.608 and log CPM 8.34, then non-myelinating, damage and
repair Schwann; its GTEx correlation with the myelinating Schwann
fraction is 0.88 and it is negative with endothelial, granulocyte and
adipocyte; it is in the truth set (CMT, tier 1, AD, 1993) and on the
three hereditary neuropathy panels. The console was clean. One trap for
the record: the first preview port had a service worker registered by
another local application, which served that application's page under
the site's address; a fresh port showed the site.

## 7.2 Citation and archive metadata

`CITATION.cff` and `.zenodo.json` follow the PNPrank files, with the
public repository name `taowis/pnatlas` and version 0.1.0 (the version
of the first archive, before the manuscript). `results/SNAPSHOT.txt`
names the commit the tables, the site and the archive are built from.

## 7.3 Data archive

What goes to Zenodo is what a reader needs to use the reference or redo
a stage without the raw downloads: the reference object (raw counts
after quality control, the final state and class, how the label was
assigned, the scANVI and scVI embeddings), the cell table, the
state-level specificity and expression tables, the cell2location
abundances and checks per section, the panel transfer tables, the
signature model's cell list and signatures, the deconvolution reference
and held-out mixtures, every result table, and the data manifest and
README. The Bhuiyan 2025 cells are left out of the cell-level object:
their portal states that cell-level derivatives are not to be
redistributed (`config/datasets.yaml`), and the object is exactly such a
derivative. The state-level tables, which are aggregates over tens to
thousands of cells, keep every dataset. The other reference datasets
are GEO or SPARC (CC BY) releases, for which processed objects are the
usual thing to archive. The staging script writes everything into one
directory on the server with an MD5 list; the upload itself needs the
Zenodo login and is done by hand, after the manuscript fixes the
version to cite.

## 7.4 Public copy

As for PNPrank, the public repository is a rewritten copy of `main`
without the manuscript and with no attribution trailer in any commit,
pushed to a fresh repository at publication so that no pull request
reference or private history goes with it. This repository has no such
trailer in its history (checked at stage 6), `private_notes/`,
`CLAUDE.md` and `.claude/` have never been tracked, and the server is
named only as "the server". The build script rehearses the rewrite and
the checks without pushing; the push is a separate decision.

The rehearsal rewrote the 110 commits of `main` in under a second and
reported no attribution line in any message, no private name in any
tree, no manuscript file in any commit and no private path tracked;
the clone was left in a temporary directory and discarded.

The staging ran under the memory guard without a warning and wrote 262
files, 2.4 GB in all, of which the reference object is 2.4 GB:
597,455 of the 644,681 reference cells (the 47,226 Bhuiyan 2025 cells
left out) by 19,420 genes, with the embeddings. Two reruns were needed
before the content was right: the first copy put the deposit README in
the wrong place and left an empty directory, and the second carried a
machine-local README from the server's results directory; the script
now excludes such files. The MD5 list of the final run is committed as
`data/metadata/zenodo_v0.1.0_md5sums.txt`.
