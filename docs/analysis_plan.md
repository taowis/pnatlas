# Analysis plan

Draft. Frozen at the end of stage 1, once data access is known. Later changes
are added as numbered amendments with date and reason; the frozen text is not
edited.

## To be fixed before freezing

1. Datasets per layer, with the inclusion and exclusion rules (`config/datasets.yaml`).
2. Unit of inference: donor. Cells are not replicates; comparisons use
   donor-level pseudobulk or mixed models.
3. Reference construction: integration method, annotation scheme
   (`config/celltype_map.tsv`), and benchmark used to choose between methods.
4. Spatial mapping method and the checks that a mapped state is real
   (marker agreement, held-out genes).
5. Deconvolution method, reference choice, and validation (pseudo-bulk
   mixtures with known proportions; agreement with histology where available).
6. Genetics: gene classes from PanelApp/OMIM with versions, enrichment test,
   eQTL colocalisation method and priors.
7. Stage 6 primary comparison: PNPrank features vs PNPrank + pnatlas
   features, same truth set, cutoff years and splits as the PNPrank
   submission; primary metric tier-1 AUROC; the difference that would count
   as an improvement, and how it is tested.
