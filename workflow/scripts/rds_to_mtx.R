#!/usr/bin/env Rscript
# Write an authors' R object out as matrix market plus CSV so the Python side
# can read it without R. Handles a Seurat object (counts of the RNA assay,
# all layers joined, plus the cell metadata) or a bare sparse/dense count
# matrix saved with saveRDS.
#
# Usage: Rscript workflow/scripts/rds_to_mtx.R <input.rds[.gz]> <outdir>
# Output: <outdir>/matrix.mtx.gz (genes x cells), features.tsv.gz,
#         barcodes.tsv.gz, metadata.csv.gz (Seurat only), object_info.txt
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: rds_to_mtx.R <input.rds> <outdir>")
inp <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# GEO gzips files that are already gzip-compressed RDS, so a second
# decompression layer is tried when the first read fails.
obj <- tryCatch(readRDS(inp), error = function(e) {
  con <- gzcon(gzfile(inp, "rb")); on.exit(close(con)); readRDS(con)
})
info <- c(paste("class", paste(class(obj), collapse = ",")))

if (inherits(obj, "Seurat")) {
  suppressPackageStartupMessages(library(Seurat))
  info <- c(info, paste("seurat_version", as.character(obj@version)),
            paste("assays", paste(Assays(obj), collapse = ",")),
            paste("default_assay", DefaultAssay(obj)))
  assay <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
  a <- obj[[assay]]
  if (inherits(a, "Assay5")) {
    layers <- Layers(a)
    info <- c(info, paste("layers", paste(layers, collapse = ",")))
    if (sum(grepl("^counts", layers)) > 1) obj[[assay]] <- JoinLayers(a)
    counts <- LayerData(obj[[assay]], layer = "counts")
  } else {
    counts <- GetAssayData(a, layer = "counts")
  }
  meta <- obj@meta.data
  meta <- cbind(barcode = rownames(meta), meta)
  for (red in names(obj@reductions)) {
    emb <- Embeddings(obj, red)
    if (ncol(emb) <= 3) {
      colnames(emb) <- paste0(red, "_", seq_len(ncol(emb)))
      meta <- cbind(meta, emb[rownames(obj@meta.data), , drop = FALSE])
    }
  }
  fwrite(meta, file.path(out, "metadata.csv.gz"))
} else {
  counts <- obj
}
if (!inherits(counts, "dgCMatrix")) counts <- as(as(as(counts, "dMatrix"), "generalMatrix"), "CsparseMatrix")
info <- c(info, paste("genes", nrow(counts)), paste("cells", ncol(counts)),
          paste("nonzero", length(counts@x)),
          paste("integer_counts", all(counts@x == round(counts@x))))
writeLines(info, file.path(out, "object_info.txt"))
writeMM(counts, file.path(out, "matrix.mtx"))
system2("gzip", c("-f", file.path(out, "matrix.mtx")))
fwrite(data.table(rownames(counts)), file.path(out, "features.tsv.gz"), col.names = FALSE)
fwrite(data.table(colnames(counts)), file.path(out, "barcodes.tsv.gz"), col.names = FALSE)
cat(paste(info, collapse = "\n"), "\n")
