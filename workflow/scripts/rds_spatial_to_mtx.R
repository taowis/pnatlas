#!/usr/bin/env Rscript
# Write an authors' spatial Seurat object out as matrix market plus CSV so
# the Python side can read it without R: the raw counts of the Spatial (or
# RNA) assay, the cell metadata, and the tissue coordinates of every image
# with its scale factors, and the tissue image stored with each image
# (the low-resolution H and E that Seurat keeps) as PNG for the anatomy check.
#
# Usage: Rscript workflow/scripts/rds_spatial_to_mtx.R <input.rds> <outdir>
# Output: <outdir>/matrix.mtx.gz (genes x spots), features.tsv.gz,
#         barcodes.tsv.gz, metadata.csv.gz, coordinates.csv.gz
#         (barcode, image, x, y, plus the image's scale factors), object_info.txt,
#         images/<image>.png
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(png)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: rds_spatial_to_mtx.R <input.rds> <outdir>")
inp <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(inp)
if (!inherits(obj, "Seurat")) stop("not a Seurat object")
info <- c(paste("seurat_version", as.character(obj@version)),
          paste("assays", paste(Assays(obj), collapse = ",")),
          paste("default_assay", DefaultAssay(obj)),
          paste("images", paste(Images(obj), collapse = ",")))
assay <- if ("Spatial" %in% Assays(obj)) "Spatial" else if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
a <- obj[[assay]]
if (inherits(a, "Assay5")) {
  layers <- Layers(a)
  info <- c(info, paste("layers", paste(layers, collapse = ",")))
  if (sum(grepl("^counts", layers)) > 1) obj[[assay]] <- JoinLayers(a)
  counts <- LayerData(obj[[assay]], layer = "counts")
} else {
  counts <- GetAssayData(a, layer = "counts")
}
if (!inherits(counts, "dgCMatrix")) counts <- as(as(as(counts, "dMatrix"), "generalMatrix"), "CsparseMatrix")

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

coords <- list()
for (img in Images(obj)) {
  tc <- GetTissueCoordinates(obj, image = img)
  # Seurat 4 names the columns imagerow/imagecol, Seurat 5 x/y with a cell column
  if (all(c("imagerow", "imagecol") %in% colnames(tc))) {
    d <- data.table(barcode = rownames(tc), image = img, x = tc$imagecol, y = tc$imagerow)
  } else {
    bc <- if ("cell" %in% colnames(tc)) tc$cell else rownames(tc)
    d <- data.table(barcode = bc, image = img, x = tc$x, y = tc$y)
  }
  sf <- tryCatch(ScaleFactors(obj[[img]]), error = function(e) NULL)
  if (!is.null(sf)) {
    d[, spot := sf$spot]; d[, fiducial := sf$fiducial]; d[, hires := sf$hires]; d[, lowres := sf$lowres]
  }
  coords[[img]] <- d
  arr <- tryCatch(GetImage(obj[[img]], mode = "raw"), error = function(e) NULL)
  if (is.array(arr) && length(dim(arr)) == 3) {
    dir.create(file.path(out, "images"), showWarnings = FALSE)
    writePNG(arr, file.path(out, "images", paste0(img, ".png")))
    info <- c(info, paste0("image_", img, " ", dim(arr)[2], "x", dim(arr)[1]))
  }
}
coords <- rbindlist(coords, fill = TRUE)
fwrite(coords, file.path(out, "coordinates.csv.gz"))

info <- c(info, paste("genes", nrow(counts)), paste("spots", ncol(counts)),
          paste("nonzero", length(counts@x)),
          paste("integer_counts", all(counts@x == round(counts@x))),
          paste("spots_with_coordinates", nrow(coords)))
writeLines(info, file.path(out, "object_info.txt"))
writeMM(counts, file.path(out, "matrix.mtx"))
system2("gzip", c("-f", file.path(out, "matrix.mtx")))
fwrite(data.table(rownames(counts)), file.path(out, "features.tsv.gz"), col.names = FALSE)
fwrite(data.table(colnames(counts)), file.path(out, "barcodes.tsv.gz"), col.names = FALSE)
cat(paste(info, collapse = "\n"), "\n")
