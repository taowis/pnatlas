#!/usr/bin/env Rscript
# Deconvolve a bulk count matrix (genes x samples, from deconv_prepare.py)
# with the nerve reference of the same directory, by BayesPrism or MuSiC,
# and write the estimated state proportions (samples x states).
# BayesPrism: the reference cells are cleaned of ribosomal, mitochondrial,
# sex chromosome and MALAT1 genes as its authors advise, then
# protein-coding genes are kept; the cell type and cell state labels are
# both the state, and the final theta is written. MuSiC: the reference
# cells go in as a SingleCellExperiment with the state as the cluster and
# the donor as the sample; the weighted estimate is written.
# Usage: Rscript workflow/scripts/deconv_run.R --prep-dir <dir> --bulk <counts.tsv.gz> \
#          --method bayesprism|music --out <proportions.tsv> [--cores N]
suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--prep-dir", type = "character"),
  make_option("--bulk", type = "character"),
  make_option("--method", type = "character"),
  make_option("--out", type = "character"),
  make_option("--cores", type = "integer", default = 8)
)))
genes <- readLines(file.path(opt$`prep-dir`, "genes.txt"))
cells <- read.delim(file.path(opt$`prep-dir`, "reference_cells.tsv"), stringsAsFactors = FALSE)
ref <- readMM(file.path(opt$`prep-dir`, "reference_counts.mtx"))   # cells x genes
colnames(ref) <- genes
rownames(ref) <- cells$cell
bulk <- as.matrix(read.delim(opt$bulk, row.names = 1, check.names = FALSE))  # genes x samples
bulk <- bulk[genes, , drop = FALSE]
cat("reference", nrow(ref), "cells,", ncol(ref), "genes; bulk", ncol(bulk), "samples\n")

if (opt$method == "bayesprism") {
  suppressPackageStartupMessages(library(BayesPrism))
  sc <- as.matrix(ref)
  sc <- cleanup.genes(input = sc, input.type = "count.matrix", species = "hs",
                      gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), exp.cells = 5)
  sc <- select.gene.type(sc, gene.type = "protein_coding")
  prism <- new.prism(reference = sc, mixture = t(bulk), input.type = "count.matrix",
                     cell.type.labels = cells$state, cell.state.labels = cells$state, key = NULL,
                     outlier.cut = 0.01, outlier.genes.cut = 0.001)
  bp <- run.prism(prism = prism, n.cores = opt$cores)
  theta <- get.fraction(bp = bp, which.theta = "final", state.or.type = "type")
} else if (opt$method == "music") {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(MuSiC)
  })
  sce <- SingleCellExperiment(assays = list(counts = t(ref)),
                              colData = data.frame(state = cells$state, donor = cells$donor, row.names = cells$cell))
  est <- music_prop(bulk.mtx = bulk, sc.sce = sce, clusters = "state", samples = "donor", verbose = FALSE)
  theta <- est$Est.prop.weighted
} else {
  stop("unknown method ", opt$method)
}
theta <- theta[, sort(colnames(theta)), drop = FALSE]
write.table(data.frame(sample = rownames(theta), theta, check.names = FALSE), opt$out,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("written", opt$out, nrow(theta), "samples x", ncol(theta), "states\n")
