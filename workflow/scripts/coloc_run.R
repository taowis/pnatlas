#!/usr/bin/env Rscript
# Colocalise every region by gene pair prepared by coloc_prepare.py for one
# GWAS with the coloc package (plan, section 6): coloc.abf on the shared
# variants for every pair (GWAS as case-control with the study's case
# fraction when known, otherwise quantitative; eQTL quantitative with sdY
# estimated from the allele frequencies and sample sizes), and, where the
# eQTL gene has a SuSiE credible set in the region, coloc.bf_bf between the
# GWAS single-signal log Bayes factors from coloc.abf and the eQTL SuSiE
# log Bayes factors, taking the largest PP.H4 over signal pairs. The
# reported PP.H4 is the SuSiE-based one where it exists and the abf one
# otherwise. Priors p1, p2, p12 from the config.
# Usage: Rscript workflow/scripts/coloc_run.R --dir <out-dir>/<gwas id> --out <results tsv> --p1 1e-4 --p2 1e-4 --p12 5e-6
suppressPackageStartupMessages({
  library(optparse)
  library(coloc)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--dir", type = "character"),
  make_option("--out", type = "character"),
  make_option("--p1", type = "double", default = 1e-4),
  make_option("--p2", type = "double", default = 1e-4),
  make_option("--p12", type = "double", default = 5e-6)
)))
d <- opt$dir
study <- read.delim(file.path(d, "study.tsv"), stringsAsFactors = FALSE)
pairs <- read.delim(file.path(d, "pairs.tsv"), stringsAsFactors = FALSE)
if (nrow(pairs) == 0) {
  write.table(data.frame(), opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("no pairs\n"); quit(status = 0)
}
is_cc <- tolower(as.character(study$case_control)) %in% c("true", "t", "1")   # written by pandas as True/False
gwas_type <- if (is_cc && !is.na(study$case_count) && !is.na(study$sample_size)) "cc" else "quant"
s_case <- if (gwas_type == "cc") study$case_count / study$sample_size else NA
rows <- list()
for (i in seq_len(nrow(pairs))) {
  rid <- pairs$region[i]; gene <- pairs$gene[i]
  g <- read.delim(file.path(d, paste0("gwas_", rid, ".tsv")), stringsAsFactors = FALSE)
  e <- read.delim(file.path(d, paste0("eqtl_", rid, "_", gene, ".tsv")), stringsAsFactors = FALSE)
  shared <- intersect(g$variant, e$variant)
  g <- g[match(shared, g$variant), ]; e <- e[match(shared, e$variant), ]
  ok <- is.finite(g$beta) & is.finite(g$se) & g$se > 0 & is.finite(e$beta) & is.finite(e$se) & e$se > 0 & is.finite(e$maf) & e$maf > 0 & e$maf < 1
  g <- g[ok, ]; e <- e[ok, ]
  if (nrow(g) < 50) next
  d1 <- list(beta = g$beta, varbeta = g$se^2, snp = g$variant, type = gwas_type, N = study$sample_size)
  if (gwas_type == "cc") d1$s <- s_case
  if (!is.na(study$sample_size) && all(is.finite(g$maf)) && all(g$maf > 0)) d1$MAF <- g$maf
  d2 <- list(beta = e$beta, varbeta = e$se^2, snp = e$variant, type = "quant", MAF = e$maf, N = round(median(e$an) / 2))
  res <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p1 = opt$p1, p2 = opt$p2, p12 = opt$p12)), error = function(err) NULL)
  if (is.null(res)) next
  s <- res$summary
  h4_abf <- unname(s["PP.H4.abf"]); h3_abf <- unname(s["PP.H3.abf"])
  top <- res$results$snp[which.max(res$results$SNP.PP.H4)]
  h4_susie <- NA; h3_susie <- NA; n_signals <- 0
  lbf_path <- file.path(d, paste0("lbf_", rid, "_", gene, ".tsv"))
  if (file.exists(lbf_path)) {
    l <- read.delim(lbf_path, stringsAsFactors = FALSE)
    l <- l[match(g$variant, l$variant), ]
    keep <- !is.na(l$variant)
    if (sum(keep) >= 50) {
      lcols <- grep("^lbf_variable", colnames(l), value = TRUE)
      bf2 <- t(as.matrix(l[keep, lcols]))
      colnames(bf2) <- l$variant[keep]
      bf2 <- bf2[apply(bf2, 1, function(r) any(abs(r) > 1e-6)), , drop = FALSE]
      bf1 <- matrix(res$results$lABF.df1[match(l$variant[keep], res$results$snp)], nrow = 1,
                    dimnames = list("gwas", l$variant[keep]))
      if (nrow(bf2) > 0 && all(is.finite(bf1))) {
        bb <- tryCatch(coloc.bf_bf(bf1, bf2, p1 = opt$p1, p2 = opt$p2, p12 = opt$p12), error = function(err) NULL)
        if (!is.null(bb)) {
          n_signals <- nrow(bf2)
          h4_susie <- max(bb$summary$PP.H4.abf); h3_susie <- bb$summary$PP.H3.abf[which.max(bb$summary$PP.H4.abf)]
        }
      }
    }
  }
  rows[[length(rows) + 1]] <- data.frame(region = rid, gene = gene, n_variants = nrow(g), gwas_type = gwas_type,
                                         pp_h4_abf = h4_abf, pp_h3_abf = h3_abf, top_variant_abf = top,
                                         eqtl_signals = n_signals, pp_h4_susie = h4_susie, pp_h3_susie = h3_susie,
                                         pp_h4 = ifelse(is.na(h4_susie), h4_abf, h4_susie),
                                         method = ifelse(is.na(h4_susie), "coloc.abf", "coloc.bf_bf (SuSiE)"),
                                         stringsAsFactors = FALSE)
}
out <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(out), "pairs written to", opt$out, "\n")
