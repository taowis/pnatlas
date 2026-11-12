#!/usr/bin/env Rscript
# Install the two deconvolution packages of stage 4 from GitHub at fixed
# commits (neither repository carries release tags) into the active R
# environment (pnatlas-r), on top of the dependencies listed in
# environment-r.yml. Prints the installed versions.
# Usage: Rscript workflow/scripts/deconv_install.R
pins <- list(
  BayesPrism = list(repo = "Danko-Lab/BayesPrism", ref = "19052052a6f30833b27e2148294459b0ba2f923e", subdir = "BayesPrism"),  # 2.2.3
  MuSiC = list(repo = "xuranw/MuSiC", ref = "f21fe67f5670d5e9fca0ad7550abaae3423eb59c", subdir = NULL)                    # 1.0.0
)
for (pkg in names(pins)) {
  p <- pins[[pkg]]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(p$repo, ref = p$ref, subdir = p$subdir, upgrade = "never", dependencies = NA)
  }
  cat(pkg, as.character(packageVersion(pkg)), p$repo, p$ref, "\n")
}
