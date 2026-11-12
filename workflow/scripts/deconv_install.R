#!/usr/bin/env Rscript
# Install the two deconvolution packages of stage 4 from GitHub at fixed
# commits into the active R environment (pnatlas-r), on top of the
# dependencies listed in environment-r.yml. Prints the installed versions.
# Usage: Rscript workflow/scripts/deconv_install.R
pins <- c(BayesPrism = "Danko-Lab/BayesPrism@v2.2.2", MuSiC = "xuranw/MuSiC@v1.0.0")
for (pkg in names(pins)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(pins[[pkg]], upgrade = "never", dependencies = NA)
  }
  cat(pkg, as.character(packageVersion(pkg)), pins[[pkg]], "\n")
}
