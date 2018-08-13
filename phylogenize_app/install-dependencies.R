### Install necessary dependencies from CRAN and BioConductor

# Recommended R version: 3.4.3 "Kite-Eating Tree"

check.cran.install <- function(pkgs) {
  for (pkg in pkgs) {
    if (!(require(pkg, character.only = TRUE))) { install.packages(pkg) }
  }
}

check.bioconductor.install <- function(pkgs) {
  for (pkg in pkgs) {
    if (!(require(pkg, character.only = TRUE))) { biocLite(pkg) }
  }
}

# CRAN

check.cran.install(c("data.table",
    "phylolm",
    "rentrez",
    "pbapply",
    "XML",
    "tidyverse",
    "phytools",
    "gridExtra",
    "igraph",
    "MASS",
    "nlme",
    "ape",
    "geiger",
    "mvMORPH",
    "fitdistrplus",
    "truncnorm",
    "seqinr",
    "kableExtra",
    "scales",
    "xml2",
    "ggplot2",
    "cowplot",
    "MCMCpack"))

# BioConductor
source("https://bioconductor.org/biocLite.R")
check.bioconductor.install(c("qvalue","ggtree","phyloseq"))

message("Complete")
