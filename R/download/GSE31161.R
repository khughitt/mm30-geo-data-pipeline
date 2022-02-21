#!/bin/env/Rscript
###############################################################################
#
# Download GSE31161
#
# Identification of multiple risk loci and regulatory mechanisms influencing
# susceptibility to multiple myeloma
#
# Went et al. (2018)
#
###############################################################################
source("R/download/common.R")

acc <- snakemake@wildcards[["acc"]]
cache_dir <- dirname(snakemake@output[[1]])
pkg_dir <- dirname(snakemake@output[[2]])

download_geo(acc, cache_dir, pkg_dir)
