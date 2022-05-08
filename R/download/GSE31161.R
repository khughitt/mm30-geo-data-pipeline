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

download_geo(acc, cache_dir, snakemake@output[[2]], snakemake@output[[3]], snakemake@output[[4]])
