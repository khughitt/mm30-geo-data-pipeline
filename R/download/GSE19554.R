#!/bin/env/Rscript
###############################################################################
#
# Download GSE19554
#
###############################################################################
source("R/download/common.R")

acc <- snakemake@wildcards[["acc"]]
cache_dir <- file.path("/data/raw", acc)

download_geo(acc, cache_dir, snakemake@output[[1]], snakemake@output[[2]], snakemake@output[[3]])
