#!/bin/env/Rscript
###############################################################################
#
# Download GSE117847
#
###############################################################################
source("R/download/common.R")

acc <- snakemake@wildcards[["acc"]]
cache_dir <- dirname(snakemake@output[[1]])

download_geo(acc, cache_dir, snakemake@output[[2]], snakemake@output[[3]], snakemake@output[[4]])
