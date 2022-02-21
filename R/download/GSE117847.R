#!/bin/env/Rscript
###############################################################################
#
# Download GSE117847
#
# A retained transcriptomic profile characterizes the CD138+ cells in the progression
# from smoldering to active multiple myeloma
#
# Storti et al. (2019)
#
###############################################################################
source("R/download/common.R")

acc <- snakemake@wildcards[["acc"]]
cache_dir <- dirname(snakemake@output[[1]])
pkg_dir <- dirname(snakemake@output[[2]])

download_geo(acc, cache_dir, pkg_dir)
