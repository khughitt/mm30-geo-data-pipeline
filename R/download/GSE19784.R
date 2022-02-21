#!/bin/env/Rscript
###############################################################################
#
# Download GSE19784
#
# Gene expression profiling for molecular classification of multiple myeloma in newly
# diagnosed patients
#
# Broyl et al. (2010)
#
# Clinical trial: http://www.hovon.nl/studies/studies-per-ziektebeeld/mm.html?action=showstudie&studie_id=5&categorie_id=3
#
###############################################################################
source("R/download/common.R")

acc <- snakemake@wildcards[["acc"]]
cache_dir <- dirname(snakemake@output[[1]])
pkg_dir <- dirname(snakemake@output[[2]])

download_geo(acc, cache_dir, pkg_dir)
