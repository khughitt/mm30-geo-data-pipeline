#!/bin/env/Rscript
###############################################################################
#
# Generate non-redundant gene expression matrices
#
# Collapses gene symbols with multiple entries into a single row by taking the
# median expression value within each sample.
#
###############################################################################
library(tidyverse)
library(iodag)

# directory to store raw and processed data
data_dir <- dirname(snakemake@output[[1]])

# load data & metadata
dat <- read_csv(snakemake@input[[1]])

row_mdata <- read_csv(snakemake@input[[2]])
col_mdata <- read_csv(snakemake@input[[3]])

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
dat_nr <- dat %>%
  group_by(symbol) %>%
  summarize_all(median)

# create a new data package, based off the old one
pkg_dir <- dirname(snakemake@output[[1]])
setwd(pkg_dir)

pkgr <- Packager$new()

updates <- list(
  data=list(
    processing="non-redundant"
  )
)

pkg <- pkgr$update_package(snakemake@input[[4]], 
                           updates, "collapse mutli-entries",
                           dat, row_mdata, col_mdata)

pkg %>%
  write_package(pkg_dir)
