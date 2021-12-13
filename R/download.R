#!/bin/env/Rscript
###############################################################################
#
# GEO Dataset downloader - general version
# KH (Dec 9, 2021)
#
# Uses GEOquery to download expression data, sample metadata, and gene annotations for a
# single dataset, and stores each as a tabular data package.
#
# For some datasets, the expression data is not stored correctly on GEO, or other issues
# are present, which necessitates a separate script to acquire the data.
#
# Scripts for each such custom data download can be found in the "custom/" folder.
#
###############################################################################
library(GEOquery)
library(tidyverse)
library(iodag)
library(yaml)

# GEO directory to store raw and processed data
cache_dir <- dirname(snakemake@output[[1]])

# download & GEO data
eset <- getGEO(snakemake@params[["accession"]], destdir = cache_dir)[[1]]

expr <- exprs(eset) %>%
  as.data.frame() %>%
  rownames_to_column("feature")

pdata <- pData(eset)
fdata <- fData(eset)

# load dataset metadata
fname <- sprintf("%s.yml", snakemake@params[["accession"]])
recipe <- normalizePath(file.path("metadata", fname))

# generate data package and write to disk
pkg_dir <- dirname(snakemake@output[[2]])
setwd(pkg_dir)

pkgr <- Packager$new()
pkg <- pkgr$build_package(recipe, expr, fdata, pdata)

pkg %>%
  write_package(pkg_dir)
