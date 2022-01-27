#!/bin/env/Rscript
###############################################################################
#
# GEO Dataset downloader - general version
# Keith Hughitt
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
library(iodat)
library(yaml)

# work-around for loading large gzip-compressed csv files
# https://github.com/r-lib/vroom/issues/361
if (snakemake@params[["accession"]] == "GSE31161") {
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
}

# GEO directory to store raw and processed data
cache_dir <- dirname(snakemake@output[[1]])

# download & GEO data
eset <- getGEO(snakemake@params[["accession"]], destdir = cache_dir)[[1]]

expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
  rownames_to_column("feature")

pdata <- pData(eset)
fdata <- fData(eset)

# load workflow-/DAG-level metadata
fname <- sprintf("%s.yml", snakemake@params[["accession"]])
dag_mdata <- normalizePath(file.path("metadata", fname))

# node-level metadata
node_mdata <- list(name = "Unprocessed Data")

# generate data package and write to disk
pkg_dir <- dirname(snakemake@output[[2]])
#setwd(pkg_dir)

resources <- list(
  "data" = expr_dat,
  "row-metadata" = fdata,
  "column-metadata" = pdata
)

pkgr <- Packager$new()

pkgr$build_package(resources, node_metadata = node_mdata, dag_metadata = dag_mdata,
                   pkg_dir = pkg_dir)
