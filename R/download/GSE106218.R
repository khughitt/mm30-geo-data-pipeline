#!/bin/env/Rscript
###############################################################################
#
# Download GSE106218 
#
# Single cell RNA sequencing of multiple myeloma I
#
# Ryu et al. (2019)
#
###############################################################################
library(GEOquery)
library(tidyverse)
library(eco)

acc <- snakemake@wildcards[["acc"]]
cache_dir <- dirname(snakemake@output[[1]])
pkg_dir <- dirname(snakemake@output[[2]])

eset <- getGEO(acc, destdir = cache_dir)[[1]]

# expression data is missing from getGEO() query result and must be downloaded
# separately
supp_file1 <- file.path(cache_dir, acc, 'GSE106218_GEO_processed_MM_raw_TPM_matrix.txt.gz')

if (!file.exists(supp_file1)) {
  getGEOSuppFiles(acc, baseDir = cache_dir)
}

# load TPM counts
#         MM02_38 MM02_48 MM02_50
# 5S_rRNA       0       0       0
# 7SK           0       0       0
# A1BG          0       0       0
expr_dat <- read.delim(gzfile(supp_file1), row.names = 1) %>%
  rownames_to_column("feature")

pdata <- pData(eset)

# save a placeholder "row metadata" dataframe with just the gene identifiers (no
# metadata is available via fData)
fdata <- data.frame("feature"=expr_dat$feature)

# load workflow-/DAG-level metadata
fname <- sprintf("%s.yml", acc)
dag_mdata <- normalizePath(file.path("metadata", fname))

# node-level metadata
node_mdata <- list(processing = "Unprocessed Data")

# generate data package and write to disk
resources <- list(
  "data" = expr_dat,
  "row-metadata" = fdata,
  "column-metadata" = pdata
)

pkgr <- Packager$new()

pkgr$build_package(resources, node_metadata = node_mdata, dag_metadata = dag_mdata,
                   pkg_dir = pkg_dir)
