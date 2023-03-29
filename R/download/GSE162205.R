#!/bin/env/Rscript
###############################################################################
#
# Download GSE162205 
#
###############################################################################
library(GEOquery)
library(tidyverse)
library(arrow)

acc <- snakemake@wildcards[["acc"]]
cache_dir <- file.path("/data/raw", acc)

# create cache dir
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, mode = "0755")
}

eset <- getGEO(acc, destdir = cache_dir)[[1]]

# expression data is missing from getGEO() query result and must be downloaded
# separately
supp_file <- file.path(cache_dir, "GSE162205_gene_raw_counts_matrix.txt.gz")

if (!file.exists(supp_file)) {
  getGEOSuppFiles(acc, baseDir = cache_dir, makeDirectory = FALSE)
}

pdata <- pData(eset)

# load raw counts
# ENSG.HGNC_symbol X20160208.Exp.13D.100nM.dBET6.18h.A X20160208.Exp.13D.100nM.dBET6.18h.B
#  ENSG00000223972|DDX11L1      1  0
#  ENSG00000227232|WASH7P       0  0
#  ENSG00000278267|MIR6859-1    0  0
expr_dat <- read.delim(gzfile(supp_file))

# save gene metadata (just symbol & biotype in this case)
fdata <- as.data.frame(str_split(expr_dat[, 1], "\\|", simplify = TRUE))
colnames(fdata) <- c("ensgene", "symbol")

colnames(expr_dat) <- c("symbol", pdata$geo_accession)

write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
