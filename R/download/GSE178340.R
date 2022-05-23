#!/bin/env/Rscript
###############################################################################
#
# Download GSE178340 
#
###############################################################################
library(GEOquery)
library(tidyverse)

acc <- snakemake@wildcards[["acc"]]
cache_dir <- file.path("/data/raw/geo", acc)

# create cache dir
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, mode = "0755")
}

eset <- getGEO(acc, destdir = cache_dir)[[1]]

# expression data is missing from getGEO() query result and must be downloaded
# separately
supp_file1 <- file.path(cache_dir, acc, 'GSE178340_fpm.tsv.gz')

if (!file.exists(supp_file1)) {
  getGEOSuppFiles(acc, baseDir = cache_dir)
}

# load FPM counts
#           gene_id SYMBOL ATRA.Cfz_S8_L001
# 1 ENSG00000000003 TSPAN6        0.2802507
# 2 ENSG00000000419   DPM1       89.3999665
# 3 ENSG00000000457  SCYL3       29.7065719
expr_dat <- read.delim(gzfile(supp_file1))
colnames(expr_dat)[1:2] <- c("ensgene", "symbol")

pdata <- pData(eset)
fdata <- expr_dat[, 1:2]

expr_dat <- expr_dat[, -1]

# sanity check
# all(make.names(pdata$title) == colnames(expr_dat)[-1])
# TRUE

# switch to geo accessions for column names
colnames(expr_dat)[-1] <- pdata$geo_accession


write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
