#!/bin/env/Rscript
###############################################################################
#
# GSE118900
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
supp_file <- file.path(cache_dir, acc, 'GSE118900_MM.scrna-seq.tpm.pass.txt.gz')

if (!file.exists(supp_file)) {
  getGEOSuppFiles(acc, baseDir = cache_dir, filter_regex = 'tpm.pass')
}

# load expr data
#    GeneID.Chr.Start IgM.MGUS1_C01 IgM.MGUS1_C04
# 1 AADACL3|chr1|12776118             0             0
# 2 AADACL4|chr1|12704566             0             0
# 3   ABCA4|chr1|94458394             0             0
expr_dat <- read.delim(gzfile(supp_file))

# separate out gene metadata
fdata <- as.data.frame(str_split(expr_dat[, 1], "\\|", simplify = TRUE))
colnames(fdata) <- c("symbol", "chr", "pos")

# drop combined gene id/position column and add symbol column
expr_dat <- expr_dat[, -1]

pdata <- pData(eset)

# normalize sample order
ind <- match(colnames(expr_dat), make.names(pdata$title))
pdata <- pdata[ind, ]

# sanity check
if (!all(colnames(expr_dat) == make.names(pdata$title))) {
  stop("Sample ID mismatch!")
}

# switch to geo accessions for column names
colnames(expr_dat) <- pdata$geo_accession

expr_dat <- expr_dat %>%
  add_column(symbol = fdata$symbol, .before = 1)

write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
