#!/bin/env/Rscript
###############################################################################
#
# GSE7039
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
expr_dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# get gene symbols
gene_symbols <- fdata$`Gene symbol`

# drop ambiguous / non-gene fields
# *Multi Hs   *Genomic sequence *Repeats containing   *Seq not verified                ESTs
#             644                 577                 567                 246                 162
mask <- !startsWith(gene_symbols, "*")

#table(mask)
# mask
# FALSE  TRUE
#  2156 10187

expr_dat <- expr_dat[mask, ]
gene_symbols <- gene_symbols[mask]

# get expression data and add gene symbol column
expr_dat <- expr_dat %>%
  add_column(symbol = gene_symbols, .before = 1) %>%
  filter(symbol != "")

# add platform
pdata$platform_type <- "Microarray"
pdata$sample_type <- "Patient"

# all patients for GSE7039 are indicated as "Newly Diagnosed MM"
pdata$disease_stage <- "MM"

if (!all(colnames(expr_dat)[-1] == pdata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# map from GRCh37 -> 38 where possible, and drop genes which could not be mapped
grch37_mask <- expr_dat$symbol %in% grch37$symbol

gene_symbols <- expr_dat$symbol[grch37_mask]
ensgenes <- grch37$ensgene[match(gene_symbols, grch37$symbol)]

mapped_mask <- ensgenes %in% grch38$ensgene
gene_symbols[mapped_mask] <- grch38$symbol[match(ensgenes[mapped_mask], grch38$ensgene)]

expr_dat$symbol[grch37_mask] <- gene_symbols

mask <- expr_dat$symbol %in% grch38$symbol
expr_dat <- expr_dat[mask, ]

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
