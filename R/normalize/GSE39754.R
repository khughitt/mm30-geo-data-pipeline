#!/bin/env/Rscript
###############################################################################
#
# GSE39754
#
# Note: GSE39754 was performed on an Affymetrix Human Exon 1.0 ST Array, with multiple
# probes for each exon.
#
# The result of this is that >95% of the probes map to multiple genes, and thus simply
# discarding multi-mapped probes will not be helpful.
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
dat <- read_feather(snakemake@input[[1]]) %>%
  column_to_rownames("feature")

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# exclude any probes with zero variance (uninformative)
dat <- dat[apply(dat, 1, var, na.rm = TRUE) > 0, ]

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id,
         diagnosis = `diagnosis:ch1`, treatment_response = `treatment_response:ch1`)

sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Patient"

sample_metadata$disease_stage <- "MM"
sample_metadata$disease_stage[grepl("Healthy", sample_metadata$diagnosis)] <- "Healthy"

# get gene symbols associated with each probe; gene symbols are stored at every
# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fdata$gene_assignment, " ///? ", simplify = TRUE)
symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
symbols <- apply(symbols, 1, function(x) {
  str_trim(x)[x != ""]
})
symbols <- unlist(lapply(symbols, paste, collapse = " // "))

# get expression data and add gene symbol column
expr_dat <- dat %>%
  add_column(symbol = symbols, .before = 1) %>%
  filter(symbol != "")

if(!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

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
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
