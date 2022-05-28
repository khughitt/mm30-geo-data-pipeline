#!/bin/env/Rscript
###############################################################################
#
# GSE83503
# 
# Note: GSE83503 was performed on an Affymetrix Human Exon 1.0 ST Array, with multiple
# probes for each exon.
#
# The result of this is that >95% of the probes map to multiple genes, and thus simply
# discarding multi-mapped probes will not be helpful. 
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE) %>%
  column_to_rownames("feature")

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# exclude control spots (4130 / 22011 probes)
mask <- fdata$SPOT_ID != 'control'
dat <- dat[mask, ]
fdata <- fdata[mask, ]

# size factor normalization
dat <- sweep(dat, 2, colSums(dat), '/') * 1E6

# exclude any probes with zero variance (uninformative)
dat <- dat[apply(dat, 1, var) > 0, ]

# columns to include
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, patient_died = `death:ch1`, pfs_event = `relapse:ch1`) %>%
  mutate(pfs_event = ifelse(pfs_event == 0, 0, 1))

sample_metadata$disease_stage <- ifelse(sample_metadata$pfs_event == 1, 'RRMM', 'MM')
sample_metadata$platform_type <- 'Microarray'
sample_metadata$sample_type <- "Patient"

# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fdata$gene_assignment, '//', simplify = TRUE)
gene_symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
gene_symbols <- apply(gene_symbols, 1, function(x) {
  str_trim(x)[x != '']
})
gene_symbols <- unlist(lapply(gene_symbols, paste, collapse = ' // '))

# get expression data and add gene symbol column
expr_dat <- dat %>%
  add_column(symbol = gene_symbols, .before = 1) %>%
  filter(symbol != '')

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# split multi-mapped symbols
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
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
