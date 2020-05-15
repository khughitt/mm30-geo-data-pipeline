#!/bin/env/Rscript
#
#	Gene Expression profiling of Multiple Myeloma
#
# Chauhan et al. (2012)
#
library(GEOquery)
library(tidyverse)
library(feather)

options(stringsAsFactors = FALSE)

# GEO accession
accession <- 'GSE39754'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/1.1', accession)

raw_data_dir <- file.path(base_dir, 'raw')
processed_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data;
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         diagnosis = `diagnosis:ch1`, treatment_response = `treatment_response:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$cell_type <- 'CD138+'

sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$disease[grepl('Healthy', sample_metadata$diagnosis)] <- 'Healthy'

# get gene symbols associated with each probe; gene symbols are stored at every
# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fData(eset)$gene_assign, ' ///? ', simplify = TRUE)
symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
symbols <- apply(symbols, 1, function(x) {
  str_trim(x)[x != '']
})
symbols <- unlist(lapply(symbols, paste, collapse = ' // '))

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
  add_column(symbol = symbols, .before = 1) %>%
  filter(symbol != '')

if(!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?") %>%
  group_by(symbol) %>%
  summarize_all(median)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
expr_outfile_nr <- sprintf('%s_gene_expr_nr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_feather(expr_dat_nr, file.path(processed_data_dir, expr_outfile_nr))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))
