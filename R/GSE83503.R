#!/bin/env/Rscript
#
# Key transcription factors altered in multiple myeloma patients revealed by logic
# programming approach combining gene expression pro ling and regulatory networks
#
# Miannay et al. (2016)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE83503'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo', accession)

raw_data_dir <- file.path(base_dir, 'raw')
processed_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data;
# result is a list with a single entry containing an ExpressionSet instance
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# exclude control spots (4130 / 22011 samples)
mask <- fData(eset)$SPOT_ID != 'control'
eset <- eset[mask, ]

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var) > 0, ]

# perform size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# get relevant sample metadata
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, patient_died = `death:ch1`, pfs_event = `relapse:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# get gene symbols associated with each probe; gene symbols are stored at every
# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fData(eset)$gene_assignment, '//', simplify = TRUE)
symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
symbols <- apply(symbols, 1, function(x) {
  str_trim(x)[x != '']
})
symbols <- unlist(lapply(symbols, paste, collapse = ' // '))

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(symbol = symbols, .after = 1)

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  filter(symbol != '') %>%
  select(-probe_id) %>%
  separate_rows(symbol, sep = " ?//+ ?") %>%
  group_by(symbol) %>%
  summarize_all(median)

# store cleaned expression data and metadata
expr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr.feather', accession))
expr_outfile_nr <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

write_feather(expr_dat, expr_outfile)
write_feather(expr_dat_nr, expr_outfile_nr)
write_tsv(sample_metadata, mdat_outfile)

sessionInfo()
