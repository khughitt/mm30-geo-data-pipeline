#!/bin/env/Rscript
#
# Expression data from bone marrow of primary multiple myeloma patients
#
# Zhou et al. (2013)
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE19554'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo', accession)

raw_data_dir <- file.path(base_dir, 'raw')
clean_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, clean_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data;
# result is a list with a single entry containing an ExpressionSet instance
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# in order to normalize downstream comparisons across datasets, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# get relevant sample metadata
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, tumor_stage = `tumor stage:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

gene_symbols <- fData(eset)$`Gene symbol`
gene_symbols <- gsub('///', ' // ', gene_symbols)

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(symbol = gene_symbols, .after = 1)

# store cleaned expression data and metadata
expr_outfile <- file.path(clean_data_dir, sprintf('%s_gene_expr.csv', accession))
mdat_outfile <- file.path(clean_data_dir, sprintf('%s_sample_metadata.csv', accession))

write_csv(expr_dat, expr_outfile)
write_csv(sample_metadata, mdat_outfile)

sessionInfo()
