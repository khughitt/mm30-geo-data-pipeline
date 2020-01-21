#!/bin/env/Rscript
#
# Identification of multiple risk loci and regulatory mechanisms influencing
# susceptibility to multiple myeloma
#
# Went et al. (2018)
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE31161'

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

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# there are three samples in GSE31161 with all missing data..
num_missing <- apply(exprs(eset), 2, function(x) {
  sum(is.na(x))
})

#table(num_missing)
# num_missing
#     0 54675
#  1035     3

eset <- eset[, num_missing == 0]

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# columns to include (GSE31161)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         treatment = `treatment:ch1`, time_of_testing = `time of testing:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'CD138+'

# get gene symbols
symbols <- fData(eset)$`Gene symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(symbol = symbols, .after = 1)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.csv', accession)
mdat_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, mdat_outfile))

sessionInfo()
