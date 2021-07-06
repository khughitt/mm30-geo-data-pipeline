#!/bin/env/Rscript
#
# A retained transcriptomic profile characterizes the CD138+ cells in the progression
# from smoldering to active multiple myeloma
#
# Storti et al. (2019)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(arrow)

# GEO accession
accession <- 'GSE117847'

# directory to store raw and processed data
raw_data_dir <- file.path('/data/raw', accession)
processed_data_dir <- sub('raw', 'clean', raw_data_dir)

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# columns to include (GSE31161)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, mm_stage = `patient diagnosis:ch1`) %>%
  mutate(mm_stage = recode(mm_stage, `active MM` = 'MM')) %>%
  mutate(disease_stage = recode(mm_stage, `progressed SMM` = 'SMM', `non-progressed SMM` = 'SMM'))

# add cell type and disease (same for all samples)
sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

# map from ensgenes to gene symbols
ensgenes <- fData(eset)$SPOT_ID
gene_symbols <- grch38$symbol[match(ensgenes, grch38$ensgene)]

ind <- which(is.na(gene_symbols))
unmapped_ensgenes <- ensgenes[ind]

gene_symbols[ind] <- grch37$symbol[match(unmapped_ensgenes, grch37$ensgene)]

#table(is.na(gene_symbols))
# 
# FALSE  TRUE 
# 53531     7 

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
  add_column(symbol = gene_symbols, .before = 1) %>%
  filter(symbol != '')

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  group_by(symbol) %>%
  summarize_all(median)

print("Final dimensions:")
print(paste0("- Num rows: ", nrow(expr_dat_nr)))
print(paste0("- Num cols: ", ncol(expr_dat_nr)))

# store cleaned expression data and metadata
expr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr.feather', accession))
expr_outfile_nr <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

write_feather(expr_dat, expr_outfile)
write_feather(expr_dat_nr, expr_outfile_nr)
write_tsv(sample_metadata, mdat_outfile)
