#!/bin/env/Rscript
#
#	Molecular prognosis in multiple myeloma (2008)
#
# Decaux et al. (2007)
#
# Note: Survival data provided by corresponding author (Stéphane Minvielle) via email
# on June 6, 2019.
#
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE7039'

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
# for GSE7039, data includes two separate esets with a small number of overlapping
# probes; Also note that the data corresponds to two platforms, only one of which
# has a corresponding .annot.gz file.

# The GEO project includes two replicates ("_A" and "_B") for each patient.
esets <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)

# load additional survival metadata provided by author
survival_dat <- read_csv('/data/human/decaux2008/MM survival time GSE7039.csv')

# combine samples from separate ExpressionSets
# survival units: days
sample_metadata <- pData(esets[[1]]) %>%
  mutate(patient_id = sub('_A', '', title)) %>%
  select(geo_accession, platform_id, patient_id) %>%
  inner_join(survival_dat, by = 'patient_id') %>%
  add_column(geo_accession2 = pData(esets[[2]])$geo_accession, .after = 1) %>%
  mutate(patient_died = ifelse(deceased == 'yes', 1, 0)) %>%
  rename(os_time = follow_up_days) %>%
  select(-deceased)

#all(sample_metadata$patient_id == sub('_B', '', pData(esets[[2]])$title))
# [1] TRUE

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'CD138+'
sample_metadata$tissue <- 'Bone Marrow'
sample_metadata$sample_type <- 'Patient'

# adjust for size separately and then combine
e1 <- exprs(esets[[1]])
e2 <- exprs(esets[[2]])

e1 <- sweep(e1, 2, colSums(e1), '/') * 1E6
e2 <- sweep(e2, 2, colSums(e2), '/') * 1E6

mask1 <- !startsWith(rownames(e1), 'AFFX-')
mask2 <- !startsWith(rownames(e2), 'AFFX-')

e1 <- e1[mask1, ]
e2 <- e2[mask2, ]

# no shared probe ids
# length(intersect(rownames(e1), rownames(e2)))
# [1] 0
expr_dat <- rbind(e1, e2)

# get gene symbols
gene_symbols <- c(fData(esets[[1]])[, 'Gene symbol'][mask1],
                  fData(esets[[2]])[, 'Gene Symbol'][mask2])

# drop ambiguous / non-gene fields
# *Multi Hs   *Genomic sequence *Repeats containing   *Seq not verified                ESTs
#             644                 577                 567                 246                 162
mask <- !startsWith(gene_symbols, '*')

#table(mask)
# mask
# FALSE  TRUE
#  2156 10187

expr_dat <- expr_dat[mask, ]
gene_symbols <- gene_symbols[mask]

# get expression data and add gene symbol column
expr_dat <- expr_dat %>%
  as.data.frame() %>%
  add_column(symbol = gene_symbols, .before = 1) %>%
  filter(symbol != '')

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
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
expr_nr_outfile <- sprintf('%s_gene_expr_nr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_feather(expr_dat_nr, file.path(processed_data_dir, expr_nr_outfile))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))
