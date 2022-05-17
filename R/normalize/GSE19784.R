#!/bin/env/Rscript
###############################################################################
#
# GSE19784
#
###############################################################################
library(tidyverse)
source("R/util/biomart.R")

# output directory to store data packages to
out_dir <- dirname(snakemake@output[[1]])

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE) %>%
  column_to_rownames("feature")

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization
dat <- sweep(dat, 2, colSums(dat), '/') * 1E6

# get relevant sample metadata
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, iss_stage = `iss:ch1`,
         patient_subgroup = `cluster:ch1`)

# add platform, cell type and disease (same for all samples)
sample_metadata$disease_stage <- 'MM'
sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

# Note; there is not sufficient information provided to link patients in table S11 to
# GSM sample identifiers; skipping.
# table_s11 <- read_csv(file.path(base_dir, 'metadata', 'broyl2010_supp_table_s11.csv'))

# load survival metadata from Kuiper et al. (2012)
survival_mdata <- read_csv('supp/clean/kuiper2012_supp_patient_survival.csv', col_types = cols())

survival_mdata <- survival_mdata %>%
  rename(geo_accession = Patient) %>%
  filter(geo_accession %in% colnames(dat))

colnames(survival_mdata) <- c('geo_accession', 'os_time', 'os_event', 'pfs_time', 'pfs_event')

# exclude samples without metadata
mask <- sample_metadata$geo_accession %in% survival_mdata$geo_accession

# ANNOT (Dec 15, 2021)
# table(mask)
# mask
# FALSE  TRUE
#    46   282

dat <- dat[, mask]
sample_metadata <- sample_metadata[mask, ]

#all(colnames(dat) == sample_metadata$geo_accession)
# [1] TRUE

# combine metadata
sample_metadata <- sample_metadata %>%
  inner_join(survival_mdata, by = 'geo_accession')

# map probes using biomart
platform <- pdata$platform_id[1]

# retrieve up-to-date probe to gene mappings
probe_mapping <- get_biomart_mapping(rownames(dat), platform, 
                                     ensembl_version = snakemake@config$biomart$ensembl_version)

# TODO: add annotations describing what is being lost here / discussing alternative approaches...
probe_mask <- rownames(dat) %in% probe_mapping$probe_id

# some probes are lost at this step of the mapping..
table(probe_mask)

# probe_mask
# FALSE  TRUE
# 12009 42666

dat <- dat[probe_mask, ]

# make sure orders match
probe_mapping <- probe_mapping[match(rownames(dat), probe_mapping$probe_id), ]

if (!all(probe_mapping$probe_id == rownames(dat))) {
  stop("Unexpected probe mapping mismatch!")
}

# get expression data and swap ensgenes for gene symbols
expr_dat <- dat %>%
  add_column(symbol = probe_mapping$symbol, .before = 1) %>%
  as_tibble()

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(probe_mapping, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
