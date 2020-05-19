#!/bin/env/Rscript
#
# Gene expression profiling for molecular classification of multiple myeloma in newly
# diagnosed patients
#
# Broyl et al. (2010)
#
# Clinical trial: http://www.hovon.nl/studies/studies-per-ziektebeeld/mm.html?action=showstudie&studie_id=5&categorie_id=3
#
library(GEOquery)
library(tidyverse)
library(arrow)

# GEO accession
accession <- 'GSE19784'

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
# result is a list with a single entry containing an ExpressionSet instance
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# exclude control sequences present in some datasets (GSE19784)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# size factor scaling has already been performed
#range(colSums(exprs(eset)))
# [1] 212082.1 239655.9

# list the columns that are neither all different or all unique
covariates <- c()

for (cname in colnames(pData(eset))) {
  num_vals <- length(table(pData(eset)[, cname]))

  if (num_vals != 1 && num_vals != ncol(eset)) {
    covariates <- c(covariates, cname)

    message(cname)
    print(table(pData(eset)[, cname]))
  }
}

#
# fields of interest:
#
# cluster:ch1
# CD-1       CD-2        CTA         HY         MF         MS    Myeloid       NFκB       None
# 13         34         22         80         32         33         40         39          9
# PR       PRL3 SOCS3/PRL3
# 15          2          9
#iss:ch1
#  I  II III  nd
#122  88  83  27
#

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, iss_stage = `iss:ch1`,
         patient_subgroup = `cluster:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# Note; there is not sufficient information provided to link patients in table S11 to
# GSM sample identifiers; skipping.
#table_s11 <- read_csv(file.path(base_dir, 'metadata', 'broyl2010_supp_table_s11.csv'))

# load survival metadata from Kuiper et al. (2012)
survival_mdata <- read_csv('/data/human/kuiper2012/kuiper2012_supp_patient_survival.csv', col_types = cols())

survival_mdata <- survival_mdata %>%
  rename(geo_accession = Patient) %>%
  filter(geo_accession %in% colnames(eset))

colnames(survival_mdata) <- c('geo_accession', 'os_time', 'os_event', 'pfs_time', 'pfs_event')

# exclude samples without metadata
mask <- sample_metadata$geo_accession %in% survival_mdata$geo_accession

#all(colnames(eset) == sample_metadata$geo_accession)
# [1] TRUE

eset <- eset[, mask]
sample_metadata <- sample_metadata[mask, ]

# combine metadata
sample_metadata <- sample_metadata %>%
  inner_join(survival_mdata, by = 'geo_accession')

# load expression data and add gene symbol column
expr_dat <- as.data.frame(exprs(eset))

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  add_column(symbol = fData(eset)$`Gene symbol`, .before = 1) %>%
  filter(symbol != '')

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

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
