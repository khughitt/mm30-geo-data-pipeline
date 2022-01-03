#!/bin/env/Rscript
#
# Gene expression profiling for molecular classification of multiple myeloma in newly
# diagnosed patients
#
# Broyl et al. (2010)
#
# Clinical trial: http://www.hovon.nl/studies/studies-per-ziektebeeld/mm.html?action=showstudie&studie_id=5&categorie_id=3
#
library(tidyverse)
library(iodag)
source("R/util/biomart.R")

# directory to store raw and processed data
data_dir <- dirname(snakemake@output[[1]])

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
probe_mapping <- get_biomart_mapping(rownames(dat), platform, ensembl_version = 105)

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

# create a new data package, based off the old one
pkgr <- Packager$new()

# resource list (consider converting values to lists and including "data_type" field for
# each resource? i.e. "data"/"row metadata"/"column metadata")
resources <- list(
  "data" = expr_dat,
  "row-metadata" = probe_mapping,
  "column-metadata" = sample_metadata
)

# changes to make to metadata
mdata <- list(
  data = list(
    processing = "reprocessed"
  ),
  rows = "symbol"
)

# annotations
annot <- list("data-prep" = read_file("annot/prepare-data/GSE19784.md"))

pkg <- pkgr$update_package(snakemake@input[[4]], mdata, "Reprocess data",
                           resources, annotations = annot)

pkg %>%
  write_package(dirname(snakemake@output[[1]]))
