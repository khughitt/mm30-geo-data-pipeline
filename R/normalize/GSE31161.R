#!/bin/env/Rscript
#
# Identification of multiple risk loci and regulatory mechanisms influencing
# susceptibility to multiple myeloma
#
# Went et al. (2018)
#
library(GEOquery)
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

# there are three samples in GSE31161 with all missing data, which are excluded below..
num_missing <- apply(dat, 2, function(x) {
  sum(is.na(x))
})

table(num_missing)
# num_missing
#     0 54675
#  1035     3

dat <- dat[, num_missing == 0]

# exclude outlier samples;
# these samples were found to have very median pairwise correlations (0.38 - 0.68)
# compared to all other samples (>0.9 for most)
exclude_samples <- c("GSM771497", "GSM772341", "GSM772335")

dat <- dat[, !colnames(dat) %in% exclude_samples]

sample_metadata <- sample_metadata %>%
  filter(geo_accession %in% colnames(dat))

# size factor normalization
dat <- sweep(dat, 2, colSums(dat), '/') * 1E6

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id,
         treatment = `treatment:ch1`, time_of_testing = `time of testing:ch1`) %>%
  mutate(relapsed = time_of_testing == 'relapse')

# add disease stage
sample_metadata$disease_stage <- ifelse(sample_metadata$relapsed, 'RRMM', 'MM')

# add platform and cell type (same for all samples)
sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

platform <- pdata$platform_id[1]

# retrieve up-to-date probe to gene mappings
probe_mapping <- get_biomart_mapping(rownames(dat), platform, ensembl_version = 105)

# TODO: add annotations describing what is being lost here / discussing alternative approaches...
probe_mask <- rownames(dat) %in% probe_mapping$probe_id

# some probes are lost at this step of the mapping..
table(probe_mask)
# probe_mask
# FALSE  TRUE 
# 11675 43000 

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
#annot <- list("data-prep" = read_file("annot/prepare-data/GSE31161.md"))

pkg <- pkgr$update_package(snakemake@input[[4]], mdata, "Reprocess data",
                           resources)
                           # resources, annotations = annot)

pkg %>%
  write_package(dirname(snakemake@output[[1]]))
