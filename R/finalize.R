###############################################################################
#
# Generate non-redundant gene expression matrices
#
# Collapses gene symbols with multiple entries into a single row by taking the
# median expression value within each sample.
#
###############################################################################
library(tidyverse)
library(arrow)

# directory to store raw and processed data
data_dir <- dirname(snakemake@output[[1]])

# load data & metadata
dat <- read_feather(snakemake@input[[1]])

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# check to make sure no missing values are present in the normalized version of the
# dataset
if (sum(is.na(dat)) > 0) {
  stop(sprintf("Missing values encountered for dataset: %s!", snakemake@output[[1]]))
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
dat <- dat %>%
  group_by(symbol) %>%
  summarize_all(median)

write_feather(dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
