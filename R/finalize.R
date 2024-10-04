###############################################################################
#
# Generate non-redundant gene expression matrices
#
# Collapses gene symbols with multiple entries into a single row by taking the
# median expression value within each sample.
#
###############################################################################
library(tidyverse)
library(annotables)
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

# collapse duplicated entries and exclude genes which could not be mapped to GRCh38
dat <- dat %>%
  filter(symbol %in% grch38$symbol) %>%
  group_by(symbol) %>%
  summarize(across(everything(), median))

write_feather(dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
