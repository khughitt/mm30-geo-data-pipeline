###############################################################################
#
# Generate non-redundant gene expression matrices
#
# Collapses gene symbols with multiple entries into a single row by taking the
# median expression value within each sample.
#
# - [ ] TODO: adjust row metadata to reflect changes from aggregation..
#
###############################################################################
library(tidyverse)

# directory to store raw and processed data
data_dir <- dirname(snakemake@output[[1]])

# load data & metadata
dat_orig <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

row_mdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
col_mdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
dat <- dat_orig %>%
  group_by(symbol) %>%
  summarize_all(median)

write_csv(dat, snakemake@output[[1]])
write_csv(row_mdata, snakemake@output[[2]])
write_csv(col_mdata, snakemake@output[[3]])
