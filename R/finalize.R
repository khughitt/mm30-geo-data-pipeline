###############################################################################
#
# Generate non-redundant gene expression matrices
#
# Collapses gene symbols with multiple entries into a single row by taking the
# median expression value within each sample.
#
###############################################################################
library(tidyverse)
library(iodat)
library(jsonlite)

# directory to store raw and processed data
data_dir <- dirname(snakemake@output[[1]])

# load data & metadata
dat_orig <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

row_mdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
col_mdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# load previous data package (used to retrieve styles metadata)
prev_pkg <- suppressMessages(read_package(snakemake@input[[4]]))

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
dat <- dat_orig %>%
  group_by(symbol) %>%
  summarize_all(median)

# perform sample pca projection
pca <- prcomp(t(as.matrix(dat[, -1])), scale = snakemake@params[["scale_pca"]])

# convert pca-projected data to a dataframe
pca_dat <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column('geo_accession')

# store summary table in data package
pca_summary <- summary(pca)$importance %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('PC')

# add plot styles to pca dataframe
style_fields <- c(color = prev_pkg$io$styles$columns$color)

if ("shape" %in% names(prev_pkg$io$styles$columns)) {
  style_fields <- c(style_fields, shape = prev_pkg$io$styles$columns$shape)
}

style_dat <- col_mdata %>%
  select(geo_accession, style_fields)

pca_dat <- pca_dat %>%
  inner_join(style_dat, by = 'geo_accession')

# load pca view
# TODO: extend with shape, if present?..
#pca_view <- jsonlite::read_json("views/sample-pca.json")

# create a new data package, based off the old one
pkg_dir <- dirname(snakemake@output[[1]])
#setwd(pkg_dir)

pkgr <- Packager$new()

# node-level metadata
node_mdata <- list(
  processing = "non-redundant"
)

resources <- list(
  "data" = dat,
  "row-metadata" = row_mdata,
  "column-metadata" = col_mdata,
  "pca-data" = pca_dat,
  "pca-summary" = pca_summary
)

pkgr$update_package(snakemake@input[[4]], resources,
                    node_metadata = node_mdata,
                    views = list("pca" = "views/sample-pca.json"),
                    pkg_dir = pkg_dir)
