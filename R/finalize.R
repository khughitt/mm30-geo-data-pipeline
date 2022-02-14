###############################################################################
#
# Generate non-redundant gene expression matrices
#
# Collapses gene symbols with multiple entries into a single row by taking the
# median expression value within each sample.
#
###############################################################################
library(tidyverse)
library(eco)
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
pca <- prcomp(t(as.matrix(dat[, -1])), scale = snakemake@config[["pca"]][["scale"]])

# convert pca-projected data to a dataframe
pca_dat <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column('geo_accession')

# store summary table in data package
pca_summary <- summary(pca)$importance %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('PC')

# limit to specified number of pc's
num_pcs <- snakemake@config[["pca"]][["num_pcs"]]

pca_dat <- pca_dat[, 1:(num_pcs + 1)]
pca_summary <- head(pca_summary, num_pcs)

# add plot styles to pca dataframe
styles <- prev_pkg$eco$metadata$styles$columns

style_fields <- as.character(styles)

if (length(style_fields) > 0) {
  style_dat <- col_mdata %>%
    select(geo_accession, all_of(style_fields))

  pca_dat <- pca_dat %>%
    inner_join(style_dat, by = 'geo_accession')
}

# load pca vegalite view
pca_view <- jsonlite::read_json("views/sample-pca.json")

# add color/shape, if specified
if ("color" %in% names(styles)) {
  pca_view$encoding$color <- list(
    "field" = styles$color,
    "type" = "nominal"
  )
}

if ("shape" %in% names(styles)) {
  pca_view$encoding$shape <- list(
    "field" = styles$shape,
    "type" = "nominal"
  )
}

# create a new data package, based off the old one
pkg_dir <- dirname(snakemake@output[[1]])

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
                    views = list("pca" = pca_view),
                    pkg_dir = pkg_dir)
