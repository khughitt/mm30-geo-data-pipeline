#!/bin/env/Rscript
###############################################################################
#
# GEO Dataset downloader - general version
# V. Keith Hughitt
#
# Uses GEOquery to download expression data, sample metadata, and gene annotations for a
# single dataset, and stores each as a tabular data package.
#
# For some datasets, the expression data is not stored correctly on GEO, or other issues
# are present, which necessitates a separate script to acquire the data.
#
# Scripts for each such custom data download can be found in the "custom/" folder.
#
###############################################################################
library(GEOquery)
library(tidyverse)

download_geo <- function(acc, cache_dir, dat_outfile, row_mdata_outfile, col_mdata_outfile, annot_gpl = FALSE) {
  # work-around for loading large gzip-compressed csv files (e.g. GSE31161)
  # https://github.com/r-lib/vroom/issues/361
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)

  # create cache dir
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, mode = "0755")
  }

  # download & GEO data
  eset <- getGEO(acc, destdir = cache_dir, AnnotGPL = annot_gpl)[[1]]

  expr_dat <- exprs(eset) %>%
    as.data.frame() %>%
    rownames_to_column("feature")

  pdata <- pData(eset)
  fdata <- fData(eset)

  # store results
  write_csv(expr_dat, dat_outfile)
  write_csv(fdata, row_mdata_outfile)
  write_csv(pdata, col_mdata_outfile)
}
