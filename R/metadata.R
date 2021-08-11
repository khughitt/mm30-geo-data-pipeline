#!/bin/env Rscript
#
# Build a dataset summary table using metadata from GEO
#
library(GEOquery)
library(readr)
options(stringAsFactors = FALSE)

ids <- list.files("/data/raw")

geo_metadata <- NULL

for (accession in ids) {
  dir_ <- file.path('/data/raw', accession)

  eset <- getGEO(accession, destdir = dir_)[[1]]

  mdat <- c(
    accession,
    eset@experimentData@title,
    eset@experimentData@name,
    eset@experimentData@abstract,
    eset@experimentData@other$overall_design,
    eset@experimentData@other$submission_date,
    eset@experimentData@other$last_update_date,
    eset@experimentData@other$platform_id,
    eset@experimentData@other$type,
    eset@experimentData@url,
    eset@experimentData@pubMedIds,
    eset@experimentData@other$supplementary_file
  )

  message(sprintf("%s: %s", accession, ncol(mdat))) 
  geo_metadata <- rbind(geo_metadata, mdat)
}

geo_metadata <- as.data.frame(geo_metadata)

colnames(geo_metadata) <- c('geo_id', 'title', 'name', 'abstract',
                            'overall_design', 'submission_date', 'last_update_date',
                            'platform_id', 'type', 'url', 'pubmed_ids',
                            'supplementary_file')

write_tsv(geo_metadata, '/data/metadata.tsv')

