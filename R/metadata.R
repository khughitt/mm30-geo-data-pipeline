#!/bin/env Rscript
#
# Build a dataset summary table using metadata from GEO
#
library(GEOquery)
library(readr)
options(stringAsFactors = FALSE)

ids <- sub(".R", "", list.files("."))
ids <- ids[startsWith(ids, "GSE")]

geo_metadata <- NULL

for (gse in ids) {
  dir_ <- file.path('/data/human/geo', gse, 'raw')

  eset <- getGEO(gse, destdir = dir_)[[1]]

  mdat <- c(
    gse,
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

  geo_metadata <- rbind(geo_metadata, mdat)
}

geo_metadata <- as.data.frame(geo_metadata)

colnames(geo_metadata) <- c('geo_id', 'title', 'name', 'abstract',
                            'overall_design', 'submission_date', 'last_update_date',
                            'platform_id', 'type', 'url', 'pubmed_ids',
                            'supplementary_file')

write_tsv(geo_metadata, '/data/human/geo/metadata.tsv')

