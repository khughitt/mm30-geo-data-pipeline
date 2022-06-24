#!/bin/env Rscript
#
# Build a dataset summary table using metadata from GEO
#
library(GEOquery)
library(readr)
options(stringAsFactors = FALSE)

ids <- list.files("/data/raw/geo")
ids <- ids[startsWith(ids, "GSE")]

geo_metadata <- NULL

for (accession in ids) {
  dir_ <- file.path("/data/raw/geo", accession)

  eset <- getGEO(accession, destdir = dir_)[[1]]

  # supplementary file field is sometimes NULL
  supp_file <- eset@experimentData@other$supplementary_file

  if (is.null(supp_file)) {
    supp_file <- ""
  }

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
    supp_file
  )

  geo_metadata <- rbind(geo_metadata, mdat)
}

geo_metadata <- as.data.frame(geo_metadata)

colnames(geo_metadata) <- c('geo_id', 'title', 'name', 'abstract',
                            'overall_design', 'submission_date', 'last_update_date',
                            'platform_ids', 'type', 'urls', 'pubmed_ids',
                            'supplementary_files')

# replace newlines to avoid issues in rendered tsv file
geo_metadata$abstract <- gsub("\n", " ", geo_metadata$abstract)
geo_metadata$overall_design <- gsub("\n", " ", geo_metadata$overall_design)
geo_metadata$platform_ids <- gsub("\n", "; ", geo_metadata$platform_ids)
geo_metadata$pubmed_ids <- gsub("\n", "; ", geo_metadata$pubmed_ids)
geo_metadata$urls <- gsub("\n", "; ", geo_metadata$urls)
geo_metadata$supplementary_files <- gsub("\n", "; ", geo_metadata$supplementary_files)

write_tsv(geo_metadata, snakemake@output[[1]])
