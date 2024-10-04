#!/bin/env/Rscript
###############################################################################
#
# GSE193531
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
expr_dat <- read_feather(snakemake@input[[1]]) %>%
  column_to_rownames("feature")

pdata <- read_feather(snakemake@input[[3]])

# load cell-level metadata from supplemental file downloaded in previous step
cell_mdat <- read_csv("/data/raw/GSE193531/GSE193531_cell-level-metadata.csv.gz",
                      show_col_types = FALSE)

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, title,  platform_id,
         disease_stage = `characteristics_ch1`) %>%
  mutate(disease_stage = recode(disease_stage,
                                `diagnosis: MGUS` = "MGUS",
                                `diagnosis: MM` = "MM",
                                `diagnosis: NBM (healthy control)` = "Healthy",
                                `diagnosis: SMM` = "SMM")) %>%
  mutate(patient_id = make.names(title))

# add cell and platform type
sample_metadata$platform_type <- "RNA-Seq"
sample_metadata$sample_type <- "Patient"

# convert sample ids to valid variable names
# sample_metadata$title <- gsub("-", ".", sample_metadata$title)

# combine cell measurements for each patient
patient_ids <- unique(sample_metadata$title)

expr_list <- list(
  "symbol" = rownames(expr_dat)
)

for (patient_id in patient_ids) {
  # patient_samples <- sample_metadata %>%
  #   filter(title == patient_id) %>%
  #   pull(geo_accession)
  cell_ids <- cell_mdat %>% 
    filter(sample_ID == patient_id) %>%
    pull(index) %>%
    make.names()

  mask <- colnames(expr_dat) %in% cell_ids

  expr_list[[patient_id]] <- rowSums(expr_dat[, mask])
}
expr_dat <- data.frame(expr_list)

# size factor normalization (ignore gene symbol column)
expr_dat[, -1] <- sweep(expr_dat[, -1], 2, colSums(expr_dat[, -1]), "/") * 1E6

# exclude any zero variance genes present
mask <- apply(expr_dat[, -1], 1, var) > 0
expr_dat <- expr_dat[mask, ]

# normalize expr data / metadata order
expr_dat <- expr_dat[, c("symbol", sample_metadata$patient_id)]

if (!all(sample_metadata$patient_id == colnames(expr_dat[, -1]))) {
  stop("patient id mismatch!")
}

sample_metadata <- sample_metadata %>%
  select(patient_id, everything())

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# create gene-level metadata
ind <- match(expr_dat$symbol, grch38$symbol)

fdata <- data.frame(list(
  symbol = expr_dat$symbol,
  ensgene = grch38$ensgene[ind],
  entrez = grch38$entrez[ind],
  chr = grch38$chr[ind],
  start = grch38$start[ind],
  end = grch38$end[ind],
  strand = grch38$strand[ind],
  biotype = grch38$biotype[ind],
  description = grch38$description[ind]
))

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
