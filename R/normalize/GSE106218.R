#!/bin/env/Rscript
###############################################################################
#
# GSE106218
#
###############################################################################
library(GEOquery)
library(annotables)
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE) %>%
  column_to_rownames("feature")

pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

acc <- snakemake@wildcards[["acc"]]

# load "clinical information" supplemental file downloaded in previous step
supp_file <- "/data/raw/GSE106218/GSE106218_GEO_clinical_information_MM.txt.gz"

# load sample metadata
clinical_metadata <- as.data.frame(t(read.delim(gzfile(supp_file), row.names = 1))) %>%
  rownames_to_column("patient_id") %>%
  select(patient_id,
         iss_stage = `ISS stage`, os_event = `Death(alive=0; death=1)`,
         os_time = `Survival time(EM;month)`,
         heavy_chain = `Heavy chain`,
         light_chain = `Light chain`) %>%
  mutate(os_event = os_event == "1",
         disease_stage = "MM")

# drop the month component of survival time
clinical_metadata$os_time <- as.numeric(str_extract(clinical_metadata$os_time, "[0-9]+"))

sample_metadata <- pdata %>%
  select(geo_accession, platform_id,
         patient_id = `patient id:ch1`,
         sex = `gender:ch1`,
         prep_site = `prep-site:ch1`)

sample_metadata <- sample_metadata %>%
  inner_join(clinical_metadata, by = "patient_id")

# clinical_metadata
#                         MM02 MM16 MM17
# Sex                        M    F    F
# ISS stage                 II    I    I
# Death(alive=0; death=1)    1    1    1

#table(rowSums(dat) == 0)
#
# FALSE  TRUE
# 35582 20278

# remove empty rows
dat <- dat[rowSums(dat) > 0, ]

sample_metadata$platform_type <- "RNA-Seq"
sample_metadata$sample_type <- "Patient"

# collapse patient samples
expr_patient_ids <- str_extract(colnames(dat), "MM[0-9]+")

# sort(table(expr_patient_ids))
# expr_patient_ids
# MM33 MM16 MM28 MM30 MM17 MM36 MM38 MM25 MM02 MM34
#   13   24   28   31   38   47   49   52   77  129

expr_dat <- NULL

for (patient_id in unique(expr_patient_ids)) {
  combined_dat <- rowSums(dat[, expr_patient_ids == patient_id])
  expr_dat <- cbind(expr_dat, combined_dat)
}
colnames(expr_dat) <- unique(expr_patient_ids)

# collapse sample metadata to match dimensions of expression data;
# for each patient, we will keep one representative geo accession to use for the
# aggregated sample id
sample_metadata <- sample_metadata %>%
  group_by(patient_id) %>%
  arrange(geo_accession) %>%
  slice(1)

if (!all(sample_metadata$patient_id == colnames(expr_dat))) {
  stop("Column mismatch!")
}

# use geo accession ids for column names
colnames(expr_dat) <- sample_metadata$geo_accession

expr_dat <- expr_dat %>%
  as.data.frame() %>%
  rownames_to_column("symbol")

# update GRCh37 symbols -> GRCh38
grch37_mask <- expr_dat$symbol %in% grch37$symbol

gene_symbols <- expr_dat$symbol[grch37_mask]
ensgenes <- grch37$ensgene[match(gene_symbols, grch37$symbol)]

mapped_mask <- ensgenes %in% grch38$ensgene
gene_symbols[mapped_mask] <- grch38$symbol[match(ensgenes[mapped_mask], grch38$ensgene)]

expr_dat$symbol[grch37_mask] <- gene_symbols

# table(expr_dat$symbol %in% grch38$symbol)
#
# FALSE  TRUE
#  1242 34340

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
