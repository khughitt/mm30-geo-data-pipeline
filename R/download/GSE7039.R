#!/bin/env/Rscript
###############################################################################
#
# Download GSE7039
#
# Notes: 
#
# 1. Survival data provided by corresponding author (Stéphane Minvielle) via email on
#    June 6, 2019.
# 2. In order to apply size factor normalization at the individual array level,
#    normalization is applied at this stage in the pipline for GSE7039 instead of at the
#    usual "normalize" pipeline step.
#
###############################################################################
library(GEOquery)
library(tidyverse)

acc <- snakemake@wildcards[["acc"]]
cache_dir <- file.path("/data/raw", acc)

# create cache dir
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, mode = "0755")
}

# load survival metadata provided by author
survival_dat <- read_csv('supp/clean/GSE7039_MM_Survival_time.csv', col_types = cols())

# download GEO data;
# for GSE7039, data includes two separate esets with a small number of overlapping
# probes; Also note that the data corresponds to two platforms, only one of which
# has a corresponding .annot.gz file.

# The GEO project includes two replicates ("_A" and "_B") for each patient.
esets <- getGEO(acc, destdir = cache_dir, AnnotGPL = TRUE)

# combine samples from separate ExpressionSets
# survival units: days
pdata <- pData(esets[[1]]) %>%
  mutate(patient_id = sub('_A', '', title)) %>%
  select(geo_accession, platform_id, patient_id) %>%
  inner_join(survival_dat, by = 'patient_id') %>%
  add_column(geo_accession2 = pData(esets[[2]])$geo_accession, .after = 1) %>%
  mutate(patient_died = ifelse(deceased == 'yes', 1, 0)) %>%
  rename(os_time = follow_up_days) %>%
  select(-deceased)

if (!all(pdata$patient_id == sub('_B', '', pData(esets[[2]])$title))) {
  stop("sample metadata mismatch!")
}

# add cell type and disease stage (same for all samples)
pdata$disease_stage <- 'MM'
pdata$platform_type <- 'Microarray'

# combine expression data from two esets
e1 <- exprs(esets[[1]])
e2 <- exprs(esets[[2]])

# remove AFFX- probes seprately and then combine
mask1 <- !startsWith(rownames(e1), 'AFFX-')
mask2 <- !startsWith(rownames(e2), 'AFFX-')

e1 <- e1[mask1, ]
e2 <- e2[mask2, ]

# perform size factor normalization separately on each microarray
e1 <- sweep(e1, 2, colSums(e1), '/') * 1E6
e2 <- sweep(e2, 2, colSums(e2), '/') * 1E6

# there are no shared probe ids across the two arrays
# length(intersect(rownames(e1), rownames(e2)))
# [1] 0

expr_dat <- rbind(e1, e2)

# exclude outlier samples;
# these samples were found to have very median pairwise correlations (0.03 - 0.17)
# compared to all other samples (average > 0.6).
exclude_samples <- c("GSM162390", "GSM162298", "GSM162472")

expr_dat <- expr_dat[, !colnames(expr_dat) %in% exclude_samples]

pdata <- pdata %>%
  filter(!geo_accession %in% exclude_samples)

# generate combined gene annotation table
fdata1 <- fData(esets[[1]])[mask1, ] %>%
  select(ID, `Gene symbol`)
fdata2 <- fData(esets[[2]])[mask2, ] %>%
  select(ID, `Gene symbol` = `Gene Symbol`)

fdata <- rbind(fdata1, fdata2)

write_csv(as.data.frame(expr_dat), snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
