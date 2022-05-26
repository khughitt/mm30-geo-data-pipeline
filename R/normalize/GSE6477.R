#!/bin/env/Rscript
###############################################################################
#
# GSE6477
#
###############################################################################
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization (ignore gene symbol column)
dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), '/') * 1E6

# add gene symbol column
expr_dat <- dat %>%
  select(-feature) %>%
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
  as_tibble()

mask <- expr_dat$symbol != ""
expr_dat <- expr_dat[mask, ]
fdata <- fdata[mask, ]

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title, 
         ploidy = characteristics_ch1,
         ch13_status = characteristics_ch1.1) %>%
  mutate(mm_stage = ifelse(grepl('Normal', title), 'Normal',
                    ifelse(grepl('New', title), 'New',
                    ifelse(grepl('MGUS', title), 'MGUS',
                    ifelse(grepl('Relapsed', title), 'Relapsed',
                    ifelse(grepl('Smoldering', title), 'Smoldering', 'Normal')))))) %>%
  select(-title)

# add disease stage
sample_metadata <- sample_metadata %>%
  mutate(disease_stage = recode(mm_stage, 
                                Normal = 'Healthy', New = 'MM',
                                Relapsed = 'RRMM', Smoldering = 'SMM'))


# add platform
sample_metadata$platform_type <- 'Microarray'

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# exclude outlier samples with low median pairwise correlations (0.63-0.64 vs. ~0.9 for
# most other samples)
exclude_samples <- c("GSM149035", "GSM149037")

expr_dat <- expr_dat[, !colnames(expr_dat) %in% exclude_samples]

sample_metadata <- sample_metadata %>%
  filter(!geo_accession %in% exclude_samples)

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])