#!/bin/env/Rscript
#
# Gene expression profiling of B lymphocytes and plasma cells from Waldenstrom's
# macroglobulinemia
#
# Gutiérrez et al. (2007)
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("../util/eset.R")

# GEO accession
accession <- 'GSE6691'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/3.1', accession)

raw_data_dir <- file.path(base_dir, 'raw')
processed_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

#head(colnames(eset))
# [1] "5.82592"  "6.799907" "5.388776" "6.630125" "7.093243" "6.426071"

# NOTE: GEOquery currently incorrectly parses the data for this experiment, skipping
# >500 lines and usign the wrong colnames; reported the issue upstream. For now, just
# fixing colnames..
colnames(eset) <- pData(eset)$geo_accession

# columns to include (GSE6699)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

# add cell type and disease
sample_metadata$cell_type <- 'CD138+'

sample_metadata$platform_type <- 'Microarray'

#table(pData(eset)$characteristics_ch1)
#
#                                                                                         BL from CLL patient, isolated using CD19-PE and CD5-APC
#                                                                                                                                               2
#                                                                                        BL from CLL patient, isolated using CD19-PE and CD5-APC.
#                                                                                                                                               9
#                                                                                             BL from healthy donors, isolated using CD19-PE-Cy7.
#                                                                                                                                               8
# BL from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE, CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7.
#                                                                                                                                               9
#    from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE, CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7.
#                                                                                                                                               1
#                                                                                                PC from healthy donors, isolated using CD38-APC.
#                                                                                                                                               5
#                                                                                                     PC from MM patient, isolated using CD38-APC
#                                                                                                                                               2
#                                                                                                    PC from MM patient, isolated using CD38-APC.
#                                                                                                                                               9
#                                                                                                     PC from MM patient,isolated using CD38-APC.
#                                                                                                                                               1
# PC from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE, CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7.
#                                                                                                                                               9
#  PC from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE,CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7.
#                                                                                                                                               1
patient_desc <- pData(eset)$characteristics_ch1

disease <- rep('MM', length(patient_desc))
disease[grepl('healthy', patient_desc)] <- "Healthy"
disease[grepl('WM', patient_desc)] <- "Waldenström's Macroglobulinemia"
disease[grepl('CLL', patient_desc)] <- "Chronic Lymphocytic Leukemia"

table(disease)
# disease
#    Chronic Lymphocytic Leukemia                         Healthy                              MM Waldenström's Macroglobulinemia 
#                              11                              13                              12                              20 
sample_metadata$mm_stage <- disease
sample_metadata$disease_stage <- disease

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# keep only healthy / myeloma samples
mask <- sample_metadata$disease_stage %in% c('Healthy', 'MM')

#table(mask)
# mask
# FALSE  TRUE
#    31    25

eset <- eset[, mask]
sample_metadata <- sample_metadata[mask, ]

# extract gene expression data
expr_dat <- process_eset(eset)

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  group_by(symbol) %>%
  summarize_all(median)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
expr_nr_outfile <- sprintf('%s_gene_expr_nr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_feather(expr_dat_nr, file.path(processed_data_dir, expr_nr_outfile))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))
