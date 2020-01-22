#!/bin/env/Rscript
#
# Single cell RNA sequencing of multiple myeloma I
#
# Ryu et al. (2019)
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE106218'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo', accession)

raw_data_dir <- file.path(base_dir, 'raw')
clean_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, clean_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# expression data is missing from getGEO() query result and must be downloaded
# separately; further, some useful clinical metadata is available as a separate file
# which we will also download using the getGEOSuppFiles() function

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

supp_file1 <- file.path(raw_data_dir, accession, 'GSE106218_GEO_processed_MM_raw_TPM_matrix.txt.gz')
supp_file2 <- file.path(raw_data_dir, accession, 'GSE106218_GEO_clinical_information_MM.txt.gz')

if (!file.exists(supp_file1)) {
  getGEOSuppFiles(accession, baseDir = raw_data_dir)
}

# load TPM counts
expr <- read.delim(gzfile(supp_file1), row.names = 1)

# load sample metadata
clinical_metadata <- as.data.frame(t(read.delim(gzfile(supp_file2), row.names = 1))) %>%
  rownames_to_column('patient_id') %>%
  select(patient_id,
         iss_stage = `ISS stage`, os_event = `Death(alive=0; death=1)`,
         os_time = `Survival time(EM;month)`,
         heavy_chain = `Heavy chain`,
         light_chain = `Light chain`) %>%
  mutate(os_event = os_event == "1")

# drop the month component of survival time
clinical_metadata$os_time <- as.numeric(str_match(clinical_metadata$os_time, '[0-9]+'))  

sample_metadata <- sample_metadata %>%
  inner_join(clinical_metadata, by = 'patient_id')

# expr
#         MM02_38 MM02_48 MM02_50
# 5S_rRNA       0       0       0
# 7SK           0       0       0
# A1BG          0       0       0

# clinical_metadata
#                         MM02 MM16 MM17
# Sex                        M    F    F
# ISS stage                 II    I    I
# Death(alive=0; death=1)    1    1    1

#table(rowSums(expr) == 0)
# 
# FALSE  TRUE 
# 35582 20278 

# remove empty rows
expr <- expr[rowSums(expr) > 0, ]

sample_metadata$disease <- "Multiple Myeloma"
sample_metadata$cell_type <- 'CD138+'

# collapse patient samples
expr_patient_ids <- str_match(colnames(expr), 'MM[0-9]+')

expr_combined <- NULL

for (patient_id in unique(expr_patient_ids)) {
  dat <- rowMeans(expr[, expr_patient_ids == patient_id])
  expr_combined <- cbind(expr_combined, dat)
}
colnames(expr_combined) <- unique(expr_patient_ids)

expr_combined <- expr_combined %>%
  as.data.frame() %>%
  rownames_to_column('symbol')

# limit sample metadata to first sample for each patient
sample_metadata <- sample_metadata %>%
  group_by(patient_id) %>%
  slice(1)

#all(colnames(expr_combined) == sample_metadata$patient_id)
# [1] TRUE

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.csv', accession)
mdat_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_combined, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, mdat_outfile))

sessionInfo()
