#!/bin/env/Rscript
#
# Molecular Signatures of Multiple Myeloma Progression through Single Cell RNA-Seq
#
# Jang et al. (2019)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(arrow)

# GEO accession
accession <- 'GSE118900'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/3.0', accession)

raw_data_dir <- file.path(base_dir, 'raw')
processed_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data
#
# eset for this dataset only include metadata; expression data is empty and will be
# loaded separately..
#
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

supp_file <- file.path(raw_data_dir, accession, 'GSE118900_MM.scrna-seq.tpm.pass.txt.gz')

if (!file.exists(supp_file)) {
  getGEOSuppFiles(accession, baseDir = raw_data_dir, filter_regex = 'tpm.pass')
}

expr_dat <- read.delim(gzfile(supp_file), row.names = 1)

# replace rownames with gene symbol

#head(rownames(expr_dat))
# [1] "AADACL3|chr1|12776118" "AADACL4|chr1|12704566" "ABCA4|chr1|94458394"   "ABCB10|chr1|229652329"
# [5] "ABCD3|chr1|94883933"   "ABL2|chr1|179068462"

gene_symbols <- str_split(rownames(expr_dat), '\\|', simplify = TRUE)[, 1]

# load GRCh38 gene symbol mapping
gene_mapping <- read_tsv('../../annot/GRCh38_alt_symbol_mapping.tsv', col_types = cols())

# mask indicating which genes are to be updated
mask <- !gene_symbols %in% grch38$symbol & gene_symbols %in% gene_mapping$alt_symbol

#table(mask)
# mask
# FALSE  TRUE 
# 21156  2242 

gene_symbols[mask] <- gene_mapping$symbol[match(gene_symbols[mask],
                                                gene_mapping$alt_symbol)]

# remove small number of duplicated gene symbol entries
# mask <- !duplicated(gene_symbols)

# expr_dat <- expr_dat[mask, ]
# rownames(expr_dat) <- gene_symbols[mask]

# size factor normalization
expr_dat <- sweep(expr_dat, 2, colSums(expr_dat), '/') * 1E6

expr_dat <- expr_dat %>%
  as.data.frame() %>%
  add_column(symbol = gene_symbols, .before = 1)

# exclude any zero variance genes present
mask <- apply(expr_dat[, -1], 1, var) > 0
expr_dat <- expr_dat[mask, ]

# table(mask)
# mask
# FALSE  TRUE 
#  3789 19609 

# columns to include (GSE118900)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         sample_name = description,
         disease_stage = `tumor stage:ch1`,
         patient = `patient:ch1`,
         mm_stage = `tumor stage:ch1`) %>%
  mutate(disease_stage = recode(disease_stage, `IgM-MGUS` = 'MGUS', `NDMM` = 'MM'))

# add cell type 
sample_metadata$cell_type <- 'CD138+'

# normalize sample ids, ex: "IgM-MGUS1_C37" -> "IgM.MGUS1_C37"
sample_metadata$sample_name <- gsub('-', '.', sample_metadata$sample_name)

# match sample order to metadata
expr_dat <- expr_dat[, c('symbol', sample_metadata$sample_name)]

if (!all(colnames(expr_dat)[-1] == sample_metadata$sample_name)) {
  stop("Sample ID mismatch!")
}

# for consistency, use GEO sample accessions as ids
colnames(expr_dat) <- c('symbol', sample_metadata$geo_accession)

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
