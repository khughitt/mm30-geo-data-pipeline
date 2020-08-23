#!/bin/env/Rscript
#
#	Gene Expression profiling of Multiple Myeloma
#
# Chauhan et al. (2012)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(arrow)

options(stringsAsFactors = FALSE)

# GEO accession
accession <- 'GSE39754'

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

# download GEO data;
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         diagnosis = `diagnosis:ch1`, treatment_response = `treatment_response:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$cell_type <- 'CD138+'

sample_metadata$disease_stage <- 'MM'
sample_metadata$disease_stage[grepl('Healthy', sample_metadata$diagnosis)] <- 'Healthy'

# Note: GSE39754 was performed on an Affymetrix Human Exon 1.0 ST Array, with multiple
# probes for each exon. The result of this is that >95% of the probes map to multiple
# genes, and thus simply discarding multi-mapped probes will not be helpful.

# get gene symbols associated with each probe; gene symbols are stored at every
# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fData(eset)$gene_assign, ' ///? ', simplify = TRUE)
symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
symbols <- apply(symbols, 1, function(x) {
  str_trim(x)[x != '']
})
symbols <- unlist(lapply(symbols, paste, collapse = ' // '))

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
  add_column(symbol = symbols, .before = 1) %>%
  filter(symbol != '')

if(!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?") %>%
  group_by(symbol) %>%
  summarize_all(median)

# for genes not already using symbols in GRCh38, attempt to map the symbols 
missing_symbols <- !expr_dat_nr$symbol %in% grch38$symbol

#table(missing_symbols)
# missing_symbols
# FALSE  TRUE 
# 17266  2666 

# load GRCh38 gene symbol mapping
gene_mapping <- read_tsv('../../annot/GRCh38_alt_symbol_mapping.tsv', col_types = cols())

# mask indicating which genes are to be updated
mask <- !expr_dat_nr$symbol %in% grch38$symbol & expr_dat_nr$symbol %in% gene_mapping$alt_symbol

# table(mask)
# mask
# FALSE  TRUE 
# 18061  1871 

expr_dat_nr$symbol[mask] <- gene_mapping$symbol[match(expr_dat_nr$symbol[mask], 
                                                      gene_mapping$alt_symbol)]

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
expr_outfile_nr <- sprintf('%s_gene_expr_nr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_feather(expr_dat_nr, file.path(processed_data_dir, expr_outfile_nr))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))
