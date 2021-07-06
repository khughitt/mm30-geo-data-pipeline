#!/bin/env/Rscript
#
# Key transcription factors altered in multiple myeloma patients revealed by logic
# programming approach combining gene expression pro ling and regulatory networks
#
# Miannay et al. (2016)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(arrow)

# GEO accession
accession <- 'GSE83503'

# directory to store raw and processed data
raw_data_dir <- file.path('/data/raw', accession)
processed_data_dir <- sub('raw', 'clean', raw_data_dir)

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data;
# result is a list with a single entry containing an ExpressionSet instance
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# exclude control spots (4130 / 22011 probes)
mask <- fData(eset)$SPOT_ID != 'control'
eset <- eset[mask, ]

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var) > 0, ]

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# get relevant sample metadata
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, patient_died = `death:ch1`, pfs_event = `relapse:ch1`) %>%
  mutate(pfs_event = ifelse(pfs_event == 0, 0, 1))

# add cell type and disease stage
sample_metadata$disease_stage <- ifelse(sample_metadata$pfs_event == 1, 'RRMM', 'MM')
sample_metadata$cell_type <- 'CD138+'

sample_metadata$platform_type <- 'Microarray'

# Note: GSE83503 was performed on an Affymetrix Human Exon 1.0 ST Array, with multiple
# probes for each exon. The result of this is that >95% of the probes map to multiple
# genes, and thus simply discarding multi-mapped probes will not be helpful.
# get gene symbols associated with each probe; gene symbols are stored at every

# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fData(eset)$gene_assignment, '//', simplify = TRUE)
gene_symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
gene_symbols <- apply(gene_symbols, 1, function(x) {
  str_trim(x)[x != '']
})
gene_symbols <- unlist(lapply(gene_symbols, paste, collapse = ' // '))

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
  add_column(symbol = gene_symbols, .before = 1) %>%
  filter(symbol != '')

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# split multi-mapped symbols
expr_dat_nr <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

# load GRCh38 gene symbol mapping
gene_mapping <- read_tsv('../annot/GRCh38_alt_symbol_mapping.tsv', col_types = cols())

# mask indicating which genes are to be updated
mask <- !expr_dat$symbol %in% grch38$symbol & expr_dat$symbol %in% gene_mapping$alt_symbol

# table(mask)
# mask
# FALSE  TRUE
# 17558    76

expr_dat$symbol[mask] <- gene_mapping$symbol[match(expr_dat$symbol[mask],
                                                   gene_mapping$alt_symbol)]

# table(duplicated(expr_dat$symbol))
# FALSE  TRUE
# 17602    32

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?") %>%
  group_by(symbol) %>%
  summarize_all(median)

# store cleaned expression data and metadata
expr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr.feather', accession))
expr_outfile_nr <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

print("Final dimensions:")
print(paste0("- Num rows: ", nrow(expr_dat_nr)))
print(paste0("- Num cols: ", ncol(expr_dat_nr)))

write_feather(expr_dat, expr_outfile)
write_feather(expr_dat_nr, expr_outfile_nr)
write_tsv(sample_metadata, mdat_outfile)
