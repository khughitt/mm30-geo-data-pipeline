#
# Helper function to map from microarray probe ID's to updated gene symbol assignments.
# Probes mapping to multiple gene ids are dropped from the dataset.
#
# V. Keith Hughitt
#

#
# process geo eset using biomart to map probe ids
#
process_eset <- function(eset, ensembl_version=104) {
  print("Calling get_biomart_mapping()...")

  # generate probe -> gene mapping
  probe_mapping <- get_biomart_mapping(eset, ensembl_version)

  # extend mapping manually
  missing_mask <- !rownames(eset) %in% probe_mapping$probe_id
  alt_mapping <- get_manual_mapping(eset[missing_mask, ])

  num_biomart <- length(unique(probe_mapping$probe_id))
  num_manual <- length(unique(alt_mapping$probe_id))

  probe_mapping <- rbind(probe_mapping, alt_mapping)

  print(sprintf("Mapped %d / %d probes to gene symbols (biomart: %d, manual: %d)",
                num_biomart + num_manual, nrow(eset),
                num_biomart, num_manual))

  # drop any duplicate entries in mapping; should be few
  probe_mapping <- probe_mapping[!duplicated(probe_mapping), ]

  # in some cases, probes may be mapped to a very large number of genes (e.g. 30-50..)
  # should a limit be set on the maximum number to allow?
  # exon arrays are an extreme example of this where ~95% of probes may map to
  # multiple genes..

  # extract gene expression data and convert to a dataframe
  dat <- exprs(eset) %>%
    as.data.frame()

  # drop any probes that could not be mapped
  dat <- dat[rownames(dat) %in% probe_mapping$probe_id, ]

  # split each multimapped probe into multiple entries with the same values; one row
  # per gene mapped to.
  res <- NULL
  res_ids <- c()

  # not the fastest approach, but good enough for now..
  for (i in seq_len(nrow(dat))) {
    probe_id <- rownames(dat)[i]

    # get ensembl gene ids mapped to probe
    ind <- which(probe_mapping$probe_id == probe_id)
    ensgenes <- probe_mapping$ensgene[ind]

    # in some cases, multiple ensembl id's map to the same symbol and can be safely
    # excluded
    gene_symbols <- annotables::grch38$symbol[match(ensgenes, annotables::grch38$ensgene)]
    gene_symbols <- unique(gene_symbols)

    # map to gene symbols and store
    res_ids <- c(res_ids, gene_symbols)

    # duplicate expression values, one for each gene and append to result dataframe
    res <- rbind(res, dat[rep(i, length(gene_symbols)), ])

    # sanity check..
    if (length(res_ids) != nrow(res)) {
      stop("Dimension mismatch!")
    }
  }

  dat <- res

  # drop probes that could not be mapped
  mask <- !is.na(res_ids)

  print(sprintf("Dropping %d / %d probes which could not be mapped to genes.",
                sum(!mask), length(mask)))

  dat <- dat[mask, ]
  res_ids <- res_ids[mask]

  # load expression data and add gene symbol column
  dat %>%
    add_column(symbol = res_ids, .before = 1)
}


#
# queries biomart for a mapping from probe to gene symbols
#
get_biomart_mapping <- function(eset, ensembl_version=104) {
  # load biomaRt
  mart <- biomaRt::useEnsembl(biomart = 'genes', 
                              dataset = 'hsapiens_gene_ensembl',
                              version = ensembl_version)

  # determine biomart platform attribute
  platform <- pData(eset)[1, ]$platform_id

  platform_mapping <- data.frame(
    geo = c("GPL96", "GPL97", "GPL570", "GPL10558", "GPL5175", "GPL6244"),
    biomart = c("affy_hg_u133a", "affy_hg_u133b", "affy_hg_u133_plus_2",
                "illumina_humanht_12_v4", "affy_huex_1_0_st_v2", "affy_hugene_1_0_st_v1")
  )
  platform_attr <- platform_mapping$biomart[match(platform, platform_mapping$geo)]

  # query biomart for mapping;
  #
  # note: in testing, ensembl gene ids were found to be more reliable to map to GRCh38
  # symbols than using the symbols returned by biomaRt
  res <- biomaRt::getBM(
    attributes = c(
      platform_attr,
      "ensembl_gene_id"),
    filters = platform_attr,
    values = rownames(eset),
    mart = mart,
    uniqueRows = TRUE
  )
  colnames(res) <- c('probe_id', 'ensgene')

  res
}

#
# generate probe -> gene symbol mapping for an eset using manual approach
#
get_manual_mapping <- function(eset) {
  # mapping from old/alt gene symbols to versions used in GRCh38
  grch38_mapping <- read_tsv('identifiers/GRCh38_alt_symbol_mapping.tsv', col_types = cols())

  # determine name of gene symbol field in fData
  fdata_fields <- c('Gene Symbol', 'Gene symbol', 'Symbol')
  fdata_field <- fdata_fields[which(fdata_fields %in% colnames(fData(eset)))]

  if (is.null(fdata_field)) {
    stop("Cannot determine which fData() field to use!")
  }

  # split multi-mapped gene symbols
  alt_mapping <- exprs(eset) %>%
    as.data.frame() %>%
    rownames_to_column('probe_id') %>%
    add_column(symbol = fData(eset)[, fdata_field], .before = 1) %>%
    separate_rows(symbol, sep = " ?//+ ?") %>%
    select(probe_id, symbol)

  # add ensgene placeholder
  alt_mapping$ensgene <- NA

  alt_mapping <- alt_mapping %>%
    select(probe_id, ensgene, symbol)

  # get ensgenes and GRCh38 symbols, where possible
  mask <- alt_mapping$symbol %in% grch38_mapping$alt_symbol
  ind <- match(alt_mapping$symbol[mask], grch38_mapping$alt_symbol)

  alt_mapping$symbol[mask] <- grch38_mapping$symbol[ind]
  alt_mapping$ensgene[mask] <- grch38_mapping$ensgene[ind]

  # in cases where symbol is present, but not ensgene, attempt to add ensgene
  mask2 <- is.na(alt_mapping$ensgene) & alt_mapping$symbol %in% annotables::grch38$symbol

  ind <- match(alt_mapping$symbol[mask2], annotables::grch38$symbol)
  alt_mapping$ensgene[mask2] <- annotables::grch38$ensgene[ind]

  # drop any probes that could not be mapped manually
  alt_mapping <- alt_mapping[!is.na(alt_mapping$ensgene), ]

  # return probe/ensgene mapping
  alt_mapping %>%
    select(probe_id, ensgene) %>%
    as.data.frame()
}
