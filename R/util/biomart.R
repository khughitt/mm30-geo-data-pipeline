#
# queries biomart for a mapping from probe to gene symbols
#
# NOTE: in cases where multiple genes are associated with a probe, the approach
# arbitrarily selects the first gene as the "true" mapping...
#
get_biomart_mapping <- function(probe_ids, platform, ensembl_version=105) {
  # load biomaRt
  mart <- biomaRt::useEnsembl(biomart = "genes", 
                              dataset = "hsapiens_gene_ensembl",
                              version = ensembl_version)

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
    values = probe_ids,
    mart = mart,
    uniqueRows = TRUE
  )
  colnames(res) <- c("probe_id", "ensgene")

  # add gene symbols from annotables
  res$symbol <- annotables::grch38$symbol[match(res$ensgene, annotables::grch38$ensgene)]

  # drop entries with missing gene symbols
  num_missing <- sum(is.na(res$symbol))
  message(sprintf("Excluding %d/%d genes due to missing gene symbols.", num_missing, nrow(res)))

  res <- res[!is.na(res$symbol), ]

  res
}
