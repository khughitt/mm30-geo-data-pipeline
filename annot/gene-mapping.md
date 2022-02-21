Gene Mapping Notes
==================

_KH (Dec21)_

Summary
-------

One step in the data processing flow where many different approaches are possible, in the normalization of gene identifiers across the different datasets.

In order to make comparisons across the datasets, it is necessary to map the datasets to a common set of identifiers.

For this analysis, the GRCh38 ensembl gene ids / gene symbols are used.

Since the gene annotations retrieved from GEO (i.e. those provided by `fData(eset)`) tend to be out-of-date, biomaRt is used to retrieve more recent annotations associated with the microarray probes.

Out of the 21 microarray datasets included in this analysis, 17 were associated with platforms for which biomaRt could be queried for updated identifiers.

For the remaining four datasets, the annotations returned from `fData()` were used directly:

1. GSE117847 (Affymetrix Human Clariom D)
2. GSE39754 (Affymetrix Human Exon 1.0 ST Array)
3. GSE7039 (UMGC-IRCNA 9k A/B)
4. GSE83503 (Affymetrix Human Exon 1.0 ST Array)

Annotables vs. biomaRt
----------------------

Rather than the gene symbols returned from biomaRt (`hgnc_symbol`), [annotables](https://github.com/stephenturner/annotables) is used to map from the biomaRt ensgenes to gene symbols.

This motivation for this is to take advantage of the fact that biomaRt results tend to have many more non-missing ensembl gene identifiers than gene symbols.

For example, for the set of all gene symbols associated with `affy_hg_u133a` probes, using the gene symbols returned from biomaRt directly would result in a total of 39,536 gene symbols being assigned to the probes.

If instead, the biomaRt provided ensgene id's are used to looking gene symbols in the annotables GRCh38 annotations, a total of 57,471 probes can be mapped:

```r
library(annotables)
library(biomaRt)

platform_attr <- "affy_hg_u133a"

# load biomaRt
mart <- biomaRt::useEnsembl(biomart = 'genes',
                            dataset = 'hsapiens_gene_ensembl',
                            version = ensembl_version)

res <- biomaRt::getBM(
  attributes = c(platform_attr, "hgnc_symbol", "ensembl_gene_id"),
  mart = mart,
  uniqueRows = TRUE
)

# unique gene symbols (biomaRt)
biomart_symbols <- unique(res$hgnc_symbol)
length(biomart_symbols)

# [1] 39536

# unique gene symbols (biomart ensgene -> annotables)
biomart_ensgenes <- unique(res$ensembl_gene_id)
annot_symbols <- unique(grch38$symbol[match(biomart_ensgenes, grch38$ensgene)])
length(annot_symbols)

# [1] 57471
```
