GEO Data Preparation Scripts
============================

A collection of R scripts to download process a set of GEO datasets in a normalized
manner.

Data are retrieved using the
[GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) package
for bioconductor.

While it would ideally be nice to generalize the functionality in the dataset-specific
scripts so that arbitrary GEO datasets could be processed this way, there was found to
be a large number of discrepencies in the quality and format of both the gene expression
data and metadata associated with each dataset, even when using GEOquery. In several
cases, missing data needed to be handled by downloading additional supplemental datasets
available either via GEO, or in some cases, as supplemental files associated with
manuscripts. As such, each dataset ultimately had to be handled separately in a
dataset-specific script.

For each given script, outputs of the following format are generated:

```
GSE2912
├── raw
│   ├── GSE2912_series_matrix.txt.gz
│   └── GPL96.annot.gz
└── processed
    ├── GSE2912_gene_expr_nr.feather
    ├── GSE2912_sample_metadata.tsv
    └── GSE2912_gene_expr.feather
```

The files in the `raw/` directory are those downloaded by GEOquery and are cached by the
script to avoid needing to re-download in the future.

The `processed/` folder contains gene expression data and sample metadata stored in a
consistent manner.

The `GSEXXXX_gene_expr.feather` contains a gene expression count matrix, stored in the
[feather file format](https://github.com/wesm/feather), which has been size-factor (CPM)
normalized and had QC-related entries removed. The original identifiers (most often
microarray probe ids) are included, as well as mapped gene symbols, when possible.

The `_nr.feather` version has been further processed so that one entry is included per
gene symbol (median expression of probes used in cases where multiple probes mapped to
the same symbol), and only HUGO gene symbol identifiers are included.

The `GSEXX_sample_metadata.tsv` files contain one row per sample in the gene expression
matrices, with fields corresponding to the sample GEO accession, platform id, and one or
more additional clinical metadata fields related to treatment, survival, etc.

