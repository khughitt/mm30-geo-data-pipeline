GEO Data Preparation Scripts
============================

Overview
--------

This repository contains a collection of R scripts to download and process a set of
publicly-available multiple myeloma GEO transcriptomics datasets.

Data are retrieved using the
[GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) package
for Bioconductor.

To enable the results to be reproduced, a Dockerfile and correspond `renv.lock` file
have been included, allowing one to construct a Docker image with the precise versions
of software used in preparation of the manuscript.

Usage
-----

To download and process the datasets as was performed in the manuscript, clone the repo and build the
included Dockerfile:

```
git clone https://github.com/khughitt/geo-data-prep
cd geo-data-prep

docker build .
```

Decide on a location on your _host_ system where you wish to store the processed
datasets and create the directory if is doesn't already exist.

Next, start a container with:

```
docker run -v /path/to/data/dir:/data -it --name geo-data <image id>
```

Where "/path/to/data/dir" refers to the directory you wish to store outputs in, and
"<image id>" is the docker image id generated during the build step.

This will start a container and open a bash shell.

Finally, to download and process the GEO datasets, use the run the `prepare_data.sh`
helper script located in the `extra/` folder:

```
./extra/prepare_data.sh
```

Outputs
-------

Running the data preparation script will result in "raw" and "clean" versions of each
dataset being generated in the specified data folder on the host system, e.g.:

```
/data/raw/GSE106218/
├── GPL16791.soft
├── GSE106218
│   ├── GSE106218_GEO_clinical_information_MM.txt.gz
│   └── GSE106218_GEO_processed_MM_raw_TPM_matrix.txt.gz
└── GSE106218_series_matrix.txt.gz

/data/clean/GSE106218/
├── GSE106218_gene_expr.feather
├── GSE106218_gene_expr_nr.feather
└── GSE106218_sample_metadata.tsv
```

The files in the `raw/` directory are those downloaded by GEOquery and are cached by the
script to avoid needing to re-download in the future.

The `clean/` folder contains gene expression data and sample metadata stored in a
consistent manner.

The `GSEXXXX_gene_expr.feather` contains a gene expression count matrix, stored in the
[feather file format](https://github.com/wesm/feather), which has been size-factor (CPM)
normalized and had QC-related entries removed. The original identifiers (most often
microarray probe ids) are included, as well as mapped gene symbols, when possible.

The `_nr.feather` version ("non-redundant") has been further processed so that one entry
is included per gene symbol (median expression of probes used in cases where multiple
probes mapped to the same symbol), and only HUGO gene symbol identifiers are included.

The `GSEXX_sample_metadata.tsv` files contain one row per sample in the gene expression
matrices, with fields corresponding to the sample GEO accession, platform id, and one or
more additional clinical metadata fields related to treatment, survival, etc.

Notes
-----

While it would ideally be nice to generalize the functionality in the dataset-specific
scripts so that arbitrary GEO datasets could be processed this way, there was found to
be a large number of discrepancies in the quality and format of both the gene expression
data and metadata associated with each dataset.

In several cases, missing data needed to be handled by downloading additional
supplemental datasets available either via GEO, or in some cases, as supplemental files
associated with manuscripts.

As such, each dataset ultimately had to be handled separately in a dataset-specific
script.
