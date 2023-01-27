MM30 GEO Data Preparation Pipeline
==================================

Overview
--------

This repository contains a [Snakemake](https://snakemake.readthedocs.io/en/stable/)
pipeline for downloading and processing a set of publicly-available multiple myeloma GEO
transcriptomics datasets.

Data are retrieved using the
[GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) package
for [Bioconductor](https://bioconductor.org/).

To enable the results to be reproduced, a Dockerfile has been provided which includes
the specific versions of all software used in the original analysis.

Pre-requisites
--------------

The only software requirement to run the analyses contained in this repo is
[Docker](https://www.docker.com/).

Because some of the datasets included useful clinical metadata that was only
available via the original publication (or in one case, correspondance with the
author), a couple additional steps must be taken to first download and prepare
these datasets for use:

**Kuiper et al. (2012)**

Download [Supplementary Document A](https://www.nature.com/articles/leu2012127#Sec14)
from Kuiper _et al._ (2012).

Next, extract the second sheet ("Survival HOVON-65GMMGHD4") from the xls file 
and save it as a csv file.

For example, using the `in2csv` command from
[csvkit](https://csvkit.readthedocs.io/en/latest/scripts/in2csv.html), one could do:

```
in2csv 41375_2012_BFleu2012127_MOESM30_ESM.xls \
  --sheet "Survival HOVON-65GMMGHD4" \
  --skip-lines 2 \
  > kuiper2012_supp_patient_survival.csv
```

**Chapman et al. (2011)**

Download the `mmrc.sample.information.xls` file from the [Broad Multiple Myeloma Genomics Portal
(MMGP)](http://portals.broadinstitute.org/mmgp/data/browseData?conversationPropagation=begin).

The single sheet contained in the xls file should similarly be extracted and converted
to csv.

In this case, `in2csv` cannot be used due to the fact that xls files do not have an
integer type, leading `in2csv` to convert some fields like "age" to floats unneccessry.

Instead, one can use Office or Libreoffice to manually extract the sheet to csv.

Place the generated csv file in the `supp/clean` folder, along-side of the Kuiper _et
al._ supplemental data file.

**Agnelli et al. (2005)**

Appendex A ("Clinical Details of the 50 MM Patients Enrolled in the Proprietary
Database") from [Molecular Classification of Multiple Myeloma: A Distinct
Transcriptional Profile Characterizes Patients Expressing CCND1 and Negative for 14q32
Translocations (Agnelli et al., 2005)](https://ascopubs.org/doi/10.1200/JCO.2005.01.3870) 
saved as csv:

```
Patient,Marrow Plasmacytosis (%),Sex,Monoclonal Component,Age (years),Stage,β2 Microglobulin,Bone Lesions
MM-004,15,F,Gκ,54,IA,P,–
...
```

**Decaux et al. (2007)**

(Currently not available online; recieved via correspondance with author)

Usage
-----

To download and process the datasets as was performed in the manuscript, clone the repo

```
git clone https://github.com/khughitt/geo-data-prep
cd geo-data-prep
```

Copy the needed manuscript supplemental files described in the step above to the
"supp/clean/" folder.

The result should look like:

```
supp/clean/
├── GSE24080_MM_UAMS565_ClinInfo_27Jun2008_LS_clean.tsv
├── GSE7039_MM_Survival_time.csv
├── kuiper2012_supp_patient_survival.csv
└── mmrc.sample.information.csv
```

Next, build the included Dockerfile by running:

```
docker build . -t mm30-geo
```

Decide on a location on your _host_ system where you wish to store the processed
datasets and create the directory if is doesn't already exist.

Next, use `docker run`  to start a container with, using the volume (`-v`) switch to specify the directory
where you wish to have the raw and processed GEO data saved to.

For example, if you would like to have the data saved to "/data/geo" on the host system, run:

```
docker run -v /data/geo:/data -it mm30-geo
```

This will start a container and open a bash shell.

Finally, to download and process the GEO datasets, run `snakemake`, specifying the
number of simultaneous jobs to run, e.g.:

```
snakemake -j 4
```

Outputs
-------

Running the data preparation pipeline will result in "raw" and "clean" versions of each
dataset being generated in the specified data folder on the host system, e.g.:

```
/data/geo/raw/GSE106218/
├── GPL16791.soft
├── GSE106218
│   ├── GSE106218_GEO_clinical_information_MM.txt.gz
│   └── GSE106218_GEO_processed_MM_raw_TPM_matrix.txt.gz
└── GSE106218_series_matrix.txt.gz

/data/geo/clean/GSE106218/
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

Trouble-shooting
----------------

Occasionally the `getGEO()` may fail due to a network error, e.g.:

```
Error in download.file(myurl, destfile, mode = mode, quiet = TRUE, method = getOption("download.file.method.GEOquery")) :
  download from 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL570&form=text&view=full'failed
Calls: getGEO ... parseGSEMatrix -> getGEO -> getGEOfile -> download.file
In addition: Warning message:
In download.file(myurl, destfile, mode = mode, quiet = TRUE, method = getOption("download.file.method.GEOquery")) :
  URL 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL570&form=text&view=full': Timeout of 60 seconds was reached
Execution halted
```

This usually indicates a temporary network issue. Snakemake is able to detect which jobs
completeted successfully, and which ones failed, so simply re-running the pipeline
is sufficient to pick up where you left off.

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
R script.
