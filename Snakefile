"""
MM25 GEO Data Preparation Pipeline
"""

# accessions = ["GSE106218", "GSE117847", "GSE118900", "GSE128251", "GSE134598",
              # "GSE13591", "GSE14519", "GSE16791", "GSE19554", "GSE24080",
              # "GSE26760", "GSE2912", "GSE31161", "GSE39754", "GSE47552", "GSE57317",
              # "GSE5900", "GSE6477", "GSE6691", "GSE68871", "GSE7039", "GSE83503",
              # "GSE9782"]

# accessions for data which can be downloaded using the generalized version of the
# download script
accessions = ["GSE117847", "GSE19784"]

rule all:
    input:
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/data.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/row-metadata.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/column-metadata.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/datapackage.json", acc=accessions),

# generates a "non-redundant" version of the datasets by averaging expression for
# genes that appear multple times.
rule finalize:
    input: 
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/datapackage.json"
    output:
        "/data/proj/mm25/4.1/geo/non-redundant/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/non-redundant/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/non-redundant/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/non-redundant/{acc}/datapackage.json"
    params:
        scale_pca=True
    script: "R/create_nonredundant.R"

rule normalize:
    input:
        "/data/proj/mm25/4.1/geo/orig/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/datapackage.json"
    output:
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/{acc}/datapackage.json"
    script: "R/prepare-data/{wildcards.acc}.R"


rule download:
    output:
        "/data/proj/mm25/4.1/geo/cache/{acc}/{acc}_series_matrix.txt.gz",
        "/data/proj/mm25/4.1/geo/orig/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/datapackage.json"
    params:
        accession="{acc}"
    script: "R/download.R"

