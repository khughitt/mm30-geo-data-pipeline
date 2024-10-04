"""
MM30 GEO Data Preparation Pipeline
"""
import os
import yaml

accessions = ['GSE106218', 'GSE117847', 'GSE118900', 'GSE128251', 'GSE134598',
              'GSE13591', 'GSE144249', 'GSE14519', 'GSE158387', 'GSE162205', 'GSE16791', 
              'GSE178340', 'GSE193531', 'GSE19554', 'GSE19784', 'GSE24080', 'GSE26760', 'GSE2912',
              'GSE31161','GSE39754', 'GSE47552', 'GSE57317', 'GSE5900', 'GSE6477', 'GSE6691',
              'GSE68871', 'GSE7039', 'GSE83503', 'GSE9782']

rule all:
    input:
        expand("/data/final/{acc}/datapackage.yml", acc=accessions),
        "/data/metadata.feather"

rule build_metadata:
    output:
        "/data/metadata.feather"
    script:
        "R/build_metadata.R"

rule build_data_package:
    input:
        "/data/final/{acc}/data.feather",
        "/data/final/{acc}/row-metadata.feather",
        "/data/final/{acc}/column-metadata.feather",
        "metadata/{acc}.yml"
    output:
        "/data/final/{acc}/datapackage.yml",
    script: "python/build_data_package.py"

rule finalize:
    input: 
        "/data/normalized/{acc}/data.feather",
        "/data/normalized/{acc}/row-metadata.feather",
        "/data/normalized/{acc}/column-metadata.feather",
    output:
        "/data/final/{acc}/data.feather",
        "/data/final/{acc}/row-metadata.feather",
        "/data/final/{acc}/column-metadata.feather",
    script: "R/finalize.R"

rule normalize:
    input:
        "/data/original/{acc}/data.feather",
        "/data/original/{acc}/row-metadata.feather",
        "/data/original/{acc}/column-metadata.feather",
    output:
        "/data/normalized/{acc}/data.feather",
        "/data/normalized/{acc}/row-metadata.feather",
        "/data/normalized/{acc}/column-metadata.feather",
    script: "R/normalize/{wildcards.acc}.R"

rule download:
    output:
        "/data/original/{acc}/data.feather",
        "/data/original/{acc}/row-metadata.feather",
        "/data/original/{acc}/column-metadata.feather",
    script: "R/download/{wildcards.acc}.R"
