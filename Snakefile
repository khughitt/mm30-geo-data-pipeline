"""
MM30 GEO Data Preparation Pipeline
"""
import os
import yaml

configfile: "config/config.yml"

out_dir = config["output_dir"]

accessions = ['GSE106218', 'GSE117847', 'GSE118900', 'GSE128251', 'GSE134598',
              'GSE13591', 'GSE144249', 'GSE14519', 'GSE158387', 'GSE162205', 'GSE16791', 
              'GSE178340', 'GSE193531', 'GSE19554', 'GSE19784', 'GSE24080', 'GSE26760', 'GSE2912',
              'GSE31161','GSE39754', 'GSE47552', 'GSE57317', 'GSE5900', 'GSE6477', 'GSE6691',
              'GSE68871', 'GSE7039', 'GSE83503', 'GSE9782']

rule all:
    input:
        expand(os.path.join(out_dir, "final/{acc}/datapackage.yml"), acc=accessions),
        os.path.join(out_dir, "metadata.tsv")

rule build_metadata:
    output:
        os.path.join(out_dir, "metadata.tsv")
    script:
        "R/build_metadata.R"

rule build_data_package:
    input:
        os.path.join(out_dir, "final/{acc}/data.csv"),
        os.path.join(out_dir, "final/{acc}/row-metadata.csv"),
        os.path.join(out_dir, "final/{acc}/column-metadata.csv"),
        os.path.join("metadata", "{acc}.yml")
    output:
        os.path.join(out_dir, "final/{acc}/datapackage.yml"),
    script: "python/build_data_package.py"

rule finalize:
    input: 
        os.path.join(out_dir, "normalized/{acc}/data.csv"),
        os.path.join(out_dir, "normalized/{acc}/row-metadata.csv"),
        os.path.join(out_dir, "normalized/{acc}/column-metadata.csv"),
    output:
        os.path.join(out_dir, "final/{acc}/data.csv"),
        os.path.join(out_dir, "final/{acc}/row-metadata.csv"),
        os.path.join(out_dir, "final/{acc}/column-metadata.csv"),
    script: "R/finalize.R"

rule normalize:
    input:
        os.path.join(out_dir, "original/{acc}/data.csv"),
        os.path.join(out_dir, "original/{acc}/row-metadata.csv"),
        os.path.join(out_dir, "original/{acc}/column-metadata.csv"),
    output:
        os.path.join(out_dir, "normalized/{acc}/data.csv"),
        os.path.join(out_dir, "normalized/{acc}/row-metadata.csv"),
        os.path.join(out_dir, "normalized/{acc}/column-metadata.csv"),
    script: "R/normalize/{wildcards.acc}.R"

rule download:
    output:
        os.path.join(out_dir, "original/{acc}/data.csv"),
        os.path.join(out_dir, "original/{acc}/row-metadata.csv"),
        os.path.join(out_dir, "original/{acc}/column-metadata.csv"),
    script: "R/download/{wildcards.acc}.R"
