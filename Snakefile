"""
MM25 GEO Data Preparation Pipeline
"""

# accessions = ["GSE106218", "GSE117847", "GSE118900", "GSE128251", "GSE134598",
              # "GSE13591", "GSE14519", "GSE16791", "GSE19554", "GSE19784", "GSE24080",
              # "GSE26760", "GSE2912", "GSE31161", "GSE39754", "GSE47552", "GSE57317",
              # "GSE5900", "GSE6477", "GSE6691", "GSE68871", "GSE7039", "GSE83503",
              # "GSE9782"]

# accessions for data which can be downloaded using the generalized version of the
# download script
accessions = ["GSE117847"]

rule all:
    input:
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/data.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/row-metadata.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/column-metadata.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/non-redundant/{acc}/datapackage.json", acc=accessions),

# rule all:
#     input: "/data/metadata.csv"

# rule all:
#     input:
#         "/data/proj/mm25/4.1/geo/reprocessed/GSE117847/datapackage.json"

rule download_general:
    output:
        "/data/proj/mm25/4.1/geo/cache/{acc}/{acc}_series_matrix.txt.gz",
        "/data/proj/mm25/4.1/geo/orig/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/{acc}/datapackage.json"
    params:
        accession="{acc}"
    script: "R/download.R"

rule reprocess_GSE117847:
    input:
        "/data/proj/mm25/4.1/geo/orig/GSE117847/data.csv",
        "/data/proj/mm25/4.1/geo/orig/GSE117847/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/GSE117847/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/orig/GSE117847/datapackage.json"
    output:
        "/data/proj/mm25/4.1/geo/reprocessed/GSE117847/data.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/GSE117847/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/GSE117847/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/reprocessed/GSE117847/datapackage.json"
    script: "R/process/GSE117847.R"

# generates a non-redundant version of the datasets by averaging expression for
# genes that appear multple times.
rule create_non_redundant:
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
        accession="{acc}"
    script: "R/create_nonredundant.R"

    # output:
    #     "/data/proj/mm25/4.1/geo/orig/GSE117847/GSE117847_series_matrix.txt.gz",
    #     "/data/proj/mm25/4.1/geo/clean/GSE117847/GSE117847_gene_expr.feather",
    #     "/data/proj/mm25/4.1/geo/clean/GSE117847/GSE117847_gene_expr_nr.feather",
    #     "/data/proj/mm25/4.1/geo/clean/GSE117847/GSE117847_sample_metadata.csv"

# rule download_GSE117847:
#     output:
#         "/data/cache/GSE117847/GSE117847_series_matrix.txt.gz",
#         "/data/orig/GSE117847/GSE117847.csv",
#         "/data/orig/GSE117847/GSE117847_pdata.csv",
#         "/data/orig/GSE117847/GSE117847_fdata.csv"
#     script: "R/download/GSE117847.R"


# rule build_metadata:
#     input:
#         rules.process_GSE106218.output,
#         rules.process_GSE117847.output,
#         rules.process_GSE118900.output,
#         rules.process_GSE128251.output,
#         rules.process_GSE134598.output,
#         rules.process_GSE13591.output,
#         rules.process_GSE14519.output,
#         rules.process_GSE16791.output,
#         rules.process_GSE19554.output,
#         rules.process_GSE19784.output,
#         rules.process_GSE24080.output,
#         rules.process_GSE26760.output,
#         rules.process_GSE2912.output,
#         rules.process_GSE31161.output,
#         rules.process_GSE39754.output,
#         rules.process_GSE47552.output,
#         rules.process_GSE57317.output,
#         rules.process_GSE5900.output,
#         rules.process_GSE6477.output,
#         rules.process_GSE6691.output,
#         rules.process_GSE68871.output,
#         rules.process_GSE7039.output,
#         rules.process_GSE83503.output,
#         rules.process_GSE9782.output
#     output: "/data/metadata.csv"
#     script: "R/build_metadata.R"
#
