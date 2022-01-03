"""
MM25 GEO Data Preparation Pipeline
"""

# accessions = ["GSE106218", "GSE117847", "GSE118900", "GSE128251", "GSE134598",
              # "GSE13591", "GSE14519", "GSE16791", "GSE19554", "GSE24080",
              # "GSE26760", "GSE2912", "GSE31161", "GSE39754", "GSE47552", "GSE57317",
              # "GSE5900", "GSE6477", "GSE6691", "GSE68871", "GSE7039", "GSE83503",
              # "GSE9782"]

# datasets for which biomart is queried to retrieve updated probe -> gene mappings
# biomart_annots = ["GSE19784", "GSE31161", "GSE57317", "GSE6477", "GSE24080",
#         "GSE9782", "GSE5900", "GSE6691", "GSE68871", "GSE128251", "GSE14519",
#         "GSE19554", "GSE47552", "GSE2912", "GSE26760", "GSE13591", "GSE16791"]

# microarray datasets for which no biomart mappings are needed or available:
# "GSE117847" "GSE39754"  "GSE7039"   "GSE83503"

# datasets for which expression data/metadata must be downloaded separately:
# - GSE19784 (metadata retrieved from manuscript)
# - GSE118900 (expr data stored in supplemental file)
# - GSE106218 (expr data and metadata stored in supplemental files)
# - GSE7039 (metadata provided by author)
# - GSE24080 (metadata retrieved from GEO ftp)
# - GSE134598 (expr data stored in supplemental file)
# - GSE2912 (metadata retrieved from manuscript)
# - GSE26760 (metadata retrieved from broad/external website)

# accessions for data which can be downloaded using the generalized version of the
# download script
accessions = ["GSE117847", "GSE19784", "GSE31161"]

# microarray datasets; includes an extra step to retrieve a probe -> gene mapping
# microarray_datasets = ['GSE117847', 'GSE128251', 'GSE13591', 'GSE14519', 'GSE16791',
#         'GSE19554', 'GSE19784', 'GSE24080', 'GSE26760', 'GSE2912', 'GSE31161',
#         'GSE39754', 'GSE47552', 'GSE57317', 'GSE5900', 'GSE6477', 'GSE6691', 'GSE68871',
#         'GSE7039', 'GSE83503', 'GSE9782']

# geo microarray platform ids
# geo_platforms = ["GPL96", "GPL97", "GPL570", "GPL10558", "GPL5175", "GPL6244"]

rule all:
    input:
        expand("/data/proj/mm25/4.1/geo/final/{acc}/data.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/final/{acc}/row-metadata.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/final/{acc}/column-metadata.csv", acc=accessions),
        expand("/data/proj/mm25/4.1/geo/final/{acc}/datapackage.json", acc=accessions)
        # expand("/data/proj/mmm25/4.1/geo/annot/{platform_id}.csv", platform_id=geo_platforms)

# rule get_biomart_mappings:
#     output:
#         "/data/proj/mmm25/4.1/geo/annot/{platform_id}.csv"
#     params:
#         platform_id="{platform_id}",
#         ensembl_version=105
#     script:
#         "R/get_biomart_mapping.R"

# generates a "final" version of the datasets by averaging expression for genes that
# appear multple times.
rule finalize:
    input: 
        "/data/proj/mm25/4.1/geo/normalized/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/normalized/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/normalized/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/normalized/{acc}/datapackage.json"
    output:
        "/data/proj/mm25/4.1/geo/final/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/final/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/final/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/final/{acc}/datapackage.json"
    params:
        scale_pca=True
    script: "R/finalize.R"

rule normalize:
    input:
        "/data/proj/mm25/4.1/geo/original/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/original/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/original/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/original/{acc}/datapackage.json"
    output:
        "/data/proj/mm25/4.1/geo/normalized/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/normalized/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/normalized/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/normalized/{acc}/datapackage.json"
    script: "R/normalize/{wildcards.acc}.R"


rule download:
    output:
        "/data/proj/mm25/4.1/geo/cache/{acc}/{acc}_series_matrix.txt.gz",
        "/data/proj/mm25/4.1/geo/original/{acc}/data.csv",
        "/data/proj/mm25/4.1/geo/original/{acc}/row-metadata.csv",
        "/data/proj/mm25/4.1/geo/original/{acc}/column-metadata.csv",
        "/data/proj/mm25/4.1/geo/original/{acc}/datapackage.json"
    params:
        accession="{acc}"
    script: "R/download.R"

