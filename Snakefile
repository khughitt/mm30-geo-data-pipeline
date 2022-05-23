"""
MM25 GEO Data Preparation Pipeline
"""
import os
import yaml

# temp/dev (jan 26, 2022)
out_dir = "/data/proj/mm25/4.1/geo"

configfile: "config/config.yml"

# accessions = ["GSE134598", "GSE24080",
              # "GSE26760", "GSE2912", , "GSE57317",
              # "GSE5900", "GSE6477", "GSE6691", "GSE68871", "GSE7039", "GSE83503",
              # "GSE9782"]

# datasets for which expression data/metadata must be downloaded separately:
# - GSE118900 (expr data stored in supplemental file)
# - GSE106218 (expr data and metadata stored in supplemental files)
# - GSE7039 (metadata provided by author)
# - GSE24080 (metadata retrieved from GEO ftp)
# - GSE134598 (expr data stored in supplemental file)
# - GSE2912 (metadata retrieved from manuscript)
# - GSE26760 (metadata retrieved from broad/external website)

accessions = ["GSE39754", "GSE178340", "GSE106218", "GSE117847", "GSE19784", "GSE31161",
              "GSE158387", "GSE162205", "GSE83503", "GSE7039", "GSE118900",
              "GSE128251", "GSE13591", "GSE14519", "GSE16791", "GSE19554", "GSE47552",
              "GSE57317"]

rule all:
    input:
        expand(os.path.join(out_dir, "final/{acc}/data.csv"), acc=accessions),
        expand(os.path.join(out_dir, "final/{acc}/row-metadata.csv"), acc=accessions),
        expand(os.path.join(out_dir, "final/{acc}/column-metadata.csv"), acc=accessions),

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

