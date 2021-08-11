"""
MM25 GEO Data Preparation Pipeline
"""
rule all:
    input: "/data/metadata.tsv"

rule GSE106218:
    output:
        "/data/raw/GSE106218/GSE106218_series_matrix.txt.gz",
        "/data/clean/GSE106218/GSE106218_gene_expr.feather",
        "/data/clean/GSE106218/GSE106218_gene_expr_nr.feather",
        "/data/clean/GSE106218/GSE106218_sample_metadata.tsv"
    script: "R/GSE106218.R"

rule GSE117847:
    output:
        "/data/raw/GSE117847/GSE117847_series_matrix.txt.gz",
        "/data/clean/GSE117847/GSE117847_gene_expr.feather",
        "/data/clean/GSE117847/GSE117847_gene_expr_nr.feather",
        "/data/clean/GSE117847/GSE117847_sample_metadata.tsv"
    script: "R/GSE117847.R"

rule GSE118900:
    output:
        "/data/raw/GSE118900/GSE118900_series_matrix.txt.gz",
        "/data/clean/GSE118900/GSE118900_gene_expr.feather",
        "/data/clean/GSE118900/GSE118900_gene_expr_nr.feather",
        "/data/clean/GSE118900/GSE118900_sample_metadata.tsv"
    script: "R/GSE118900.R"

rule GSE128251:
    output:
        "/data/raw/GSE128251/GSE128251_series_matrix.txt.gz",
        "/data/clean/GSE128251/GSE128251_gene_expr.feather",
        "/data/clean/GSE128251/GSE128251_gene_expr_nr.feather",
        "/data/clean/GSE128251/GSE128251_sample_metadata.tsv"
    script: "R/GSE128251.R"

rule GSE134598:
    output:
        "/data/raw/GSE134598/GSE134598_series_matrix.txt.gz",
        "/data/clean/GSE134598/GSE134598_gene_expr.feather",
        "/data/clean/GSE134598/GSE134598_gene_expr_nr.feather",
        "/data/clean/GSE134598/GSE134598_sample_metadata.tsv"
    script: "R/GSE134598.R"

rule GSE13591:
    output:
        "/data/raw/GSE13591/GSE13591_series_matrix.txt.gz",
        "/data/clean/GSE13591/GSE13591_gene_expr.feather",
        "/data/clean/GSE13591/GSE13591_gene_expr_nr.feather",
        "/data/clean/GSE13591/GSE13591_sample_metadata.tsv"
    script: "R/GSE13591.R"

rule GSE14519:
    output:
        "/data/raw/GSE14519/GSE14519_series_matrix.txt.gz",
        "/data/clean/GSE14519/GSE14519_gene_expr.feather",
        "/data/clean/GSE14519/GSE14519_gene_expr_nr.feather",
        "/data/clean/GSE14519/GSE14519_sample_metadata.tsv"
    script: "R/GSE14519.R"

rule GSE16791:
    output:
        "/data/raw/GSE16791/GSE16791_series_matrix.txt.gz",
        "/data/clean/GSE16791/GSE16791_gene_expr.feather",
        "/data/clean/GSE16791/GSE16791_gene_expr_nr.feather",
        "/data/clean/GSE16791/GSE16791_sample_metadata.tsv"
    script: "R/GSE16791.R"

rule GSE19554:
    output:
        "/data/raw/GSE19554/GSE19554_series_matrix.txt.gz",
        "/data/clean/GSE19554/GSE19554_gene_expr.feather",
        "/data/clean/GSE19554/GSE19554_gene_expr_nr.feather",
        "/data/clean/GSE19554/GSE19554_sample_metadata.tsv"
    script: "R/GSE19554.R"

rule GSE19784:
    output:
        "/data/raw/GSE19784/GSE19784_series_matrix.txt.gz",
        "/data/clean/GSE19784/GSE19784_gene_expr.feather",
        "/data/clean/GSE19784/GSE19784_gene_expr_nr.feather",
        "/data/clean/GSE19784/GSE19784_sample_metadata.tsv"
    script: "R/GSE19784.R"

rule GSE24080:
    output:
        "/data/raw/GSE24080/GSE24080_series_matrix.txt.gz",
        "/data/clean/GSE24080/GSE24080_gene_expr.feather",
        "/data/clean/GSE24080/GSE24080_gene_expr_nr.feather",
        "/data/clean/GSE24080/GSE24080_sample_metadata.tsv"
    script: "R/GSE24080.R"

rule GSE26760:
    output:
        "/data/raw/GSE26760/GSE26760_series_matrix.txt.gz",
        "/data/clean/GSE26760/GSE26760_gene_expr.feather",
        "/data/clean/GSE26760/GSE26760_gene_expr_nr.feather",
        "/data/clean/GSE26760/GSE26760_sample_metadata.tsv"
    script: "R/GSE26760.R"

rule GSE2912:
    output:
        "/data/raw/GSE2912/GSE2912_series_matrix.txt.gz",
        "/data/clean/GSE2912/GSE2912_gene_expr.feather",
        "/data/clean/GSE2912/GSE2912_gene_expr_nr.feather",
        "/data/clean/GSE2912/GSE2912_sample_metadata.tsv"
    script: "R/GSE2912.R"

rule GSE31161:
    output:
        "/data/raw/GSE31161/GSE31161_series_matrix.txt.gz",
        "/data/clean/GSE31161/GSE31161_gene_expr.feather",
        "/data/clean/GSE31161/GSE31161_gene_expr_nr.feather",
        "/data/clean/GSE31161/GSE31161_sample_metadata.tsv"
    script: "R/GSE31161.R"

rule GSE39754:
    output:
        "/data/raw/GSE39754/GSE39754_series_matrix.txt.gz",
        "/data/clean/GSE39754/GSE39754_gene_expr.feather",
        "/data/clean/GSE39754/GSE39754_gene_expr_nr.feather",
        "/data/clean/GSE39754/GSE39754_sample_metadata.tsv"
    script: "R/GSE39754.R"

rule GSE47552:
    output:
        "/data/raw/GSE47552/GSE47552_series_matrix.txt.gz",
        "/data/clean/GSE47552/GSE47552_gene_expr.feather",
        "/data/clean/GSE47552/GSE47552_gene_expr_nr.feather",
        "/data/clean/GSE47552/GSE47552_sample_metadata.tsv"
    script: "R/GSE47552.R"

rule GSE57317:
    output:
        "/data/raw/GSE57317/GSE57317_series_matrix.txt.gz",
        "/data/clean/GSE57317/GSE57317_gene_expr.feather",
        "/data/clean/GSE57317/GSE57317_gene_expr_nr.feather",
        "/data/clean/GSE57317/GSE57317_sample_metadata.tsv"
    script: "R/GSE57317.R"

rule GSE5900:
    output:
        "/data/raw/GSE5900/GSE5900_series_matrix.txt.gz",
        "/data/clean/GSE5900/GSE5900_gene_expr.feather",
        "/data/clean/GSE5900/GSE5900_gene_expr_nr.feather",
        "/data/clean/GSE5900/GSE5900_sample_metadata.tsv"
    script: "R/GSE5900.R"

rule GSE6477:
    output:
        "/data/raw/GSE6477/GSE6477_series_matrix.txt.gz",
        "/data/clean/GSE6477/GSE6477_gene_expr.feather",
        "/data/clean/GSE6477/GSE6477_gene_expr_nr.feather",
        "/data/clean/GSE6477/GSE6477_sample_metadata.tsv"
    script: "R/GSE6477.R"

rule GSE6691:
    output:
        "/data/raw/GSE6691/GSE6691_series_matrix.txt.gz",
        "/data/clean/GSE6691/GSE6691_gene_expr.feather",
        "/data/clean/GSE6691/GSE6691_gene_expr_nr.feather",
        "/data/clean/GSE6691/GSE6691_sample_metadata.tsv"
    script: "R/GSE6691.R"

rule GSE68871:
    output:
        "/data/raw/GSE68871/GSE68871_series_matrix.txt.gz",
        "/data/clean/GSE68871/GSE68871_gene_expr.feather",
        "/data/clean/GSE68871/GSE68871_gene_expr_nr.feather",
        "/data/clean/GSE68871/GSE68871_sample_metadata.tsv"
    script: "R/GSE68871.R"

rule GSE7039:
    output:
        "/data/raw/GSE7039/GSE7039-GPL4819_series_matrix.txt.gz",
        "/data/raw/GSE7039/GSE7039-GPL4820_series_matrix.txt.gz",
        "/data/clean/GSE7039/GSE7039_gene_expr.feather",
        "/data/clean/GSE7039/GSE7039_gene_expr_nr.feather",
        "/data/clean/GSE7039/GSE7039_sample_metadata.tsv"
    script: "R/GSE7039.R"

rule GSE83503:
    output:
        "/data/raw/GSE83503/GSE83503_series_matrix.txt.gz",
        "/data/clean/GSE83503/GSE83503_gene_expr.feather",
        "/data/clean/GSE83503/GSE83503_gene_expr_nr.feather",
        "/data/clean/GSE83503/GSE83503_sample_metadata.tsv"
    script: "R/GSE83503.R"

rule GSE9782:
    output:
        "/data/raw/GSE9782/GSE9782-GPL96_series_matrix.txt.gz",
        "/data/raw/GSE9782/GSE9782-GPL97_series_matrix.txt.gz",
        "/data/clean/GSE9782/GSE9782_gene_expr.feather",
        "/data/clean/GSE9782/GSE9782_gene_expr_nr.feather",
        "/data/clean/GSE9782/GSE9782_sample_metadata.tsv"
    script: "R/GSE9782.R"

rule build_metadata:
    input: 
        rules.GSE106218.output,
        rules.GSE117847.output,
        rules.GSE118900.output,
        rules.GSE128251.output,
        rules.GSE134598.output,
        rules.GSE13591.output,
        rules.GSE14519.output,
        rules.GSE16791.output,
        rules.GSE19554.output,
        rules.GSE19784.output,
        rules.GSE24080.output,
        rules.GSE26760.output,
        rules.GSE2912.output,
        rules.GSE31161.output,
        rules.GSE39754.output,
        rules.GSE47552.output,
        rules.GSE57317.output,
        rules.GSE5900.output,
        rules.GSE6477.output,
        rules.GSE6691.output,
        rules.GSE68871.output,
        rules.GSE7039.output,
        rules.GSE83503.output,
        rules.GSE9782.output
    output: "/data/metadata.tsv"
    script: "R/metadata.R"

