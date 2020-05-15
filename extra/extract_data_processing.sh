#!/bin/env bash
#
# Extracts GEO dataset data processing information from series matrix files
#
BASE_DIR="/data/human/geo/1.1"

for x in $(ls $BASE_DIR/GSE*/raw/GSE*_series_matrix.txt.gz); do
    outfile=$(basename $x)
    outfile=${outfile/_series_matrix.txt.gz/_data_processing.txt}

    zgrep "^!Sample_data_processing" $x > $outfile
done
