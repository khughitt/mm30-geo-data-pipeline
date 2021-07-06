#!/bin/bash
cd /geo/R

for x in GSE*.R; do
    Rscript $x;
done
