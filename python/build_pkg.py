#!/bin/env python
"""
Generate datapackages for finalized datasets
"""
import frictionless
import os
import yaml

# load dataset metadata
with open(snakemake.input[3]) as fp:
    mdata = yaml.load(fp, Loader=yaml.FullLoader)

# switch to data dir
out_dir = os.path.dirname(snakemake.input[0])
os.chdir(out_dir)

pkg = frictionless.describe("*.csv")

# add dataset metadata
for key in mdata.keys():
    # skip "profile" for now (frictionless sets to "datapackage")
    if key == 'profile':
        continue

    pkg.custom[key] = mdata[key]

# write datapackage out
pkg.to_yaml("datapackage.yml")
