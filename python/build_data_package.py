#!/bin/env python
"""
Generate datapackages for finalized datasets
"""
import datetime
import os
import yaml
from frictionless import describe

snek = snakemake

# load dataset metadata
with open(snek.input[3]) as fp:
    mdata = yaml.load(fp, Loader=yaml.FullLoader)

# add entry to provenance section
mdata["provenance"].append({
    "action": "package",
    "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%m:%s"),
    "urls": ["https://github.com/khughitt/mm30-geo-data-pipeline"],
    "description": "Refer to manuscript / github url above for description of post-processing and packaging steps."
})

# switch to data dir
out_dir = os.path.dirname(snek.input[0])
os.chdir(out_dir)

pkg = describe("*.feather", stats=True)

# add dataset metadata
for key in mdata.keys():
    # skip "profile" for now (frictionless sets to "datapackage")
    if key == 'profile':
        continue

    pkg.custom[key] = mdata[key]

# write datapackage out
pkg.to_yaml("datapackage.yml")
