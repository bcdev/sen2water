#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import json
import pyximport

from polymer.copernicus_dem import CopernicusDEM
pyximport.install()
from polymer.level1_netcdf import Level1_NETCDF
from polymer.level2 import Level2
from polymer.main import run_atm_corr

def run():
    with open(sys.argv[1]) as f: parameters = json.load(f)
    if sys.argv[2] == "nasa":
        l1 = Level1_NETCDF(sys.argv[-2])
    else:
        l1 = Level1_NETCDF(sys.argv[-2],
                           ancillary='ECMWFT',
                           altitude=CopernicusDEM(directory=sys.argv[3]))
    l2 = Level2(fmt="netcdf4", filename=sys.argv[-1])
    run_atm_corr(l1, l2, **parameters)

if __name__ == "__main__":
    run()
