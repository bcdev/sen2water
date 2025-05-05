#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Converts Jpeg2000 B01 into dummy granule with geo-coding"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "proprietary, TBD"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import sys
import rioxarray as rio
import xarray as xr
import numpy as np

input = sys.argv[1]
granule = input[input.rfind('/')+2:input.rfind('/')+7]
#da = rio.open_rasterio('/windows/tmp/sen2water/S2A_MSIL1C_20230611T103631_N0509_R008_T32UME_20230611T142059.SAFE/GRANULE/L1C_T32UME_A041618_20230611T103626/IMG_DATA/T32UME_20230611T103631_B01.jp2')
da = rio.open_rasterio(input)
da2 = da.isel({'band': 0}, drop=True)
mask = xr.zeros_like(da2, dtype=np.uint8)
#mask.rio.to_raster('32UME-zeros.tif', tiled=True, compress='LZW', blockxsize=1024, blockysize=1024)
mask.rio.to_raster(f'{granule}-zeros.tif', tiled=True, compress='LZW', blockxsize=512, blockysize=512)
print(f"{granule}-zeros.tif")
