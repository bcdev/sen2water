# -*- coding: utf-8 -*-

"""Combines Idepix flags and additional flags and selected bands into pixel class"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

import xarray as xr
import numpy as np
import dask.array as da

from sen2water.eoutils.eologging import logger

ACOLITE_BANDS_S2A = [
    "rhos_443",
    "rhos_492",
    "rhos_560",
    "rhos_665",
    "rhos_704",
    "rhos_740",
    "rhos_783",
    "rhos_833",
    "rhos_865",
    "rhos_1614",
    "rhos_2202",
]
ACOLITE_BANDS_S2B = [
    "rhos_442",
    "rhos_492",
    "rhos_559",
    "rhos_665",
    "rhos_704",
    "rhos_739",
    "rhos_780",
    "rhos_833",
    "rhos_864",
    "rhos_1610",
    "rhos_2186",
]
ACOLITE_BANDS_S2C = [
    "rhos_444",
    "rhos_489",
    "rhos_561",
    "rhos_667",
    "rhos_707",
    "rhos_741",
    "rhos_785",
    "rhos_835",
    "rhos_866",
    "rhos_1612",
    "rhos_2191",
]

class AcoliteDataset(xr.Dataset):

    __slots__ = "acolite", "dummy_band"

    def __init__(self, acolite):
        self.acolite = acolite
        self.dummy_band = None

    def __get__(self, instance, owner):
        return self.acolite.get(instance, owner)

    def __getitem__(self, key):
        if key in self.acolite.variables:
            return self.acolite.__getitem__(key)
        else:
            if self.dummy_band is None:
                logger.warning(f"adding all-nan band {key} to ACOLITE result")
                for some_band in self.acolite.variables:
                    if some_band.startswith("rhos_"):
                        break
                some_variable = self.acolite[some_band]
                dummy_data = da.empty(some_variable.data.shape, chunks=some_variable.data.chunks, dtype=some_variable.data.dtype)
                dummy_data[:] = np.nan
                attrs = {"wavelength": key[len("rhos_"):]}
                self.dummy_band = xr.DataArray(dummy_data, dims=some_variable.dims, attrs=attrs)
            return self.dummy_band

    def acolite_bands(self):
        return ACOLITE_BANDS_S2A if ACOLITE_BANDS_S2A[0] in self.acolite.variables \
            else ACOLITE_BANDS_S2B if ACOLITE_BANDS_S2B[0] in self.acolite.variables \
            else ACOLITE_BANDS_S2C