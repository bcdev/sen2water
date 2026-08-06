# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Olaf Danne (Brockmann Consult GmbH)"
__copyright__ = "Copyright 2026, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"


from typing import Optional, List

import numpy as np
import xarray as xr
from eopf.computing import EOProcessingUnit, MappingDataType, MappingAuxiliary

from sen2water.msiidepix.constants import IdepixMsiConstants
from msiidepix.pixelclassification import PixelClassification


class IdepixClassificationPU(EOProcessingUnit):
    """
    Processing Unit for classification by threshold tests.
    See table 'Candidates for Processing Units' presented at PM1.
    """
    def run(
            self,
            inputs: MappingDataType,
            adfs: Optional[MappingAuxiliary] = None,
            mode: Optional[str] = None,
            *,
            resolution: int = 60,
            chunksize_in_meters=36600,
            ancillary: Optional[List[str]] = None,
            **kwargs,
    ) -> MappingDataType:

        # check input

        if not "l1c" in inputs:
            raise KeyError("No 'l1c' in input of ResamplingMainPU.")
        l1c: xr.DataTree = inputs["l1c"]
        if not isinstance(l1c, xr.DataTree):
            raise TypeError("Input 'l1c' of ResamplingMainPU is not an xarray.DataTree.")

        dims = {
            "y": l1c[f"measurements/reflectance/resampled"].coords.sizes["y"],
            "x": l1c[f"measurements/reflectance/resampled"].coords.sizes["x"]
        }

        toa_data = [l1c[f"measurements/reflectance/resampled/{band}"].data for band in IdepixMsiConstants.bands]
        pixel_classif = PixelClassification().apply(
            *toa_data,
            thresh_cw=0.007,
            thresh_gcl=-0.11,
            thresh_cl=0.007,
            dtype=np.int32,
            **kwargs
        )

        result = xr.DataTree()
        result["pixel_classif_flag"] = xr.DataArray(
            pixel_classif,
            dims=dims,
            attrs={"long_name": "pixel_classif_flag"}
        )

        return {"pixel_classif_flag": result}

