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

from msiidepix.slopeaspectorientation import SlopeAspectOrientation


class SlopeAspectOrientationPU(EOProcessingUnit):
    """
    Processing Unit for computation of Slope, Aspect and Orientation for S2 MSI..
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

        saa_data = l1c["conditions/geometry/resampled/saa"].data
        sza_data = l1c["conditions/geometry/resampled/sza"].data
        vaa_data = l1c["conditions/geometry/resampled/vaa_mean"].data
        vza_data = l1c["conditions/geometry/resampled/vza_mean"].data

        geometry_data = [saa_data, sza_data, vaa_data, vza_data]

        # TODO: we also need elevation

        sao = SlopeAspectOrientation().apply(
            *geometry_data,
            dtype=np.float32,
            **kwargs
        )

        result = xr.DataTree()
        result["slope_aspect_orientation"] = xr.DataArray(
            sao,
            dims=dims,
            attrs={"long_name": "slope_aspect_orientation"}
        )

        return {"slope_aspect_orientation": result}

