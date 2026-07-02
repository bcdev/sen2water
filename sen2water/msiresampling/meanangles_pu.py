# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

import warnings
from typing import Tuple

import numpy as np
import dask.array as da
import xarray as xr
from eopf.computing import MappingDataType, EOProcessingUnit


class MeanAnglesPU(EOProcessingUnit):

    def run(
        self,
        inputs: MappingDataType,
        band_prefix: str = None,
        long_name: str = None,
        is_azimuth_angle=False,
        dims: Tuple[int, int] = None,
        dtype: str = None,
        **kwargs,
    ) -> MappingDataType:
        result_data = da.map_blocks(
            self.mean_angles,
            *[ inputs[band_name][band_name].data for band_name in inputs ],
            is_azimuth_angle=is_azimuth_angle,
            dtype=dtype,
            meta=np.array((), dtype=dtype),
            **kwargs,
        )
        result = xr.DataTree()
        result[band_prefix + "mean"] = xr.DataArray(
            result_data,
            dims=dims,
            attrs={
                "long_name": long_name,
                "units": "degrees",
                "_FillValue": np.nan,
                "coordinates": "crs y x",
            },
        )
        return {"interpolated": result}

    def mean_angles(
        self,
        *angles: np.ndarray,
        is_azimuth_angle=False,
    ) -> np.ndarray:
        """
        Block-wise mean of a stack of angles
        
        Parameters
        ----------
        angles: List[np.ndarray]
            list of blocks of angles, blocks may contain NaN
        is_azimuth_angle: bool
            whether angle may cross 360 degrees

        Returns
        -------
        np.ndarray
            block of angles.
        """
        
        angles_stack = np.stack(angles)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if is_azimuth_angle:
                angles_mask = np.nanmax(angles_stack, axis=0) > 270.0
                angles_stack_mask = np.tile(angles_mask, (angles_stack.shape[0], 1, 1))
                angles_stack[angles_stack_mask] = (angles_stack[angles_stack_mask] + 180.0) % 360.0
            angles_mean = np.nanmean(angles_stack, axis=0)
        if is_azimuth_angle:
            angles_mean[angles_mask] = (angles_mean[angles_mask] + 180.0) % 360.0
        return angles_mean

    func = mean_angles
