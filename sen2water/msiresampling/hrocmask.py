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
import dask.array
import numpy as np
import xarray as xr
import dask.array as da

from sen2water.eoutils.eoprocessingifc import Operator, BlockAlgorithm

BANDS = [
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B8",
    "B8A",
    "B9",
    "B10",
    "B11",
    "B12",
]

class TecquaAlgorithm(BlockAlgorithm):

    def tecquamask(
        self,
        *quality_masks: np.ndarray,
    ) -> np.ndarray:
        return np.any((np.stack(quality_masks) & 15) != 0, axis=0)

    func = tecquamask

class HrocMask(Operator):

    def run(
        self,
        resampled: xr.Dataset,
        hroc_mask: xr.Dataset,
    ) -> xr.Dataset:
        tecqua_data = TecquaAlgorithm().apply(
            *[resampled[f"quality_flags_{band}"].data for band in BANDS],
            dtype=np.int8
        )

        resampled["tecqua_mask"] = xr.DataArray(
            tecqua_data,
            dims=resampled["B1"].dims,
            attrs={
                "long_name": "tecqua mask, a combination of the quality flags of all bands",
            }
        )
        resampled["hroc_watermask"] = xr.DataArray(
            hroc_mask["band_data"].data[0].astype(np.uint8),
            dims=resampled["B1"].dims,
            attrs={
                "long_name": "static HR-OC water mask",
                "flag_values": [np.int32(0), np.int32(1), np.int32(2), np.int32(3), np.int32(4)],
                "flag_meanings": ["invalid", "land", "ocean", "coastal", "bottom_reflection"],
            }
        )
        return resampled

    def run_with_constant(
        self,
        resampled: xr.Dataset,
        value: np.uint8,
    ) -> xr.Dataset:
        tecqua_data = TecquaAlgorithm().apply(
            *[resampled[f"quality_flags_{band}"].data for band in BANDS],
            dtype=np.int8
        )

        resampled["tecqua_mask"] = xr.DataArray(
            tecqua_data,
            dims=resampled["B1"].dims,
            attrs={
                "long_name": "tecqua mask, a combination of the quality flags of all bands",
            }
        )
        ocean_data = da.full_like(tecqua_data, value, chunks=tecqua_data.chunks, dtype=np.uint8)
        resampled["hroc_watermask"] = xr.DataArray(
            ocean_data,
            dims=resampled["B1"].dims,
            attrs={
                "long_name": "static HR-OC water mask",
                "flag_values": [np.int32(0), np.int32(1), np.int32(2), np.int32(3), np.int32(4)],
                "flag_meanings": ["invalid", "land", "ocean", "coastal", "bottom_reflection"],
            }
        )
        return resampled
