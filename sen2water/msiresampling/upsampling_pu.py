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

from typing import Tuple, Union, Literal

import numpy as np
from eopf.computing import EOProcessingUnit, MappingDataType

import xarray as xr
import dask.array as da

from sen2water.eoutils.eoprocessingifc import OverlapAlgorithm


class UpsamplingPU(EOProcessingUnit):

    def run(self,
            inputs: MappingDataType,
            var_name: str,
            band_chunksize: int,
            src_image_shape: Tuple[int, int],
            src_image_chunksize: Tuple[int, int],
            depth: int=0,
            mode: Literal[
                "mean", "median", "min", "max", "first",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean", "masterdetfoo"] = 'mean',
            factor: int = 6,
            is_azimuth_angle: bool = False,
            is_reflectance: bool = False,
            **kwargs
            ) -> MappingDataType:
        input_data = inputs['data'][var_name].data.rechunk(band_chunksize)
        dtype = input_data.dtype
        result_data = da.map_overlap(self.downsample, input_data, dtype=dtype, meta=np.array((), dtype=dtype),
                                     depth=depth,
                                     mode=mode, factor=factor,
                                     src_image_shape=src_image_shape,
                                     src_image_chunksize=src_image_chunksize,
                                     is_azimuth_angle=is_azimuth_angle,
                                     is_reflectance=is_reflectance)
        result = xr.DataTree()
        result[var_name] = xr.DataArray(result_data, dims=inputs['data'][var_name].dims)
        return result

    def upsample(
            self,
            src_data: np.ndarray,
            *,
            mode: Union["nearest", "bilinear", "bicubic"],
            factor: Tuple[int, int],
            src_image_shape: Tuple[int, int],
            src_image_chunksize: Tuple[int, int],
            is_azimuth_angle: bool = False,
            is_reflectance: bool = False,
            block_id: Tuple[int, int],
    ) -> np.ndarray:
        """Upsampling from e.g. 60m to 10m by filling or interpolation"""
        y_factor, x_factor = factor
        if mode == 'nearest':
            # src_data is not extended by map_overlap

            target = np.empty((src_data.shape[0] * y_factor, src_data.shape[1] * x_factor), dtype=src_data.dtype)
            for j in range(y_factor):
                for i in range(x_factor):
                    target[j::y_factor, i::x_factor] = src_data
            return target

        if mode == 'bilinear':
            # map_overlap has extended data by one row and column unless we are at the border.

            # # scale and offset
            # if scale_attrs is not None:
            #     src_data = scale_attrs.scale(src_data)

            # extend tp grid by one column and one row for interpolation unless done by map_overlap
            if block_id[0] == 0:
                # extend to the top
                src_dummy_row = 2 * src_data[0] - src_data[1]
                src_data = np.vstack([src_dummy_row, src_data])
            if (block_id[0] + 1) * src_image_chunksize[0] >= src_image_shape[0]:
                # extend to the bottom
                src_dummy_row = 2 * src_data[-1] - src_data[-2]
                src_data = np.vstack([src_data, src_dummy_row])
            if block_id[1] == 0:
                # extend to the left
                src_dummy_column = 2 * src_data[:, :1] - src_data[:, 1:2]
                src_data = np.hstack([src_dummy_column, src_data])
            if (block_id[1] + 1) * src_image_chunksize[1] >= src_image_shape[1]:
                # extend to the right
                src_dummy_column = 2 * src_data[:, -1:] - src_data[:, -2:-1]
                src_data = np.hstack([src_data, src_dummy_column])

            # The buffered tp band's first pixel is factor/2 left of/above the target border
            # and another 1/2 pixel from the first target pixel centre. src_offset is positive.
            # Example: if factor is 6 then the offset is 3.5 target pixels
            # Example: if factor is 5 then offset is 3 target pixels
            src_y_offset = y_factor / 2 + 0.5
            src_x_offset = x_factor / 2 + 0.5

            # pixel coordinates of the block's corners
            y1 = block_id[0] * src_image_chunksize[0] * y_factor
            x1 = block_id[1] * src_image_chunksize[1] * x_factor
            y2 = min(y1 + src_image_chunksize[0] * y_factor, src_image_shape[0] * y_factor)
            x2 = min(x1 + src_image_chunksize[1] * x_factor, src_image_shape[1] * x_factor)

            # 1-D pixel coordinates, corresponding tp pixel coordinates, weights
            y = np.arange(y2 - y1).reshape((y2 - y1, 1))
            x = np.arange(x2 - x1)
            y_tp = (y + src_y_offset).astype(int) // y_factor
            x_tp = (x + src_x_offset).astype(int) // x_factor
            wy = (y - y_tp * y_factor + src_y_offset) / y_factor
            wx = (x - x_tp * x_factor + src_x_offset) / x_factor

            if is_azimuth_angle:
                reference = src_data[1, 1]
                src_data = (src_data - reference + 540.0) % 360.0 - 180.0

            if is_reflectance:
                src_data = src_data.astype(np.float32)
                src_data[src_data == 0] = np.nan
                result = np.nansum(np.stack([
                    (1 - wy) * (1 - wx) * src_data[y_tp, x_tp],
                    (1 - wy) * wx * src_data[y_tp, x_tp + 1],
                    wy * (1 - wx) * src_data[y_tp + 1, x_tp],
                    + wy * wx * src_data[y_tp + 1, x_tp + 1]
                ]), axis=0).astype(src_data.dtype)

            else:
                # 2-D interpolation using numpy broadcasting
                result = (
                        (1 - wy) * (1 - wx) * src_data[y_tp, x_tp]
                        + (1 - wy) * wx * src_data[y_tp, x_tp + 1]
                        + wy * (1 - wx) * src_data[y_tp + 1, x_tp]
                        + wy * wx * src_data[y_tp + 1, x_tp + 1]
                ).astype(src_data.dtype)

            if is_azimuth_angle:
                result = (result + reference + 540.0) % 360.0 - 180.0
            return result

        if mode == "bicubic":

            # extend tp grid by two columns and two rows for interpolation unless done by map_overlap
            if block_id[0] == 0:
                # extend to the top
                src_dummy_row1 = 3 * src_data[0] - 2 * src_data[1]
                src_dummy_row2 = 2 * src_data[0] - src_data[1]
                src_data = np.vstack([src_dummy_row1, src_dummy_row2, src_data])
            if (block_id[0] + 1) * src_image_chunksize[0] >= src_image_shape[0]:
                # extend to the bottom
                src_dummy_row1 = 3 * src_data[-1] - 2 * src_data[-2]
                src_dummy_row2 = 2 * src_data[-1] - src_data[-2]
                src_data = np.vstack([src_data, src_dummy_row2, src_dummy_row1])
            if block_id[1] == 0:
                # extend to the left
                src_dummy_column1 = 3 * src_data[:, :1] - 2 * src_data[:, 1:2]
                src_dummy_column2 = 2 * src_data[:, :1] - src_data[:, 1:2]
                src_data = np.hstack([src_dummy_column1, src_dummy_column2, src_data])
            if (block_id[1] + 1) * src_image_chunksize[1] >= src_image_shape[1]:
                # extend to the right
                src_dummy_column1 = 3 * src_data[:, -1:] - 2 * src_data[:, -2:-1]
                src_dummy_column2 = 2 * src_data[:, -1:] - src_data[:, -2:-1]
                src_data = np.hstack([src_data, src_dummy_column2, src_dummy_column1])

            y_start = (3 + 1 / y_factor) / 2
            x_start = (3 + 1 / x_factor) / 2
            y_count = (src_data.shape[0] - 4) * y_factor
            x_count = (src_data.shape[1] - 4) * x_factor
            y_end = y_start + (y_count - 1) / y_factor
            x_end = x_start + (x_count - 1) / x_factor
            y_target = np.linspace(y_start, y_end, y_count)
            x_target = np.linspace(x_start, x_end, x_count)
            target_grid = np.stack(np.meshgrid(y_target, x_target, indexing="ij"))
            import scipy
            result = scipy.ndimage.map_coordinates(src_data, target_grid, order=3)
            return result

        raise ValueError(f"unknown upsampling mode {mode}")

    func = upsample
