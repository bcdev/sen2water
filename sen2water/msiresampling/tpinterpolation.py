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

from typing import Tuple
import numpy as np
from sen2water.eoutils.eoprocessingifc import BlockAlgorithm


class TpInterpolation(BlockAlgorithm):

    def interpolate_tp(
        self,
        dummy_data: np.ndarray,
        *,
        tp_data: np.array=None,
        resolution: int=60,
        tp_resolution: int=5000,
        image_shape: Tuple[int, int],
        image_chunksize: Tuple[int, int],
        block_id: Tuple[int, int],
    ) -> np.ndarray:
        """
        Block-wise interpolation of tie-point data to the image grid
        using linear interpolation based on pixel coordinates.

        Parameters
        ----------
        dummy_data: np.ndarray
            ignored, a dummy array
        tp_data: np.ndarray
            complete tie-point data array with full extent,
            anchored with its first point at the upper left corner of the target image
        resolution: int
            resolution of target grid in meters
        tp_resolution: int
            resolution of tp_data in meters
        tp_step: Tuple[int, int]
            tuple of tie point step in y and x, e.g. (1, 64)
        image_shape: Tuple[int, int]
            tuple of image shape, e.g. (1830, 1830)
        image_chunksize: Tuple[int, int]
            tuple of image chunksize, e.g. (610, 610)
        block_id: Tuple[int, int]
            tuple of block index

        Returns
        -------
        np.ndarray
            array of interpolated tp_data with the extent of the image block,
            which is usually image_chunksize except for the right and lower border.
        """

        # corner pixel coordinates of the block
        y1 = block_id[0] * image_chunksize[0]
        x1 = block_id[1] * image_chunksize[1]
        y2 = min(y1 + image_chunksize[0], image_shape[0])
        x2 = min(x1 + image_chunksize[1], image_shape[1])

        # # scale and offset
        # if scale_attrs is not None:
        #     tp_data = scale_attrs.scale(tp_data)

        # 1-D pixel coordinates, reference tp pixel coordinates, weights
        # It is the target pixel centre that needs to be in the tp cell.
        y = np.arange(y1, y2).reshape((y2 - y1, 1))
        x = np.arange(x1, x2)
        y_tp = (y * resolution + resolution // 2) // tp_resolution
        x_tp = (x * resolution + resolution // 2) // tp_resolution
        wy = (y * resolution + resolution // 2 - y_tp * tp_resolution) / tp_resolution
        wx = (x * resolution + resolution // 2 - x_tp * tp_resolution) / tp_resolution

        # 2-D interpolation using numpy broadcasting
        result = (
            (1 - wy) * (1 - wx) * tp_data[y_tp, x_tp]
            + (1 - wy) * wx * tp_data[y_tp, x_tp + 1]
            + wy * (1 - wx) * tp_data[y_tp + 1, x_tp]
            + wy * wx * tp_data[y_tp + 1, x_tp + 1]
        ).astype(np.float32)

        return result

    func = interpolate_tp
