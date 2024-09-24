# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.51"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 0.51:
# change criterion when to rotate angles to avoid discontinuity at ante-meridian

import numpy as np
from sen2water.eoutils.eoprocessingifc import BlockAlgorithm

class AncillaryInterpolation(BlockAlgorithm):

    def interpolate_anc(
        self,
        lat: np.ndarray,
        lon: np.ndarray,
        *,
        anc_lat: np.array=None,
        anc_lon: np.array=None,
        anc_data: np.array=None,
        variable: str=None,
    ) -> np.ndarray:
        """
        Block-wise interpolation of geographically projected data to the image grid
        using linear interpolation based on geo-coordinates.

        Parameters
        ----------
        lat: np.ndarray
            latitudes of target image. 2-d
        lon: np.ndarray
            longitudes of target image, 2-d
        anc_lat: np.ndarray
            latitudes of anc_data, 1-d
        anc_lon: np.ndarray
            longitudes of anc_data, 1-d
        anc_data: np.ndarray
            complete ancillary data array with full extent, 2-d

        Returns
        -------
        np.ndarray
            array of interpolated anc_data with the extent of the lat and lon block
        """

        # shift by 180 degrees if anc may cross antimeridian
        # assumption: anc_lon spans less than 90 degrees
        if anc_lon[0] > 90.0:
            anc_lon = (anc_lon + 360.0) % 360.0 - 180.0
            lon = (lon + 360.0) % 360.0 - 180.0

        lat_start = anc_lat[0]
        lon_start = anc_lon[0]
        lat_end = anc_lat[-1]
        lon_end = anc_lon[-1]
        lat_step = (lat_end - lat_start) / (len(anc_lat) - 1)
        lon_step = (lon_end - lon_start) / (len(anc_lon) - 1)

        # we do not fully trust coverage of the ancillary
        lat_eps = 0.0001 * lat_step
        lon_eps = 0.0001 * lon_step
        lat[lat>lat_start] = lat_start
        lat[lat<=lat_end] = lat_end - lat_eps
        lon[lon<lon_start] = lon_start
        lon[lon>=lon_end] = lon_end - lon_eps

        y = (lat - lat_start) / lat_step
        x = (lon - lon_start) / lon_step
        yi = y.astype(np.int32)
        xi = x.astype(np.int32)
        wy = y - yi
        wx = x - xi

        # 2-D bi-linear interpolation
        result = (
                (1 - wy) * (1 - wx) * anc_data[yi, xi]
                + (1 - wy) * wx * anc_data[yi, xi + 1]
                + wy * (1 - wx) * anc_data[yi + 1, xi]
                + wy * wx * anc_data[yi + 1, xi + 1]
        ).astype(np.float32)

        return result

    func = interpolate_anc
