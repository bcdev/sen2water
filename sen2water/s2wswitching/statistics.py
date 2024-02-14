# -*- coding: utf-8 -*-

"""Counts pixels distinguished by ocean/inland water/land and by clear/snow/cloud"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import dask.array as da


class S2wStatistics(object):

    @staticmethod
    def attributes():
        return [
            "stat_clear_ocean_count",
            "stat_clear_inland_water_count",
            "stat_clear_land_count",
            "stat_snow_ice_ocean_count",
            "stat_snow_ice_inland_water_count",
            "stat_snow_ice_land_count",
            "stat_cloud_ocean_count",
            "stat_cloud_inland_water_count",
            "stat_cloud_land_count",
            "stat_valid_ocean_count",
            "stat_valid_inland_water_count",
            "stat_valid_land_count",
            "stat_valid_count",
        ]

    @staticmethod
    def count_pixels(pixel_class: da.array, s2wmask: da.array):
        return [
            da.count_nonzero((pixel_class == 2)),
            da.count_nonzero((pixel_class == 3)),
            da.count_nonzero((pixel_class == 1)),
            da.count_nonzero(
                ((pixel_class == 4) & ((s2wmask == 64) | (s2wmask == 128)))
            ),
            da.count_nonzero(((pixel_class == 4) & (s2wmask > 64) & (s2wmask <= 96))),
            da.count_nonzero(((pixel_class == 4) & (s2wmask >= 159))),
            da.count_nonzero(
                (
                    (pixel_class >= 5)
                    & (pixel_class <= 8)
                    & ((s2wmask == 64) | (s2wmask == 128))
                )
            ),
            da.count_nonzero(
                (
                    (pixel_class >= 5)
                    & (pixel_class <= 8)
                    & (s2wmask > 64)
                    & (s2wmask <= 96)
                )
            ),
            da.count_nonzero(
                ((pixel_class >= 5) & (pixel_class <= 8) & (s2wmask >= 159))
            ),
            da.count_nonzero(
                (
                    (pixel_class == 2)
                    | ((pixel_class >= 4) & ((s2wmask == 64) | (s2wmask == 128)))
                )
            ),
            da.count_nonzero(
                (
                    (pixel_class == 3)
                    | ((pixel_class >= 4) & (s2wmask > 64) & (s2wmask <= 96))
                )
            ),
            da.count_nonzero(
                ((pixel_class == 1) | ((pixel_class >= 4) & (s2wmask >= 159)))
            ),
            da.count_nonzero((pixel_class > 0)),
        ]
