# -*- coding: utf-8 -*-

"""MSI band names and auxiliary information"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

import numpy as np


class IdepixMsiConstants(object):

    B3B11_THRESH = 1.0


    GCW_THRESH = -0.1
    TCW_TC_THRESH = -0.08
    TCW_NDWI_THRESH = 0.4

    GCL_THRESH_DEFAULT = -0.11
    VISBRIGHT_THRESH = 0.12

    FLAG_MASK_UNTRUE: int = np.uint8(1)
    """The flag value of the invalid pixel mask."""
    FLAG_MASK_CLOUD: int = np.uint8(2)
    """The flag value of the cloud mask."""
    FLAG_MASK_CLOUD_AMBIGUOUS: int = np.uint8(4)
    """The flag value of the cloud ambiguous mask."""
    FLAG_MASK_CLOUD_SURE: int = np.uint8(8)
    """The flag value of the cloud sure mask."""
    FLAG_MASK_CLOUD_BORDER: int = np.uint8(16)
    """The flag value of the cloud border mask."""
    FLAG_MASK_CLOUD_SHADOW: int = np.uint8(32)
    """The flag value of the cloud shadow mask."""
    FLAG_MASK_CLEAR_SNOW: int = np.uint8(64)
    """The flag value of the clear-sky snow and ice mask."""
    FLAG_MASK_MOUNTAIN_SHADOW: int = np.uint16(2048)
    """The flag value of the mountain shadow mask"""

    bands = [
        "b01",
        "b02",
        "b03",
        "b04",
        "b05",
        "b06",
        "b07",
        "b08",
        "b8a",
        "b09",
        "b10",
        "b11",
        "b12",
    ]
    resolutions = {
        "b01": 60,
        "b02": 10,
        "b03": 10,
        "b04": 10,
        "b05": 20,
        "b06": 20,
        "b07": 20,
        "b08": 10,
        "b8a": 20,
        "b09": 60,
        "b10": 60,
        "b11": 20,
        "b12": 20,
    }

    flag_band_prefixes = [
        "B_ancillary_lost",
        "B_ancillary_degraded",
        "B_msi_lost",
        "B_msi_degraded",
        "B_defective",
        "B_nodata",
        "B_partially_corrected_crosstalk",
        "B_saturated_l1a",
    ]

    cloud_ice_flags = [
        "B_opaque_clouds",
        "B_cirrus_clouds",
        "B_snow_and_ice_areas"
    ]

    ancillary_bands = [
        "z",
        "aod469",
        "aod550",
        "aod670",
        "aod865",
        "aod1240",
        "bcaod550",
        "duaod550",
        "omaod550",
        "ssaod550",
        "suaod550",
        "tcwv",
        "msl",
        "tco3",
        "u10",
        "v10",
        "r",
    ]

    sun_angles = [
        "sun_zenith",
        "sun_azimuth",
    ]

    detector_colour = [
        [  0,  0,255,255],
        [  0,  0,255,255],
        [  0,255,  0,255],
        [  0,255,  0,255],
        [255,  0,  0,255],
        [255,  0,  0,255],
        [255,175,175,255],
        [255,250,250,255],
        [255,200,  0,255],
        [255,255,  0,255],
        [255,  0,255,255],
        [255,  0,255,255],
    ]
    quality_flags_colour = [
        [255,200,  0,255],
        [255,200,  0,255],
        [255,  0,255,255],
        [255,  0,  0,255],
        [255,  0,  0,255],
        [  0,255,255,255],
        [255,175,175,255],
        [178,  0,  0,255],
    ]
    cloud_ice_colour = [
        [  0,  0,255,255],
        [255,200,  0,255],
        [255,255,255,255],
    ]
