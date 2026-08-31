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


class C2rccMsiConstants(object):

    B3B11_THRESH = 1.0


    GCW_THRESH = -0.1
    TCW_TC_THRESH = -0.08
    TCW_NDWI_THRESH = 0.4

    GCL_THRESH_DEFAULT = -0.11
    VISBRIGHT_THRESH = 0.12
    TCL_THRESH = 0.085

    FLAG_MASK_INVALID: int = np.int32(1)
    """The flag value of the invalid pixel mask."""
    FLAG_MASK_CLOUD: int = np.int32(2)
    """The flag value of the cloud mask."""
    FLAG_MASK_CLOUD_AMBIGUOUS: int = np.int32(4)
    """The flag value of the cloud ambiguous mask."""
    FLAG_MASK_CLOUD_SURE: int = np.int32(8)
    """The flag value of the cloud sure mask."""
    FLAG_MASK_CLOUD_BORDER: int = np.int32(16)
    """The flag value of the cloud border mask."""
    FLAG_MASK_CLOUD_SHADOW: int = np.int32(32)
    """The flag value of the cloud shadow mask."""
    FLAG_MASK_CLEAR_SNOW: int = np.int32(64)
    """The flag value of the clear-sky snow and ice mask."""
    FLAG_MASK_BRIGHT: int = np.int32(128)
    """The flag value of the bright mask."""
    FLAG_MASK_WHITE: int = np.int32(256)
    """The flag value of the white mask."""
    FLAG_MASK_COASTLINE: int = np.int32(512)
    """The flag value of the coastline mask."""
    FLAG_MASK_LAND: int = np.int32(1024)
    """The flag value of the land mask."""
    FLAG_MASK_CIRRUS_SURE: int = np.int32(2048)
    """The flag value of the cirrus sure mask."""
    FLAG_MASK_CIRRUS_AMBIGUOUS: int = np.int32(4096)
    """The flag value of the cirrus ambiguous mask."""
    FLAG_MASK_CLEAR_LAND: int = np.int32(8192)
    """The flag value of the clear land mask."""
    FLAG_MASK_CLEAR_WATER: int = np.int32(16384)
    """The flag value of the clear water mask."""
    FLAG_MASK_WATER: int = np.int32(32768)
    """The flag value of the water mask."""
    FLAG_MASK_BRIGHTWHITE: int = np.int32(65536)
    """The flag value of the brightwhite mask."""
    FLAG_MASK_VEG_RISK: int = np.int32(131072)
    """The flag value of the veg risk mask."""
    FLAG_MASK_MOUNTAIN_SHADOW: int = np.int32(262144)
    """The flag value of the mountain shadow mask"""
    FLAG_MASK_POTENTIAL_SHADOW: int = np.int32(524288)
    """The flag value of the potential shadow mask"""
    FLAG_MASK_CLUSTERED_CLOUD_SHADOW: int = np.int32(524288)
    """The flag value of the clustered cloud shadow mask"""

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

    wvls = {
        443.0,     # B1
        490.0,     # B2
        560.0,     # B3
        665.0,     # B4
        705.0,     # B5
        740.0,     # B6
        783.0,     # B7
        842.0,     # B8
        865.0,     # B8A
        945.0,     # B9
        1375.0,    # B10
        1610.0,    # B11
        2190.0     # B12
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
