# -*- coding: utf-8 -*-

"""MSI band names and auxiliary information"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

class MsiConstants(object):

    bands = [
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
    resolutions = {
        "B1": 60,
        "B2": 10,
        "B3": 10,
        "B4": 10,
        "B5": 20,
        "B6": 20,
        "B7": 20,
        "B8": 10,
        "B8A": 20,
        "B9": 60,
        "B10": 60,
        "B11": 20,
        "B12": 20,
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
