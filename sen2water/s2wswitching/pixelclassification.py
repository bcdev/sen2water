# -*- coding: utf-8 -*-

"""Combines Idepix flags and additional flags and selected bands into pixel class"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.51"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 0.51:
# add parameters scale_factor and add_offset and apply before comparing against thresholds, fixes misclassification of
import warnings

import xarray as xr
import numpy as np
from sen2water.eoutils.eoprocessingifc import Operator, BlockAlgorithm
from sen2water.eoutils.eoutils import copy_variable
from sen2water.s2wswitching.acolitebands import AcoliteDataset

C2RCC_BANDS = [
    "rhow_B1",
    "rhow_B2",
    "rhow_B3",
    "rhow_B4",
    "rhow_B5",
    "rhow_B6",
    "rhow_B7",
    "rhow_B8A",
]

OCEAN_WAVELENGTHS = [443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190]

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

"""
valids of pixel_classif_flags:
IDEPIX_INVALID = 0;
IDEPIX_CLOUD = 1;
IDEPIX_CLOUD_AMBIGUOUS = 2;
IDEPIX_CLOUD_SURE = 3;
IDEPIX_CLOUD_BUFFER = 4;
IDEPIX_CLOUD_SHADOW = 5;
IDEPIX_SNOW_ICE = 6;
IDEPIX_BRIGHT = 7;
IDEPIX_WHITE = 8;
IDEPIX_COASTLINE = 9;
IDEPIX_LAND = 10;
IDEPIX_CIRRUS_SURE = 11;
IDEPIX_CIRRUS_AMBIGUOUS = 12;
IDEPIX_CLEAR_LAND = 13;
IDEPIX_CLEAR_WATER = 14;
IDEPIX_WATER = 15;
IDEPIX_BRIGHTWHITE = 16;
IDEPIX_VEG_RISK = 17;
IDEPIX_MOUNTAIN_SHADOW = 18;
IDEPIX_POTENTIAL_SHADOW = 19;
IDEPIX_CLUSTERED_CLOUD_SHADOW = 20;

valids of s2wmask:
val binary   zone (colour)
=== ======== =============
192 11000000 land (green)
191 10111111 coastal transition (only near transition water)
... ...      ...
161 10100001 coastal transition (-"-)
160 10100000 coastal land (yellow)
159 10011111 coastal land candidate for coastal inland water (?)
128 10000000 bottom reflection (red)
096 01100000 inland water (lightblue)
095 01011111 transition water
... ...      ...
065 01000001 transition water
064 01000000 ocean (darkblue)
000 00000000 invalid (black)

valids for pixel_class:
0=NO_DATA
1=CLEAR_LAND_OR_VEGETATION
2=CLEAR_OCEAN_WATER
3=CLEAR_INLAND_WATER
4=SNOW_ICE
5=CIRRUS
6=CLOUD_OR_MOUNTAIN_SHADOW
7=AMBIGUOUS_CLOUD
8=CLOUD
9=AC_OUT_OF_BOUNDS
"""


class PixelClassificationAlgorithm(BlockAlgorithm):

    def classify_pixels(
        self,
        pixel_classif_flags: np.ndarray,
        c2rcc_flags: np.ndarray,
        s2wmask: np.ndarray,
        Rrs443_a: np.ndarray,
        Rrs492_a: np.ndarray,
        Rrs560_a: np.ndarray,
        Rrs665_a: np.ndarray,
        Rrs865_a: np.ndarray,
        B2: np.ndarray,
        B4: np.ndarray,
        B11: np.ndarray,
        *quality_masks: np.ndarray,
        scale_factor: np.float32 = None,
        add_offset: np.float32 = None,
    ) -> np.ndarray:
        tecqua = np.any((np.stack(quality_masks) & 15) != 0, axis=0)
        block_shape = pixel_classif_flags.shape

        pixel_class = np.zeros(block_shape, dtype=np.uint8)
        # initialise classes with static s2wmask
        pixel_class[(s2wmask >= 159)] = (
            1  # land, coastal, candidate for ocean if not dynamic land
        )
        pixel_class[(s2wmask == 64)] = 2  # ocean
        pixel_class[((s2wmask > 64) & (s2wmask <= 96))] = 3  # inland water
        pixel_class[(s2wmask == 128)] = (
            6  # bottom reflection, represented as "mountain shadow"
        )

        # update by dynamic land water discrimination
        # coastal with dynamic water gets ocean
        pixel_class[
            ((s2wmask == 160) & (pixel_classif_flags & 0b100000000000000000 == 0))
        ] = 2
        # coastal transition or land with dynamic water gets inland water
        pixel_class[
            ((s2wmask > 160) & (pixel_classif_flags & 0b100000000000000000 == 0))
        ] = 3
        pixel_class[
            ((s2wmask == 159) & (pixel_classif_flags & 0b100000000000000000 == 0))
        ] = 3
        # ocean with dynamic land gets land
        pixel_class[
            ((s2wmask == 64) & (pixel_classif_flags & 0b100000000000000000 != 0))
        ] = 1
        # inland water with dynamic land gets land
        pixel_class[
            (
                (s2wmask <= 96)
                & (s2wmask > 64)
                & (pixel_classif_flags & 0b100000000000000000 != 0)
            )
        ] = 1

        # update with Idepix flags
        # cirrus
        pixel_class[((pixel_classif_flags & 0b1100000000000) != 0)] = 5
        # cloud ambiguous
        pixel_class[((pixel_classif_flags & 0b1100) == 0b0100)] = 7
        # ice
        pixel_class[((pixel_classif_flags & 0b000001000000) != 0)] = 4
        # TODO move to switching between ocean and inland, only apply to ocean (2x)
        vis_a = np.reshape(
            np.dstack((Rrs443_a, Rrs492_a, Rrs560_a, Rrs665_a)), (*Rrs443_a.shape, 4)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            vis_a_mean = np.nanmean(vis_a, axis=-1)
            vis_a_max = np.nanmax(vis_a, axis=-1)
        del vis_a
        # Rrs are in fact rhow here. Thresholds need to be multiplied by PI.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pixel_class[
                (
                    (vis_a_mean > 0.04 * np.pi)
                    & (vis_a_max > 0.05 * np.pi)
                    & (B11 * scale_factor + add_offset < 0.12)
                )
            ] = 4
        pixel_class[
            (
                (Rrs443_a > 0.025 * np.pi)
                & (Rrs865_a > 0.015 * np.pi)
                & (Rrs443_a > Rrs865_a)
                & (B11 * scale_factor + add_offset < 0.12)
            )
        ] = 4
        del vis_a_mean, vis_a_max
        # cloud or mountain shadow
        pixel_class[((pixel_classif_flags & 0b1000000000000100000) != 0)] = 6
        # cloud
        pixel_class[((pixel_classif_flags & 0b11000) != 0)] = 8
        # only apply additional HR-OC cloud test to ocean
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            pixel_class[
                (
                    (pixel_class == 2)
                    & ((c2rcc_flags & 0b1000) != 0)
                    & (
                        ((B2 - B4) / (B2 + B4) > 0.2)
                        | (B11 * scale_factor + add_offset > 0.04)  # TODO times PI?
                        | ((pixel_classif_flags & 0b1000000000000) != 0)
                    )
                )
            ] = 8
        # invalid
        pixel_class[(tecqua != 0)] = 0
        pixel_class[((pixel_classif_flags & 0b1) != 0)] = 0

        return pixel_class

    func = classify_pixels


class PixelClassification(Operator):
    """
    Converts mainly Idepix flags into simple pixel class
    """

    def run(
        self,
        resampled: xr.Dataset,
        idepix: xr.Dataset,
        c2rcc: xr.Dataset,
        acolite: AcoliteDataset,
        s2wmask: xr.Dataset,
    ) -> xr.Dataset:
        """
        Creates pixelclass Dataset with pixel_class

        Parameters
        ----------
        resampled:
            L1B input, resampled to 60m
        idepix:
            Idepix output with pixel_classif_flags
        c2rcc:
            C2RCC output with rhow_B1, .., rhow_B8A, c2rcc_flags
        acolite: str
            ACOLITE output with rhos_443, .., rhos_2202
        s2wmask: str
            Static Sen2Water ocean/inland water/land mask

        Returns
        -------
        xr.Dataset
            Output product with ocean_B1, .., ocean_B8A, status, s2w_status, pixel_classif_flags
        """

        dims = {"y": idepix.sizes["y"], "x": idepix.sizes["x"]}

        acolite_bands = acolite.acolite_bands()
        if "tecqua_mask" in resampled:
            pixel_class_data = PixelClassificationAlgorithm().apply(
                idepix["pixel_classif_flags"].data,
                c2rcc["c2rcc_flags"].data,
                s2wmask["s2wmask"].data[0],
                acolite[acolite_bands[0]].data,
                acolite[acolite_bands[1]].data,
                acolite[acolite_bands[2]].data,
                acolite[acolite_bands[3]].data,
                acolite[acolite_bands[8]].data,
                resampled["B2"].data,
                resampled["B4"].data,
                resampled["B11"].data,
                resampled["tecqua_mask"].data,
                scale_factor=resampled["B2"].attrs["scale_factor"],
                add_offset=resampled["B2"].attrs["add_offset"],
                dtype=np.int8
            )
        else:
            pixel_class_data = PixelClassificationAlgorithm().apply(
                idepix["pixel_classif_flags"].data,
                c2rcc["c2rcc_flags"].data,
                s2wmask["s2wmask"].data[0],
                acolite[acolite_bands[0]].data,
                acolite[acolite_bands[1]].data,
                acolite[acolite_bands[2]].data,
                acolite[acolite_bands[3]].data,
                acolite[acolite_bands[8]].data,
                resampled["B2"].data,
                resampled["B4"].data,
                resampled["B11"].data,
                *[resampled[f"quality_flags_{band}"].data for band in BANDS],
                scale_factor=resampled["B2"].attrs["scale_factor"],
                add_offset=resampled["B2"].attrs["add_offset"],
                dtype=np.int8
            )
        coordinate_bands = {}
        coordinate_bands["y"] = copy_variable(resampled["y"])
        coordinate_bands["x"] = copy_variable(resampled["x"])
        coordinate_bands["crs"] = copy_variable(resampled["crs"])

        target_bands = {}
        target_bands["pixel_class"] = xr.DataArray(
            pixel_class_data,
            attrs={
                "long_name": "Pixel classification flags",
                "flag_meanings": "NO_DATA CLEAR_LAND_OR_VEGETATION CLEAR_OCEAN_WATER "
                + "CLEAR_INLAND_WATER SNOW_ICE CIRRUS CLOUD_OR_MOUNTAIN_SHADOW "
                + "AMBIGUOUS_CLOUD CLOUD AC_OUT_OF_BOUNDS",
                "flag_values": [
                    np.uint8(0),
                    np.uint8(1),
                    np.uint8(2),
                    np.uint8(3),
                    np.uint8(4),
                    np.uint8(5),
                    np.uint8(6),
                    np.uint8(7),
                    np.uint8(8),
                    np.uint8(9),
                ],
                "_FillValue": 0,
            },
            dims=dims,
        )
        target_attrs = {}
        return xr.Dataset(target_bands, coordinate_bands, target_attrs)

