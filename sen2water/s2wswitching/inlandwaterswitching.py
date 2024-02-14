# -*- coding: utf-8 -*-

"""Inland water switching of ocean and POLYMER"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import xarray as xr
import numpy as np
from sen2water.eoutils.eoprocessingifc import Operator, BlockAlgorithm
from sen2water.eoutils.eoutils import copy_variable

OCEAN_WAVELENGTHS = [443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190]
POLYMER_BANDS = [
    "Rw443",
    "Rw490",
    "Rw560",
    "Rw665",
    "Rw705",
    "Rw740",
    "Rw783",
    "Rw842",
    "Rw865",
    "Rw945",
    "Rw1375",
    "Rw1610",
    "Rw2190",
]

"""
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
"""


class InlandWaterSwitchingStep(BlockAlgorithm):

    def switch_inland_water(self, *inputs: np.ndarray) -> np.ndarray:
        # disassemble the stack of inputs
        water = inputs[0 : len(OCEAN_WAVELENGTHS)]
        polymer = inputs[len(OCEAN_WAVELENGTHS) : -4]
        polymer_bitmask = inputs[-4].astype(np.uint16)
        pixel_class = inputs[-3].astype(np.uint8)
        sen2water_flags = inputs[-2].astype(np.uint8)
        s2wmask = inputs[-1].astype(np.uint8)

        block_shape = polymer_bitmask.shape

        # insert polymer results into inland water pixels including dynamic water except transition zone
        inland_water_zone = (
            (pixel_class == 3)
            & ((s2wmask == 96) | (s2wmask == np.uint(192)) | (s2wmask == np.uint(159)))
        ).astype(np.bool_)
        for i in range(len(water)):
            water[i][inland_water_zone] = polymer[i][inland_water_zone]
        sen2water_flags[inland_water_zone] |= 32  # with_polymer
        # cancel c2rcc_oor and acolite_negatives where only polymer is used
        sen2water_flags[inland_water_zone] &= np.uint8(255 - 16 - 8 - 2 - 1)
        del inland_water_zone

        # switch in transition areas by distance from ocean
        # We encode distance in s2wmask as 161 .. 191.
        # 160 is distance 0, 192 would be distance 32.
        # Only distances from 1 to 31 will be looked up.
        # We use a sinuid function to weight ocean and polymer.
        # We force distance 0 to scale 0 and distance 32 to scale 1
        d = np.arange(33)
        x = np.tanh(4 * d / 33 - 2) / 2 * np.tanh(2) + 0.5
        scale = x - (32 - d) / 32 * x[0] + d / 32 * (1 - x[32])
        del d, x

        # update water with weighted combination
        distance_from_ocean = np.zeros(block_shape, dtype=np.uint8)
        water_transition_zone = (
            (pixel_class == 3) & (s2wmask > 64) & (s2wmask < 96)
        ).astype(np.bool_)
        coastal_transition_zone = (
            (pixel_class == 3) & (s2wmask > np.uint(160)) & (s2wmask < np.uint(192))
        ).astype(np.bool_)
        distance_from_ocean[water_transition_zone] = s2wmask[
            water_transition_zone
        ] - np.uint(64)
        distance_from_ocean[coastal_transition_zone] = s2wmask[
            coastal_transition_zone
        ] - np.uint(160)
        for b in range(len(water)):
            water[b][water_transition_zone] = (
                scale[distance_from_ocean[water_transition_zone]]
                * polymer[b][water_transition_zone]
                + (1.0 - scale[distance_from_ocean[water_transition_zone]])
                * water[b][water_transition_zone]
            )
            water[b][coastal_transition_zone] = (
                scale[distance_from_ocean[coastal_transition_zone]]
                * polymer[b][coastal_transition_zone]
                + (1.0 - scale[distance_from_ocean[coastal_transition_zone]])
                * water[b][coastal_transition_zone]
            )
        sen2water_flags[
            (water_transition_zone | coastal_transition_zone)
        ] |= 32  # with_polymer
        del water_transition_zone, coastal_transition_zone, distance_from_ocean

        # set polymer bitmask
        sen2water_flags[
            (((sen2water_flags & 32) != 0) & ((polymer_bitmask & 1023) != 0))
        ] |= 4
        # set out of range for all ACs after switching
        pixel_class[((pixel_class > 0) & (sen2water_flags & 7 != 0))] = 9

        outputs = np.concatenate(
            (
                water,
                pixel_class.reshape((1, *block_shape)),
                sen2water_flags.reshape((1, *block_shape)),
            )
        )

        return outputs

    func = switch_inland_water


class InlandWaterSwitching(Operator):
    """
    switches between ocean and POLYMER
    """

    def run(
        self,
        ocean: xr.Dataset,
        polymer: xr.Dataset,
        pixelclass: xr.Dataset,
        s2wmask: xr.Dataset,
    ) -> xr.Dataset:
        """
        Creates water Dataset with Rw*, pixel_class, sen2water_flags

        Parameters
        ----------
        ocean:
            dataset with Rw*, sen2water_flags
        polymer:
            POLYMER output with Rw*, bitmask
        pixelclass
            simple pixel classification with pixel_class
        idepix:
            IDEPIX output with pixel_classif_flags
        s2wmask:
            Sen2Water static mask with ocean, inland water, land, transition zones

        Returns
        -------
        xr.Dataset
            Output product with Rw*, pixel_class, sen2water_mask
        """

        block_shape = polymer["bitmask"].data.shape
        water_bands_stack = InlandWaterSwitchingStep().apply(
            *[ocean[f"Rw{wavelength}"].data for wavelength in OCEAN_WAVELENGTHS],
            *[polymer[b].data for b in POLYMER_BANDS],
            polymer["bitmask"].data,
            pixelclass["pixel_class"].data,
            ocean["sen2water_flags"].data,
            s2wmask["s2wmask"].data[0],
            dtype=np.float32,
            new_axis=0,
            #meta=np.array((len(OCEAN_WAVELENGTHS) + 2, *block_shape), dtype=np.float32),
        )
        water_bands_stack.compute_chunk_sizes()
        dims = ocean.dims

        coordinate_bands = {}
        coordinate_bands["y"] = copy_variable(ocean["y"])
        coordinate_bands["x"] = copy_variable(ocean["x"])

        target_bands = {}
        for i in range(len(OCEAN_WAVELENGTHS)):
            wavelength = OCEAN_WAVELENGTHS[i]
            data = water_bands_stack[i]
            target_bands[f"Rw{wavelength}"] = xr.DataArray(
                data,
                attrs={
                    "long_name": f"Water leaving reflectance at {wavelength} nm",
                    "units": "1",
                    "radiation_wavelength": float(wavelength),
                    "wavelength_unit": "nm",
                    "_FillValue": np.nan,
                },
                dims=dims,
            )
        data = water_bands_stack[-2].astype(np.uint8)
        target_bands["pixel_class"] = xr.DataArray(
            data,
            attrs={
                "long_name": "Pixel classification flags",
                "flag_meanings": "NO_DATA CLEAR_LAND CLEAR_OCEAN_WATER "
                + "CLEAR_INLAND_WATER SNOW_ICE CIRRUS CLOUD_OR_MOUNTAIN_SHADOW "
                + "AMBIGUOUS_CLOUD CLOUD OUT_OF_BOUNDS_SATURATED",
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
        data = water_bands_stack[-1].astype(np.uint8)
        target_bands["sen2water_flags"] = xr.DataArray(
            data,
            attrs={
                "long_name": "quality and algorithm flags",
                "flag_meanings": "c2rcc_oor acolite_negatives polymer_invalid "
                + "with_c2rcc with_acolite with_polymer",
                "flag_masks": [
                    np.uint8(1),
                    np.uint8(2),
                    np.uint8(4),
                    np.uint8(8),
                    np.uint8(16),
                    np.uint8(32),
                ],
            },
            dims=dims,
        )
        target_attrs = {}
        return xr.Dataset(target_bands, coordinate_bands, target_attrs)
