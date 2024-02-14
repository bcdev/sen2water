# -*- coding: utf-8 -*-

"""Ocean water switching of C2RCC and ACOLITE"""

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
ACOLITE_BANDS = [
    "rhos_443",
    "rhos_492",
    "rhos_560",
    "rhos_665",
    "rhos_704",
    "rhos_740",
    "rhos_783",
    "rhos_833",
    "rhos_865",
    "rhos_1614",
    "rhos_2202",
]
OCEAN_WAVELENGTHS = [443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190]
C2RCC_INDICES = [0, 1, 2, 3, 4, 5, 6, 8]
ACOLITE_INDICES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11]


class CoastalWaterSwitchingAlgorithm(BlockAlgorithm):

    def switch_ocean(self, *inputs: np.ndarray) -> np.ndarray:
        # disassemble the stack of inputs
        c2rcc = inputs[0 : len(C2RCC_BANDS)]
        acolite = inputs[len(C2RCC_BANDS) : -2]
        c2rcc_flags = inputs[-2].astype(np.uint16)
        pixel_class = inputs[-1].astype(np.uint8)

        block_shape = c2rcc[0].shape

        # introduce the names used in HR-OC (we use Rrs for rhow to avoid changing formulas).
        # remember that we have to divide by PI (Rrs=rhow/np.pi) or multiply thresholds by PI.
        Rrs_ref = c2rcc[C2RCC_BANDS.index("rhow_B8A")]
        Rrs_ref[Rrs_ref <= 0] = np.nan
        SGLlow = 10
        SGLhigh = 20
        Rrslimit = 0.0005 * np.pi  # Rrs865_c
        with np.errstate(all="ignore"):
            ratio = (
                c2rcc[C2RCC_BANDS.index("rhow_B3")]
                / c2rcc[C2RCC_BANDS.index("rhow_B8A")]
            )
            ratio2 = (
                c2rcc[C2RCC_BANDS.index("rhow_B4")]
                / c2rcc[C2RCC_BANDS.index("rhow_B3")]
            )
        ratio[ratio <= 0] = np.nan
        f1 = np.log(ratio / SGLlow) / np.log(SGLhigh / SGLlow)
        f2 = np.log(SGLhigh / ratio) / np.log(SGLhigh / SGLlow)

        sen2water_flags = np.zeros(block_shape, dtype=np.uint8)
        sen2water_flags[
            c2rcc_flags & 32775 != 0
        ] |= 1  # c2rcc_oor, Rtosa_OOS=1 Rtosa_OOR=2 Rhow_OOR=4 Rhow_OOS=32768
        del c2rcc_flags
        # cancel over invalid, land, cloud
        # not_water = ((pixel_class == 0) | (pixel_class == 1) | (pixel_class == 8)).astype(np.bool_)
        not_water = ((pixel_class != 2) & (pixel_class != 3)).astype(np.bool_)
        del pixel_class

        outputs = np.empty((len(OCEAN_WAVELENGTHS) + 1, *block_shape))
        dummy_band = np.empty(block_shape)
        dummy_band[:] = np.nan
        for i in range(len(OCEAN_WAVELENGTHS)):
            # skip missing bands B8, B9, B10, B11, B12 for c2rcc, B10 for acolite
            try:
                c2rcc_i = C2RCC_INDICES.index(i)  # 7 for B8A
                Rrs1 = c2rcc[c2rcc_i]
            except:
                Rrs1 = dummy_band
            try:
                acolite_i = ACOLITE_INDICES.index(i)
                Rrs2 = acolite[i]
            except:
                Rrs2 = dummy_band
            # if i < 7 or i == 8:  # restrict masking to B1..B7, B8A
            # if True or i < 7 or i == 8:  # restrict masking to B1..B7, B8A
            if i <= 8:  # restrict masking to B1..B8A
                sen2water_flags[(Rrs2 <= 0)] |= 2  # acolite_negatives
            Rrs1[(Rrs1 <= 0)] = np.nan
            Rrs2[(Rrs2 <= 0)] = np.nan
            Rrs12 = Rrs2 * f2 + Rrs1 * f1
            if i == 8:
                Rrs12_865 = Rrs12
            outputs[i] = np.where(
                not_water,
                np.nan,
                np.where(
                    (ratio <= SGLlow) & (Rrs_ref >= Rrslimit),
                    Rrs2,
                    np.where(
                        (ratio <= SGLlow) & (ratio2 > 1.1),
                        Rrs2,
                        np.where(
                            (ratio > SGLlow)
                            & (ratio < SGLhigh)
                            & (Rrs_ref >= Rrslimit)
                            & np.isfinite(Rrs12),
                            Rrs12,
                            np.where(
                                ratio >= SGLhigh,
                                Rrs1,
                                np.where(
                                    (Rrs_ref < Rrslimit) & (ratio2 < 1.1), Rrs1, np.nan
                                ),
                            ),
                        ),
                    ),
                ),
            )
        del f1, f2, Rrs1, Rrs2

        walgo = np.where(
            not_water,
            np.int8(0),
            np.where(
                (ratio <= SGLlow) & (Rrs_ref >= Rrslimit),
                np.int8(3),
                np.where(
                    (ratio <= SGLlow) & (ratio2 > 1.1),
                    np.int8(4),
                    np.where(
                        (ratio > SGLlow)
                        & (ratio < SGLhigh)
                        & (Rrs_ref >= Rrslimit)
                        & np.isfinite(Rrs12_865),
                        np.int8(2),
                        np.where(
                            ratio >= SGLhigh,
                            np.int8(1),
                            np.where(
                                (Rrs_ref < Rrslimit) & (ratio2 < 1.1),
                                np.int8(1),
                                np.int8(0),
                            ),
                        ),
                    ),
                ),
            ),
        )
        del Rrs12_865

        sen2water_flags[((walgo == 0) | (walgo > 2))] &= (
            255 - 1
        )  # cancel c2rcc_oor where only acolite is used
        sen2water_flags[(walgo < 2)] &= (
            255 - 2
        )  # cancel acolite_negatives where only c2rcc is used
        sen2water_flags[((walgo > 0) & (walgo <= 2))] |= 8  # with_c2rcc
        sen2water_flags[(walgo >= 2)] |= 16  # with_acolite

        outputs[-1] = sen2water_flags

        return outputs

    func = switch_ocean


class CoastalWaterSwitching(Operator):
    """
    Switches between C2RCC and ACOLITE
    """

    def run(
        self, c2rcc: xr.Dataset, acolite: xr.Dataset, pixelclass: xr.Dataset
    ) -> xr.Dataset:
        """
        Creates ocean Dataset with Rw* and sen2water_flags

        Parameters
        ----------
        c2rcc:
            C2RCC output with rhow_B1, .., rhow_B8A, c2rcc_flags
        acolite:
            ACOLITE output with rhos_443, .., rhos_2202
        pixelclass:
            simple pixel classification with pixel_class

        Returns
        -------
        xr.Dataset
            Output product with ocean_B1, .., ocean_B8A, sen2water_flags
        """

        block_shape = c2rcc["c2rcc_flags"].data.shape
        ocean_bands_stack = CoastalWaterSwitchingAlgorithm().apply(
            *[c2rcc[b].data for b in C2RCC_BANDS],
            *[acolite[b].data for b in ACOLITE_BANDS],
            c2rcc["c2rcc_flags"].data,
            pixelclass["pixel_class"].data,
            dtype=np.float32,
            new_axis=0,
            #meta=np.array((len(OCEAN_WAVELENGTHS) + 2, *block_shape), dtype=np.float32),
        )
        dims = {"y": c2rcc.sizes["y"], "x": c2rcc.sizes["x"]}
        ocean_bands_stack.compute_chunk_sizes()

        coordinate_bands = {}
        coordinate_bands["y"] = copy_variable(c2rcc["y"])
        coordinate_bands["x"] = copy_variable(c2rcc["x"])
        #coordinate_bands["crs"] = copy_variable(c2rcc["crs"])

        target_bands = {}
        for i in range(len(OCEAN_WAVELENGTHS)):
            wavelength = OCEAN_WAVELENGTHS[i]
            data = ocean_bands_stack[i]
            target_bands[f"Rw{wavelength}"] = xr.DataArray(
                data,
                attrs={
                    "long_name": f"Water leaving reflectance at {wavelength} nm",
                    "units": "1",
                    "radiation_wavelength": float(wavelength),
                    "radiation_wavelength_unit": "nm",
                    "_FillValue": np.nan,
                },
                dims=dims,
            )
        data = ocean_bands_stack[-1].astype(np.uint8)
        target_bands["sen2water_flags"] = xr.DataArray(
            data,
            attrs={
                "long_name": "quality and algorithm flags",
                "flag_meanings": "c2rcc_oor acolite_negatives polymer_invalid with_c2rcc with_acolite with_polymer",
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
