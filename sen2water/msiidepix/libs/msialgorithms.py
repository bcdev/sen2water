from pathlib import Path
from typing import Any

# noinspection PyPackageRequirements
import numpy as np
# noinspection PyPackageRequirements
from dask import array as da

from sen2water.msiidepix.constants import IdepixMsiConstants as ic
from sen2water.msiidepix.interface.algorithm import BlockAlgorithm

# noinspection PyPackageRequirements
# import numba

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

_F = np.uint8(0)
"""Numeric representation of Boolean value 'false'. """

_T = np.uint8(1)
"""Numeric representation of Boolean value 'true'. """

CLOUD_HEIGHT_MAX = 12_000
"""Maximum cloud height to consider for cloud shadow"""

MEAN_EARTH_RADIUS_METERS = 6372000
"""Average radius of the earth (m)"""

DELTA_RHO_TOA_442_THRESHOLD = 0.03;
RHO_TOA_442_THRESHOLD = 0.03;
WATER_MASK_SOUTH_BOUND = -58.0;

UNCERTAINTY_VALUE = 0.5;
LAND_THRESH = 0.9;
WATER_THRESH = 0.9;

BRIGHTWHITE_THRESH = 1.5;
NDSI_THRESH = 0.6;
BRIGHT_THRESH = 0.25;
BRIGHT_THRESH = 0.45;      # JW, GK 20220825
BRIGHT_FOR_WHITE_THRESH = 0.8;
WHITE_THRESH = 0.9;
NDVI_THRESH = 0.5;

B3B11_THRESH = 1.0;

GCW_THRESH = -0.1;
TCW_TC_THRESH = -0.08;
TCW_NDWI_THRESH = 0.4;
ELEVATION_THRESH = 2000.0;


class MsiPixelClassification(BlockAlgorithm):
    """The algorithm to compute the pixel classification for S2 MSI.

    @author Olaf Danne, Martin Böttcher
    """

    def __init__(self, dtype: np.dtype = np.uint8, neural_network_path: Path | str | None = None):
        """
        """
        super().__init__(dtype)
        # TODO add more if needed

    def chunks(self, *inputs: da.Array) -> tuple[int, ...] | None:
        return None

    @property
    def created_axes(self) -> list[int] | None:
        return None

    @property
    def dropped_axes(self) -> list[int]:
        return []

    # noinspection PyMethodMayBeStatic
    def do_classif(
        self,
        *toa: np.ndarray,
        thresh_cw: Any = 0.007,
        thresh_gcl: Any = -0.11,
        thresh_cl: Any = 0.007,
    ) -> np.ndarray:
        """Computes S2 MSI pixel classif flags. Implementation logig follows the ESA SNAP Java implementation.

        :param toa: The observed top of atmosphere spectral reflectance
        for all 21 OLCI channels (400 nm to 1020 nm).

        :param thresh_cw: A 443 nm reflectance threshold value
        to mark bright land (exclusive).
        :param thresh_gcl: A 443 nm reflectance threshold value
        to mark bright land (exclusive).
        :param thresh_cl: A 865 nm reflectance threshold value
        to mark bright sea (exclusive).

        :return: The classif flag.
        """
        # classification logic
        tc4_cirrus_value = -0.8239 * toa[1] + 0.0849 * toa[2] + 0.4396 * toa[3] - 0.058 * toa[8] + 0.2013 * toa[
            11] - 0.2773 * toa[12] - toa[10]

        tc4_value = -0.8239 * toa[1] + 0.0849 * toa[2] + 0.4396 * toa[3] - 0.058 * toa[8] + 0.2013 * toa[
            11] - 0.2773 * toa[12]

        ndwi_value = (toa[8] - toa[11]) / (toa[8] + toa[11])

        b3_b11_value = (toa[2] / toa[11])

        vis_bright_value = (toa[1] + toa[2] + toa[3]) / 3.0

        is_b3_b11 = b3_b11_value > ic.B3B11_THRESH

        is_invalid = np.zeros(toa[0].shape, dtype=bool)  # TODO implement
        is_clear_snow = np.zeros(toa[0].shape, dtype=bool)  # TODO implement

        gcw = tc4_cirrus_value < ic.GCW_THRESH
        tcw = self._and(tc4_value < ic.TCW_TC_THRESH, ndwi_value < ic.TCW_NDWI_THRESH)
        acw = self._and(is_b3_b11, (self._or(gcw, tcw)))
        gcl = self._and(self._not(is_b3_b11), tc4_cirrus_value < ic.GCL_THRESH_DEFAULT,
                        vis_bright_value > ic.VISBRIGHT_THRESH)
        is_cloud_sure = self._and(self._not(is_invalid), self._not(is_clear_snow), self._or(acw, gcl))

        is_cloud_ambiguous = b3_b11_value > ic.B3B11_THRESH

        is_cloud = self._or(is_cloud_sure, is_cloud_ambiguous)

        # final classification
        result = np.zeros(shape=toa[0].shape, dtype=np.int32)

        result[is_cloud] = ic.FLAG_MASK_CLOUD
        result[is_cloud_sure] = ic.FLAG_MASK_CLOUD | ic.FLAG_MASK_CLOUD_SURE
        result[is_cloud_ambiguous] = ic.FLAG_MASK_CLOUD | ic.FLAG_MASK_CLOUD_AMBIGUOUS

        return result


    func = do_classif

    @staticmethod
    def _and(a: np.ndarray, b: np.ndarray, *other: np.ndarray) -> np.ndarray:
        r = np.logical_and(a, b)
        for o in other:
            r = np.logical_and(r, o)
        return r

    @staticmethod
    def _not(c: np.ndarray) -> np.ndarray:
        return np.logical_not(c)

    @staticmethod
    def _or(a: np.ndarray, b: np.ndarray, *other: np.ndarray) -> np.ndarray:
        r = np.logical_or(a, b)
        for o in other:
            r = np.logical_or(r, o)
        return r

    @property
    def name(self) -> str:
        return "cloud_masks"


