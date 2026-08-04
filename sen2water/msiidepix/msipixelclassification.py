from typing import Any

# noinspection PyPackageRequirements
import numpy as np
from numpy import dtype, ndarray

from sen2water.eoutils.eoprocessingifc import BlockAlgorithm
from sen2water.msiidepix.constants import IdepixMsiConstants as ic

# noinspection PyPackageRequirements

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

DELTA_RHO_TOA_442_THRESHOLD = 0.03
RHO_TOA_442_THRESHOLD = 0.03
WATER_MASK_SOUTH_BOUND = -58.0

NODATA_VALUE = -1
UNCERTAINTY_VALUE = 0.5
LAND_THRESH = 0.9
WATER_THRESH = 0.9

BRIGHTWHITE_THRESH = 1.5
NDSI_THRESH = 0.6
# BRIGHT_THRESH = 0.25
BRIGHT_THRESH = 0.45  # JW, GK 20220825
BRIGHT_FOR_WHITE_THRESH = 0.8
WHITE_THRESH = 0.9
NDVI_THRESH = 0.5

B3B11_THRESH = 1.0

GCW_THRESH = -0.1
TCW_TC_THRESH = -0.08
TCW_NDWI_THRESH = 0.4
ELEVATION_THRESH = 2000.0


class MsiPixelClassification(BlockAlgorithm):
    """The algorithm to compute the pixel classification for S2 MSI.

    @author Olaf Danne, Martin Böttcher
    """

    def do_classification(
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

        result = np.zeros(shape=toa[0].shape, dtype=np.int32)

        result[self.is_invalid(toa)] = ic.FLAG_MASK_INVALID
        result[self.is_cloud_sure(toa)] = ic.FLAG_MASK_CLOUD_SURE
        result[self.is_cloud_ambiguous(toa)] = ic.FLAG_MASK_CLOUD_AMBIGUOUS
        result[self.is_cloud(toa)] = ic.FLAG_MASK_CLOUD
        result[self.is_cirrus(toa)] = ic.FLAG_MASK_CIRRUS_SURE
        result[self.is_cirrus_ambiguous(toa)] = ic.FLAG_MASK_CIRRUS_AMBIGUOUS
        result[self.is_bright_white(toa)] = ic.FLAG_MASK_BRIGHTWHITE
        result[self.is_clear_snow(toa)] = ic.FLAG_MASK_CLEAR_SNOW
        result[self.is_clear_land(toa)] = ic.FLAG_MASK_CLEAR_LAND
        result[self.is_clear_water(toa)] = ic.FLAG_MASK_CLEAR_WATER

        return result

    func = do_classification


    # Functions building the pixel_classif_flag...

    def is_cloud_sure(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        tc4_cirrus_value = self.tc4_cirrus_value(toa)
        tc4_value = self.tc4_value(toa)
        ndwi_value = self.ndwi_value(toa)

        gcw = tc4_cirrus_value < ic.GCW_THRESH
        tcw = self._and(tc4_value < ic.TCW_TC_THRESH, ndwi_value < ic.TCW_NDWI_THRESH)
        is_b3_b11_water = self.is_b3_b11_water(toa)
        acw = self._and(is_b3_b11_water, (self._or(gcw, tcw)))
        gcl = self._and(self._not(is_b3_b11_water),
                        tc4_cirrus_value < ic.GCL_THRESH_DEFAULT,
                        self.vis_bright_value(toa) > ic.VISBRIGHT_THRESH)

        is_invalid = self.is_invalid(toa)
        is_clear_snow = self.is_clear_snow(toa)

        return self._and(self._not(is_invalid), self._not(is_clear_snow), self._or(acw, gcl))

    def is_cloud_ambiguous(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        tc4_cirrus_value = self.tc4_cirrus_value(toa)
        vis_bright_value = self.vis_bright_value(toa)
        is_b3_b11_water = self.is_b3_b11_water(toa)
        tcl = self._and(self._not(is_b3_b11_water),
                        tc4_cirrus_value < ic.TCL_THRESH,
                        vis_bright_value > ic.VISBRIGHT_THRESH)

        is_invalid = self.is_invalid(toa)
        is_clear_snow = self.is_clear_snow(toa)
        is_cloud_sure = self.is_cloud_sure(toa)

        return self._and(self._not(is_invalid), self._not(is_clear_snow), self._not(is_cloud_sure), tcl)

    def is_cloud(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return self._or(self.is_cloud_sure(toa), self.is_cloud_ambiguous(toa))

    @staticmethod
    def is_cirrus(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        # TODO implement: we need elevation
        return np.zeros(shape=toa[0].shape, dtype=np.uint8)

    @staticmethod
    def is_cirrus_ambiguous(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        # TODO implement: we need elevation
        return np.zeros(shape=toa[0].shape, dtype=np.uint8)

    def is_bright_white(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        is_invalid = self.is_invalid(toa)
        bright_value = self.bright_value(toa)
        white_value = self.white_value(toa)

        return self._and(self._not(is_invalid), bright_value + white_value > BRIGHTWHITE_THRESH)

    @staticmethod
    def is_clear_snow(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        # TODO implement: we need elevation
        return np.zeros(shape=toa[0].shape, dtype=np.uint8)

    def is_clear_land(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        is_invalid = self.is_invalid(toa)
        is_cloud_sure = self.is_cloud_sure(toa)
        is_cloud_ambiguous = self.is_cloud_ambiguous(toa)
        is_cirrus = self.is_cirrus(toa)
        is_cirrus_ambiguous = self.is_cirrus_ambiguous(toa)

        radiometric_land_value = self.radiometric_land_value(toa)
        land_value = radiometric_land_value

        return self._and(self._not(is_invalid),
                         self._not(is_cloud_sure),
                         self._not(is_cloud_ambiguous),
                         self._not(is_cirrus),
                         self._not(is_cirrus_ambiguous),
                         land_value > LAND_THRESH)

    def is_clear_water(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        is_invalid = self.is_invalid(toa)
        is_cloud_sure = self.is_cloud_sure(toa)
        is_cloud_ambiguous = self.is_cloud_ambiguous(toa)
        is_cirrus = self.is_cirrus(toa)
        is_cirrus_ambiguous = self.is_cirrus_ambiguous(toa)
        is_clear_snow = self.is_clear_snow(toa)
        is_bright_white = self.is_bright_white(toa)

        radiometric_water_value = self.radiometric_water_value(toa)
        water_value = radiometric_water_value

        return self._and(self._not(is_invalid),
                         self._not(is_cloud_sure),
                         self._not(is_cloud_ambiguous),
                         self._not(is_cirrus),
                         self._not(is_cirrus_ambiguous),
                         self._not(is_clear_snow),
                         self._not(is_bright_white),
                         water_value > WATER_THRESH)

    @staticmethod
    def is_invalid(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        # TODO implement
        return np.zeros(shape=toa[0].shape, dtype=np.uint8)

    def is_b3_b11_water(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        b3_b11_value = self.b3_b11_value(toa)
        cond = b3_b11_value < ic.B3B11_THRESH
        return  np.where(cond, np.ones(shape=toa[0].shape, dtype=np.int32), np.zeros(shape=toa[0].shape, dtype=np.int32))


    # Functions providing spectral quantities used for classification...

    @staticmethod
    def test_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return toa[0] + toa[1]

    @staticmethod
    def tc4_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return -0.8239 * toa[1] + 0.0849 * toa[2] + 0.4396 * toa[3] - 0.058 * toa[8] + 0.2013 * toa[
            11] - 0.2773 * toa[12]

    def tc4_cirrus_value(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return self.tc4_value(toa) - toa[10]

    @staticmethod
    def tcl_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return 0.3029 * toa[1] + 0.2786 * toa[2] + 0.4733 * toa[3] + 0.5599 * toa[8] + 0.508 * toa[11] + 0.1872 * toa[
            12]

    @staticmethod
    def vis_bright_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return (toa[1] + toa[2] + toa[3]) / 3.0

    @staticmethod
    def b3_b11_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return toa[2] / toa[11]

    @staticmethod
    def ndwi_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return (toa[8] - toa[11]) / (toa[8] + toa[11])

    @staticmethod
    def ndsi_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return (toa[2] - toa[11]) / (toa[2] + toa[11])

    @staticmethod
    def ndvi_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        ndvi = 0.5 * (1.0 + (toa[8] - toa[11]) / (toa[8] + toa[11]))
        return np.max(0.0, np.min(1.0, ndvi))

    @staticmethod
    def radiometric_land_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        cond = toa[7] >= toa[3]
        return np.where(cond, np.ones(shape=toa[0].shape, dtype=np.int32), 0.5 * np.ones(shape=toa[0].shape, dtype=np.int32))

    @staticmethod
    def radiometric_water_value(toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        cond = toa[7] < toa[3]
        return np.where(cond, np.ones(shape=toa[0].shape, dtype=np.int32),
                        0.5 * np.ones(shape=toa[0].shape, dtype=np.int32))

    def a_priori_land_value(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return self.radiometric_land_value(toa)

    def a_priori_water_value(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        return self.radiometric_water_value(toa)

    def spectral_flatness_value(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        slope0 = self.spectral_slope(toa[1], toa[0], list(ic.wvls)[1], list(ic.wvls)[0])
        slope1 = self.spectral_slope(toa[2], toa[3], list(ic.wvls)[2], list(ic.wvls)[3])
        slope2 = self.spectral_slope(toa[4], toa[6], list(ic.wvls)[4], list(ic.wvls)[6])

        flatness = 1.0 - np.abs(1000.0 * (slope0 + slope1 + slope2) / 3.0)
        flatness[np.where(flatness < 0.0)] = 0.0

        return flatness

    def bright_value(self, toa: tuple[np.ndarray, ...]) -> np.ndarray:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        nodata_arr = -1.0 * np.ones(shape=toa[0].shape)
        cond_nodata = self._or(toa[1] <= 0.0, toa[2] <= 0.0, toa[3] <= 0.0, toa[8] <= 0.0, toa[11] <= 0.0, toa[12] <= 0.0)

        bright_value = np.where(cond_nodata, nodata_arr,
                                0.3029 * toa[1] + 0.2786 * toa[2] + 0.4733 * toa[3] + 0.5599 * toa[8] +
                                0.508 * toa[11] + 0.1872 * toa[12])
        return bright_value

    def white_value(self, toa: tuple[np.ndarray, ...]) -> ndarray[tuple[Any, ...], dtype[Any]]:
        """

        Parameters
        ----------
        toa

        Returns
        -------

        """
        nodata_arr = np.zeros(shape=toa[0].shape)
        cond_nodata = self.bright_value(toa) <= BRIGHT_FOR_WHITE_THRESH
        spectral_flatness_value = self.spectral_flatness_value(toa)
        white_value = np.where(cond_nodata, nodata_arr, spectral_flatness_value)
        return white_value


    # Helper functions...

    @staticmethod
    def spectral_slope(toa1: ndarray, toa2: ndarray, wvl1: float, wvl2: float) -> ndarray:
        """

        Parameters
        ----------
        toa1
        toa2
        wvl1
        wvl2

        Returns
        -------

        """
        return (toa2 - toa1) / (wvl2 - wvl1)

    @staticmethod
    def _and(a: np.ndarray, b: np.ndarray, *other: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        a
        b
        other

        Returns
        -------

        """
        r = np.logical_and(a, b)
        for o in other:
            r = np.logical_and(r, o)
        return r

    @staticmethod
    def _not(c: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        c

        Returns
        -------

        """
        return np.logical_not(c)

    @staticmethod
    def _or(a: np.ndarray, b: np.ndarray, *other: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        a
        b
        other

        Returns
        -------

        """
        r = np.logical_or(a, b)
        for o in other:
            r = np.logical_or(r, o)
        return r
