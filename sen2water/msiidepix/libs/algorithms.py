#  Copyright (c) 2025 EUMETSAT (contract EUM/CO/19/4600002254/JCh)
#
#  Code authored 2025 by Brockmann Consult GmbH under EUMETSAT Contract
#  EUM/CO/19/4600002254/JCh "Sentinel-3 Synergy Cloud Mask Development".

from importlib import resources
from typing import Any, Tuple
from typing import Literal
from pathlib import Path

# noinspection PyPackageRequirements
import numpy as np
from numpy.typing import NDArray

# noinspection PyPackageRequirements
import scipy.ndimage

# noinspection PyPackageRequirements
from dask import array as da

# noinspection PyPackageRequirements
from dask_image.ndinterp import affine_transform

# noinspection PyPackageRequirements
from scipy.ndimage import generic_filter

from skimage.morphology import binary_dilation
from skimage.draw import line

import numba

from sen2water.msiidepix.interface.algorithm import Algorithm
from sen2water.msiidepix.interface.algorithm import Algorithm
from sen2water.msiidepix.interface.algorithm import BlockAlgorithm
from sen2water.msiidepix.interface.algorithm import OverlapAlgorithm
from sen2water.msiidepix.libs.ann import Ann
from sen2water.msiidepix.libs.datasetadapter import is_dask_array
from sen2water.msiidepix.libs.filters import filter_max

from sen2water.eoutils.eologging import logger

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


class CloudMasks(BlockAlgorithm):
    """The algorithm to compute the cloud masks for OLCI. The computation
    of cloud masks is described in detail in Wevers et al. (2022, hereafter
    referred to as ATBD). This algorithm does not compute cloud border and
    cloud shadow masks.

    Wevers, J., Müller, D., Kirches, G., Quast, R., & Brockmann, C. (2022).
    IdePix for Sentinel-3 OLCI Algorithm Theoretical Basis Document (1.0).
    Zenodo. https://doi.org/10.5281/zenodo.6517333

    @author Ralf Quast
    """

    def __init__(self, dtype: np.dtype = np.uint8, neural_network_path: Path | str | None = None):
        """
        :param neural_network_path: Path to the neural network model file. If None, use default network model.
        """
        super().__init__(dtype)
        package = "idepix_msi.nets"
        name = "class-sequential-i21x11x10x4x3x2o1-0000.net.md"
        if neural_network_path is not None:
            neural_network_path = neural_network_path.resolve()
            if not neural_network_path.exists():
                raise ValueError(f"Could not find the model at {neural_network_path}")
            self._ann = Ann(str(neural_network_path))
            return

        with resources.path(package, name) as resource:
            with open(resource) as r:
                self._ann = Ann(r.name)

    def chunks(self, *inputs: da.Array) -> tuple[int, ...] | None:
        return None

    @property
    def created_axes(self) -> list[int] | None:
        return None

    @property
    def dropped_axes(self) -> list[int]:
        return []

    # noinspection PyMethodMayBeStatic
    def cloud_masks(
        self,
        inv: np.ndarray,
        lnd: np.ndarray,
        gli: np.ndarray,
        nci: np.ndarray,
        psn: np.ndarray,
        eos: np.ndarray,
        *toa: np.ndarray,
        thr_bright_lnd_1: Any = 0.30,
        thr_bright_lnd_2: Any = 0.25,
        thr_bright_sea_1: Any = 0.20,
        thr_bright_sea_2: Any = 0.08,
        thr_bright_swir_1: Any = 0.76,
        thr_bright_swir_2: Any = 0.50,
        thr_clear_snow_min: Any = 0.00,
        thr_clear_snow_max: Any = 1.10,
        thr_cloud_opaque_min: Any = 1.10,
        thr_cloud_opaque_max: Any = 2.75,
        thr_cloud_semitr_min: Any = 2.75,
        thr_cloud_semitr_max: Any = 3.50,
        thr_cloud_mixlnd_min: Any = 3.50,
        thr_cloud_mixlnd_max: Any = 3.85,
        thr_cloud_mixsea_min: Any = 3.50,
        thr_cloud_mixsea_max: Any = 3.75,
    ) -> np.ndarray:
        """Computes cloud masks from (basically) the output of a trained
        artificial neural network fed with top of atmosphere reflectance
        spectra and threshold values supplied as keyword arguments. Some
        custom threshold tests are applied in addition.

        :param inv: A Boolean invalid pixel mask (e.g., a combination of
        OLCI L1B quality flags).
        :param lnd: A Boolean land mask (e.g., a combination of OLCI L1B
        quality flags).
        :param gli: A Boolean sun-glint mask (e.g., the corresponding
        OLCI L1B quality flag).
        :param nci: A Boolean permanent or climatological non sea-ice
        mask. Indicates the absence of sea ice. The mask may be false
        everywhere, in which case it has no effect.
        The mask can be true everywhere, in which case it has no effect.
        :param psn: A Boolean permanent (or climatological) snow mask.
        Indicates the presence of snow-covered land surface. The mask
        may be false everywhere, in which case it has no effect.
        :param eos: A Boolean excess of scattering (EOS) mask [ATBD,
        pp. 10]. The mask is used over permanent (or climatological) snow
        only. The mask may be false everywhere, in which case it has no
        effect.
        :param toa: The observed top of atmosphere spectral reflectance
        for all 21 OLCI channels (400 nm to 1020 nm).
        :param thr_bright_lnd_1: A 443 nm reflectance threshold value
        to mark bright land (exclusive).
        :param thr_bright_lnd_2: A 443 nm reflectance threshold value
        to mark bright land (exclusive).
        :param thr_bright_sea_1: A 865 nm reflectance threshold value
        to mark bright sea (exclusive).
        :param thr_bright_sea_2: A 865 nm reflectance threshold value
        to mark bright sea (exclusive).
        :param thr_bright_swir_1: A 1020 nm reflectance threshold value
        to mark clouds over snow (exclusive). The threshold value is applied
        over permanent snow areas only.
        :param thr_bright_swir_2: A 1020 nm reflectance threshold value
        to mark clouds over snow (exclusive). The threshold value is applied
        over permanent snow areas only.
        :param thr_clear_snow_min: The lower neural network threshold
        value to mark snow and ice (inclusive).
        :param thr_clear_snow_max: The upper neural network threshold
        value to mark snow and ice (exclusive).
        :param thr_cloud_opaque_min: The lower neural network threshold
        value to mark opaque clouds (inclusive).
        :param thr_cloud_opaque_max: The upper neural network threshold
        value to mark opaque clouds (exclusive).
        :param thr_cloud_semitr_min: The lower neural network threshold
        value to mark semi-transparent clouds (inclusive).
        :param thr_cloud_semitr_max: The upper neural network threshold
        value to mark semi-transparent clouds (exclusive).
        :param thr_cloud_mixlnd_min: The lower neural network threshold
        value to mark spatially mixed clouds over land (inclusive).
        :param thr_cloud_mixlnd_max: The upper neural network threshold
        value to mark spatially mixed clouds over land (exclusive).
        :param thr_cloud_mixsea_min: The lower neural network threshold
        value to mark spatially mixed clouds over sea (inclusive).
        :param thr_cloud_mixsea_max: The upper neural network threshold
        value to mark spatially mixed clouds over sea (exclusive).
        :return: The cloud masks.
        """
        # the final cloud masks, initially all zero
        res = np.zeros(lnd.shape, self.dtype)
        # the Boolean sea mask, logical negation of the land mask
        sea = np.logical_not(lnd)
        # the Boolean no-glint mask, logical negation of the glint mask
        ngl = np.logical_not(gli)

        # Boolean bright land mask [ATBD, p. 9]
        bright_lnd_1 = toa[2] > thr_bright_lnd_1
        # Boolean bright land mask [ATBD, p. 9]
        bright_lnd_2 = toa[2] > thr_bright_lnd_2
        # Boolean bright water mask [ATBD, p. 9]
        bright_sea_1 = toa[16] > thr_bright_sea_1
        # Boolean bright water mask [ATBD, p. 9]
        bright_sea_2 = toa[16] > thr_bright_sea_2

        # the neural network values [ATBD, pp. 13]
        nnv = self._neural_network_values(np.stack(toa, axis=2))
        # Boolean invalid mask [ATBD, p. 9]
        inv = self._or(inv, self._not(np.isfinite(nnv)))
        # Boolean clear-sky snow/ice mask [ATBD, p. 14]
        clear_snow = self._and(
                nnv >= thr_clear_snow_min, nnv < thr_clear_snow_max
            )
        # Boolean opaque cloud mask [ATBD, p. 14]
        cloud_opaque = self._and(
            nnv >= thr_cloud_opaque_min, nnv < thr_cloud_opaque_max
        )
        # Boolean semi-transparent cloud mask [ATBD, p. 14]
        cloud_semitr = self._and(
            nnv >= thr_cloud_semitr_min, nnv < thr_cloud_semitr_max
        )
        # Boolean spatially mixed cloud mask over land [ATBD, p. 14]
        cloud_mixlnd = self._and(
            nnv >= thr_cloud_mixlnd_min, nnv < thr_cloud_mixlnd_max
        )
        # Boolean spatially mixed cloud mask over sea [ATBD, p. 14]
        cloud_mixsea = self._and(
            nnv >= thr_cloud_mixsea_min, nnv < thr_cloud_mixsea_max
        )

        # Boolean cloud sure mask over water [ATBD, p. 15]
        cloud_sure_sea = self._and(sea, bright_sea_1, cloud_opaque)
        # Boolean cloud ambiguous mask over water [ATBD, p. 15]
        cloud_ambi_sea = self._and(
            sea,
            self._or(
                self._and(gli, bright_sea_2, cloud_semitr),
                self._and(ngl, bright_sea_2, self._or(cloud_semitr, cloud_mixsea)),
            ),
        )
        # Boolean clear-sky snow/ice mask over water [ATBD, p. 15]
        clear_snow[self._and(sea, nci)] = False

        # Boolean cloud sure mask over land [ATBD, p. 15]
        cloud_sure_lnd = self._and(lnd, bright_lnd_1, cloud_opaque)
        # Boolean cloud ambiguous mask over land [ATBD, p. 15]
        cloud_ambi_lnd = self._and(
            lnd, bright_lnd_2, self._or(cloud_semitr, cloud_mixlnd)
        )

        # Boolean cloud sure mask [ATBD, pp. 15]
        cloud_sure = self._or(cloud_sure_sea, cloud_sure_lnd)
        # Boolean cloud ambiguous mask [ATBD, pp. 15]
        cloud_ambi = self._or(cloud_ambi_sea, cloud_ambi_lnd)

        # Boolean masks modified over permanent snow [ATBD, pp. 15]
        cloud_over_snow = self._or(
            toa[20] > thr_bright_swir_1,
            self._and(toa[20] > thr_bright_swir_2, eos),
        )
        cloud_over_permanent_snow = self._and(psn, cloud_over_snow)
        # restore cloud sure mask, decline cloud ambiguous and clear masks
        cloud_sure[cloud_over_permanent_snow] = True
        cloud_ambi[cloud_over_permanent_snow] = False
        clear_snow[cloud_over_permanent_snow] = False
        clear_over_permanent_snow = self._and(psn, ~cloud_over_snow)
        # restore clear mask, decline cloud masks
        cloud_sure[clear_over_permanent_snow] = False
        cloud_ambi[clear_over_permanent_snow] = False
        clear_snow[clear_over_permanent_snow] = True

        # final cloud, snow/ice and invalid masks [ATBD, pp. 15]
        res[cloud_sure] = FLAG_MASK_CLOUD | FLAG_MASK_CLOUD_SURE
        res[cloud_ambi] = FLAG_MASK_CLOUD | FLAG_MASK_CLOUD_AMBIGUOUS
        res[clear_snow] = FLAG_MASK_CLEAR_SNOW
        res[inv] = FLAG_MASK_UNTRUE
        return res

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

    def _neural_network_values(self, toa: np.ndarray) -> np.ndarray:
        x = toa.reshape((-1, 21))
        y = self._ann.predict(x)
        return y.reshape(toa.shape[:-1])

    compute_block = cloud_masks

    @property
    def name(self) -> str:
        return "cloud_masks"


