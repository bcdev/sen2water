# noinspection PyPackageRequirements
import numpy as np
from xarray import DataTree

from sen2water.eoutils.eoprocessingifc import OverlapAlgorithm


# noinspection PyPackageRequirements


class SlopeAspectOrientation(OverlapAlgorithm):
    """The algorithm to compute the Slope, Aspect and Orientation for S2 MSI.
    Needs elevation and geometries from source.
    We need to consider overlaps because of computation of 3x3 macro pixels

    @author Olaf Danne, Martin Böttcher
    """

    def compute_sao(
            self,
            *geometries: DataTree,
    ) -> type[NotImplementedError]:
        """
        Method to compute the Slope, Aspect and Orientation for S2 MSI.
        Called by SlopeAspectOrientation.().apply(...)

        Parameters
        ----------
        geometries

        Returns
        -------

        """
        # TODO implement
        return NotImplementedError

    func = compute_sao


    # Functions building the pixel_classif_flag...

    def compute_slope_aspect(self, geometries: tuple[np.ndarray, ...]) -> type[NotImplementedError]:
        """

        Parameters
        ----------
        geometries

        Returns
        -------

        """

        # TODO implement
        return NotImplementedError

    def compute_orientation(self, geometries: tuple[np.ndarray, ...]) -> type[NotImplementedError]:
        """

        Parameters
        ----------
        geometries

        Returns
        -------

        """

        # TODO implement
        return NotImplementedError

