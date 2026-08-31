# noinspection PyPackageRequirements
import numpy as np

from sen2water.eoutils.eoprocessingifc import BlockAlgorithm

# noinspection PyPackageRequirements

class C2rccMsiAlgorithm(BlockAlgorithm):
    """The algorithm to compute C2RCC quantities for S2 MSI.

    @author Olaf Danne, Martin Böttcher
    """

    def do_c2rcc(
            self,
            *toa: np.ndarray,
    ) -> np.ndarray:
        """Computes S2 MSI C2RCC quantities. Implementation logic follows the ESA SNAP Java implementation.

        :param toa: The observed top of atmosphere spectral reflectance
        for all 21 OLCI channels (400 nm to 1020 nm).

        :return: The classif flag.
        """

        result = np.zeros(shape=toa[0].shape, dtype=np.int32)

        # TODO implement

        return result

    func = do_c2rcc
