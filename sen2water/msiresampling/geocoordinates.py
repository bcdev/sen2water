# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

import numpy as np
from sen2water.eoutils.eoprocessingifc import BlockAlgorithm
from pyproj import Transformer


class GeoCoordinates(BlockAlgorithm):

    def geo_coordinates(
        self,
        x,
        y,
        *,
        transformer: Transformer=None,
    ) -> np.ndarray:
        """
        Transforms UTM coordinates into geographic coordinates

        x: np.ndarray
        """

        lat, lon = transformer.transform(x, y)
        return np.stack((lat, lon))

    func = geo_coordinates
