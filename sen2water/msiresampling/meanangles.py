# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...
import warnings

import numpy as np
from sen2water.eoutils.eoprocessingifc import BlockAlgorithm


class MeanAngles(BlockAlgorithm):

    def mean_angles(
        self,
        *angles: np.ndarray,
        is_azimuth_angle=False,
    ) -> np.ndarray:
        """
        Block-wise mean of a stack of angles
        
        Parameters
        ----------
        angles: List[np.ndarray]
            list of blocks of angles, blocks may contain NaN
        is_azimuth_angle: bool
            whether angle may cross 360 degrees

        Returns
        -------
        np.ndarray
            block of angles.
        """
        
        angles_stack = np.stack(angles)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if is_azimuth_angle:
                angles_mask = np.nanmax(angles_stack, axis=0) > 270.0
                angles_stack_mask = np.tile(angles_mask, (angles_stack.shape[0], 1, 1))
                angles_stack[angles_stack_mask] = (angles_stack[angles_stack_mask] + 180.0) % 360.0
            angles_mean = np.nanmean(angles_stack, axis=0)
        if is_azimuth_angle:
            angles_mean[angles_mask] = (angles_mean[angles_mask] + 180.0) % 360.0
        return angles_mean

    func = mean_angles
