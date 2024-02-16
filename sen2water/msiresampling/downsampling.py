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

from typing import Union
import numpy as np
from sen2water.eoutils.eoprocessingifc import BlockAlgorithm


class Downsampling(BlockAlgorithm):

    def downsample(
        self,
        *inputs: np.ndarray,
        mode: Union[
            "mean", "median", "min", "max",
            "flagand", "flagor", "flagmedianand", "flagmedianor",
            "majority", "detectormean"]='mean',
        factor: int=6,
        is_azimuth_angle: bool=False,
        is_reflectance: bool=False,
    ) -> np.ndarray:
        """Downsampling from e.g. 10m to 60m by mean or other aggregators"""
        data = inputs[0]
        if len(inputs) > 1:
            detector_index = inputs[1]
        if mode == 'first':
            j = (factor - 1) // 2
            i = (factor - 1) // 2
            return data[j::factor,i::factor]
        if is_azimuth_angle:
            reference = self._reference_per_detector(data, detector_index)
            data = self._normalise_angles_180_180(data, reference[detector_index])
        pixelcontributions = np.stack([data[j::factor,i::factor] for i in range(factor) for j in range(factor)])
        if mode == 'detectormean':
            selected_detector = self._select_majority_detector(detector_index, factor, pixelcontributions)
            result = np.nanmean(pixelcontributions, axis=0)
        elif mode == 'majority':  # for detector index
            result = np.round(np.mean(pixelcontributions, axis=0)).astype(np.uint8)
        elif mode == 'mean':
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0] = np.nan
                result = np.nanmean(pixelcontributions, axis=0).astype(np.uint16)
            else:
                result = np.mean(pixelcontributions, axis=0)
        elif mode == 'median':
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0] = np.nan
                result = np.nanmedian(pixelcontributions, axis=0).astype(np.uint16)
            else:
                result = np.median(pixelcontributions, axis=0)
        elif mode == 'min' or mode == 'flagand':
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0] = np.nan
                result = np.nanmin(pixelcontributions, axis=0).astype(np.uint16)
            else:
                result = np.min(pixelcontributions, axis=0)
        elif mode == 'max' or mode == 'flagor':
            result = np.max(pixelcontributions, axis=0)
        elif mode == 'flagmedianand':
            raise ValueError(f"downsampling mode {mode} not yet implemented")
        elif mode == 'flagmedianor':
            raise ValueError(f"downsampling mode {mode} not yet implemented")
        else:
            raise ValueError(f"unknown mode {mode} in downsample")
        if is_azimuth_angle:
            result = self._normalise_angles_0_360(result, reference[selected_detector])
        return result

    func = downsample

    def _reference_per_detector(self, data:np.ndarray, detector_index:np.ndarray) -> np.ndarray:
        """The first not-NaN angle per detector, 13 values"""
        return np.array([np.nan] +
                        [Downsampling._first(data[((detector_index==detector) & (np.isnan(data)==0))])
                         for detector in range(1, 13)])

    def _select_majority_detector(
            self, detector_index: np.ndarray, factor: int, pixelcontributions: np.ndarray) -> np.ndarray:
        """We make use of the fact that only adjacent detectors overlap such that the mean determines majority."""
        detector_contributions = np.stack([detector_index[j::factor,i::factor] for i in range(factor)
                                           for j in range(factor)])
        selected_detector = np.round(np.mean(detector_contributions, axis=0)).astype(np.uint8)
        pixelcontributions[detector_contributions != selected_detector] = np.nan
        return selected_detector

    @staticmethod
    def _normalise_angles_180_180(angles: np.ndarray, reference: np.ndarray) -> np.ndarray:
        return (angles - reference + 540.0) % 360.0 - 180.0

    @staticmethod
    def _normalise_angles_0_360(angles: np.ndarray, reference: np.ndarray) -> np.ndarray:
        return (angles + reference + 360.0) % 360.0

    @staticmethod
    def _first(angles: np.ndarray) -> np.ndarray:
        return angles[0] if len(angles) > 0 else np.nan

