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
from typing import Union
import numpy as np
from sen2water.eoutils.eoprocessingifc import BlockAlgorithm


class Downsampling(BlockAlgorithm):

    def downsample(
        self,
        *inputs: np.ndarray,
        mode: Union[
            "mean", "median", "min", "max", "first",
            "flagand", "flagor", "flagmedianand", "flagmedianor",
            "majority", "detectormean", "masterdetfoo"]='mean',
        factor: int=6,
        is_azimuth_angle: bool=False,
        is_reflectance: bool=False,
    ) -> np.ndarray:
        """Downsampling from e.g. 10m to 60m by mean or other aggregators"""
        data = inputs[0]
        if len(inputs) > 1:
            target_detector_index = inputs[1]
        if len(inputs) > 2:
            detector_index = inputs[2]
        # quick return for mode first
        if mode == 'first':
            j = (factor - 1) // 2
            i = (factor - 1) // 2
            return data[j::factor,i::factor]
        # turn angles by a reference per detector to have them close to +/- 0
        if is_azimuth_angle and len(inputs) > 2:
            reference = self._reference_per_detector(data, detector_index)
            data = self._normalise_angles_180_180(data, reference[detector_index])
        # stack data according to factor
        pixelcontributions = np.stack([data[j::factor,i::factor] for i in range(factor) for j in range(factor)])
        # distinguish modes and aggregate stacked data
        if mode == 'masterdetfoo':
            # for detector index, masks pixels that do not have at least one contrib from master_detfoo
            target_detector_index[np.all(pixelcontributions != target_detector_index, axis=0)] = 0
            result = target_detector_index
        elif mode == 'majority':  # for detector index if not masterdetfoo
            result = np.round(np.mean(pixelcontributions, axis=0)).astype(np.uint8)
        elif mode == 'detectormean' and len(inputs) > 2:
            # for reflectances
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0.0] = np.nan
            self._mask_contributions_of_non_target_detector(target_detector_index, detector_index, factor,
                                                            pixelcontributions)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = np.nanmean(pixelcontributions, axis=0)
            if is_reflectance:
                result[np.isnan(result)] = 0.0
                result = result.astype(np.uint16)
        elif mode == 'mean':
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0] = np.nan
                mean = np.nanmean(pixelcontributions, axis=0)
                mean[np.isnan(mean)] = 0.0
                result = mean.astype(np.uint16)
            else:
                result = np.mean(pixelcontributions, axis=0)
        elif mode == 'median':
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0] = np.nan
                median = np.nanmedian(pixelcontributions, axis=0)
                median[np.isnan(median)] = 0.0
                result = median.astype(np.uint16)
            else:
                result = np.median(pixelcontributions, axis=0)
        elif mode == 'min' or mode == 'flagand':
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0] = np.nan
                min = np.nanmin(pixelcontributions, axis=0)
                min[np.isnan(min)] = 0.0
                result = min.astype(np.uint16)
            else:
                result = np.min(pixelcontributions, axis=0)
        elif mode == 'max' or mode == 'flagor':
            if is_reflectance:
                pixelcontributions = pixelcontributions.astype(np.float32)
                pixelcontributions[pixelcontributions==0] = np.nan
                max = np.nanmax(pixelcontributions, axis=0)
                max[np.isnan(max)] = 0.0
                result = max.astype(np.uint16)
            else:
                result = np.max(pixelcontributions, axis=0)
        elif mode == 'flagmedianand':
            raise ValueError(f"downsampling mode {mode} not yet implemented")
        elif mode == 'flagmedianor':
            raise ValueError(f"downsampling mode {mode} not yet implemented")
        else:
            raise ValueError(f"unknown mode {mode} in downsample")
        # turn back angles per detector
        if is_azimuth_angle and len(inputs) > 2:
            result = self._normalise_angles_0_360(result, reference[target_detector_index])
        return result

    func = downsample

    def _reference_per_detector(self, data:np.ndarray, detector_index:np.ndarray) -> np.ndarray:
        """The first not-NaN angle per detector, 13 values"""
        return np.array([np.nan] +
                        [Downsampling._first(data[((detector_index==detector) & (np.isnan(data)==0))])
                         for detector in range(1, 13)])

    def _mask_contributions_of_non_target_detector(
            self, target_detector_index: np.ndarray, detector_index: np.ndarray, factor: int,
            pixelcontributions: np.ndarray
    ) -> np.ndarray:
        detector_contributions = np.stack([detector_index[j::factor,i::factor] for i in range(factor)
                                           for j in range(factor)])
        pixelcontributions[detector_contributions != target_detector_index] = np.nan

    @staticmethod
    def _normalise_angles_180_180(angles: np.ndarray, reference: np.ndarray) -> np.ndarray:
        return (angles - reference + 540.0) % 360.0 - 180.0

    @staticmethod
    def _normalise_angles_0_360(angles: np.ndarray, reference: np.ndarray) -> np.ndarray:
        return (angles + reference + 360.0) % 360.0

    @staticmethod
    def _first(angles: np.ndarray) -> np.ndarray:
        return angles[0] if len(angles) > 0 else np.nan

