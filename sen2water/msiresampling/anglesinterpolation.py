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
from typing import Tuple, Any
import numpy as np
from sen2water.eoutils.eoprocessingifc import BlockAlgorithm


class AnglesInterpolation(BlockAlgorithm):

    def interpolate_angles(
        self,
        detector_index: np.ndarray,
        *,
        detectors:np.ndarray=None,
        detector_angles:np.ndarray=None,
        resolution:int=60,
        angles_resolution:int=5000,
        band:str=None,
        is_azimuth_angle=False,
        image_chunksize: Tuple[int, int]=(610,610),
        block_id: Tuple[int, int],
    ):
        """
        Block-wise interpolation of detector-dependent angles provided in a tie-point 5 km grid
        to an image with the extent of the detector index.

        Parameters
        ----------
        detector_index: np.ndarray
            block of detector numbers per pixel in target grid, 1-based, 0 for out-of-swath areas
        detectors: np.ndarray
            array of detector numbers of the layers of detector_angles
        detector_angles: np.ndarray
            3-D array of angles with dims detector_i, row, col
            NaN where there is no acquisition by the respective detector
        resolution: int
            resolution of target grid in meters, 10, 20, or 60
        angles_resolution: int
            resolution of detector_angles in meters, usually 5000
        band: str
            band name, for debugging purposes
        is_azimuth_angle: bool
            whether the angle is an azimuth angle that requires transformation before aggregation
        image_chunksize: Tuple[int, int]
            tuple of image chunksize, e.g. (610, 610)
        block_id: Tuple[int, int]
            tuple of block index

        Returns
        -------
        np.ndarray
            array of interpolated angles with the extent of the detector_index block,
            which is usually image_chunksize except for the right and lower border.
        """

        # # scale and offset
        # if scale_attrs is not None:
        #     tp_data = scale_attrs.scale(tp_data)

        if is_azimuth_angle:
            reference = self._reference_per_detector(detector_angles, detectors)
            detector_angles = AnglesInterpolation._normalise_angles_180_180(detector_angles, reference[detectors].reshape(len(detectors), 1, 1))
        else:
            reference = None

        # corner pixel coordinates of the block
        y1 = block_id[0] * image_chunksize[0]
        x1 = block_id[1] * image_chunksize[1]
        y2 = y1 + detector_index.shape[0]
        x2 = x1 + detector_index.shape[1]

        # 1-D pixel coordinates, reference tp pixel coordinates, weights
        # It is the target pixel centre that needs to be in the tp cell.
        y = np.arange(y1, y2).reshape((y2 - y1, 1))
        x = np.arange(x1, x2)
        y_tp = (y * resolution + resolution // 2) // angles_resolution
        x_tp = (x * resolution + resolution // 2) // angles_resolution
        wy = (y * resolution + resolution // 2 - y_tp * angles_resolution) / angles_resolution
        wx = (x * resolution + resolution // 2 - x_tp * angles_resolution) / angles_resolution

        result = np.empty(detector_index.shape, dtype=np.float32)
        result[:] = np.nan
        for detector_i, detector in enumerate(detectors):
            self._interpolate_angles_of_detector(
                detector_i,
                detector,
                detector_angles,
                detector_index,
                wy,
                wx,
                x_tp,
                y_tp,
                result)

        if is_azimuth_angle:
            result = AnglesInterpolation._normalise_angles_0_360(result, reference[detector_index])

        return result

    func = interpolate_angles

    def _interpolate_angles_of_detector(
            self,
            detector_i: np.ndarray,
            detector: np.ndarray,
            detector_angles: np.ndarray,
            detector_index: np.ndarray,
            wy: np.ndarray,
            wx: np.ndarray,
            x_tp: np.ndarray,
            y_tp: np.ndarray,
            result: np.ndarray
    ) -> np.ndarray:
        angles = detector_angles[detector_i]
        # 2-D interpolation using numpy broadcasting
        result[detector_index == detector] = ((
                (1 - wy) * (1 - wx) * angles[y_tp, x_tp]
                + (1 - wy) * wx * angles[y_tp, x_tp + 1]
                + wy * (1 - wx) * angles[y_tp + 1, x_tp]
                + wy * wx * angles[y_tp + 1, x_tp + 1])[detector_index == detector]).astype(np.float32)

    def _reference_per_detector(self, detector_angles:np.ndarray, detectors:np.ndarray) -> np.ndarray:
        reference = np.empty((13), dtype=np.float32)
        reference[:] = np.nan
        for detector_i, detector in enumerate(detectors):
            reference[detector] = AnglesInterpolation._first(
                detector_angles[detector_i][(np.isnan(detector_angles[detector_i]) == 0)])
        return reference

    def expand_angles_per_detector(self, vz: np.ndarray, va: np.ndarray) -> np.ndarray:
        """
        Fills NaN neighbours of the area with angles per detector, extrapolates angles provided,
        updates and returns vz, va
        """
        d2r = np.pi / 180.0
        r2d = 180.0 / np.pi
        target_shape = vz.shape
        num_detectors, num_rows, num_cols = target_shape

        # x and y are distances at unit altitude of 1 meter
        # x and y have the extent of the angles grid with dims detector_i, row, col
        x = np.empty(target_shape)
        y = np.empty(target_shape)
        x[:] = np.nan
        y[:] = np.nan
        angles_area = (np.isnan(vz) == 0) & (np.isnan(va) == 0)
        x[angles_area] = np.tan(vz[angles_area] * d2r) * np.cos(va[angles_area] * d2r)
        y[angles_area] = np.tan(vz[angles_area] * d2r) * np.sin(va[angles_area] * d2r)

        # x = dxr * row + dxc * col + cx  plane in rows and columns for x part of angles
        # y = dyr * row + dyc * col + cy  ... and for y part of angles

        # we shift row and col independent for y and x to determine dy and dx per detector
        # we patch 0.0 where there is no neighbour
        # dxc etc. have dim detector_i
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            dxc = np.nanmean(x[:, :, 1:] - x[:, :, :-1], axis=(1, 2)).reshape((num_detectors, 1, 1))
            dxr = np.nanmean(x[:, 1:, :] - x[:, :-1, :], axis=(1, 2)).reshape((num_detectors, 1, 1))
            dyc = np.nanmean(y[:, :, 1:] - y[:, :, :-1], axis=(1, 2)).reshape((num_detectors, 1, 1))
            dyr = np.nanmean(y[:, 1:, :] - y[:, :-1, :], axis=(1, 2)).reshape((num_detectors, 1, 1))
        dxc[np.isnan(dxc)] = 0.0
        dxr[np.isnan(dxr)] = 0.0
        dyc[np.isnan(dyc)] = 0.0
        dyr[np.isnan(dyr)] = 0.0

        cols = np.arange(num_cols).reshape((1, 1, num_cols))
        rows = np.arange(num_rows).reshape((1, num_rows, 1))
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            cx = np.nanmean(x - dxc * cols - dxr * rows, axis=(1, 2)).reshape((num_detectors, 1, 1))
            cy = np.nanmean(y - dyc * cols - dyr * rows, axis=(1, 2)).reshape((num_detectors, 1, 1))
        del x, y

        x_fill = dxc * cols + dxr * rows + cx
        y_fill = dyc * cols + dyr * rows + cy

        #fill_area = angles_area==0
        #va[fill_area] = np.arctan2(y_fill, x_fill)[fill_area]
        #vz[fill_area] = np.arctan(np.sqrt(x_fill*x_fill+y_fill*y_fill))[fill_area]
        # fill 8-neighbours of the angles area
        fill_area = np.zeros(target_shape, dtype=np.bool_)
        fill_area[:,:,:-1] |= (angles_area[:,:,:-1] == 0) & angles_area[:,:,1:]
        fill_area[:,:,1:] |= (angles_area[:,:,1:] == 0) & angles_area[:,:,:-1]
        fill_area[:,:-1,:] |= (angles_area[:,:-1,:] == 0) & angles_area[:,1:,:]
        fill_area[:,1:,:] |= (angles_area[:,1:,:] == 0) & angles_area[:,:-1,:]
        fill_area[:,:-1,:-1] |= (angles_area[:,:-1,:-1] == 0) & angles_area[:,1:,1:]
        fill_area[:,1:,1:] |= (angles_area[:,1:,1:] == 0) & angles_area[:,:-1,:-1]
        fill_area[:,:-1,1:] |= (angles_area[:,:-1,1:] == 0) & angles_area[:,1:,:-1]
        fill_area[:,1:,:-1] |= (angles_area[:,1:,:-1] == 0) & angles_area[:,:-1,1:]
        va[fill_area] = np.arctan2(y_fill[fill_area], x_fill[fill_area]) * r2d
        vz[fill_area] = np.arctan(np.sqrt(x_fill[fill_area]*x_fill[fill_area]+y_fill[fill_area]*y_fill[fill_area])) * r2d

        return vz, va

    @staticmethod
    def _first(angles: np.ndarray) -> Any:
        return angles[0] if len(angles) > 0 else np.nan

    @staticmethod
    def _normalise_angles_180_180(angles: np.ndarray, reference: np.ndarray) -> np.ndarray:
        return (angles - reference + 540.0) % 360.0 - 180.0

    @staticmethod
    def _normalise_angles_0_360(angles: np.ndarray, reference: np.ndarray) -> np.ndarray:
        return (angles + reference + 360.0) % 360.0

