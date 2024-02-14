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

import numpy as np

from sen2water.msiresampling.anglesinterpolation import AnglesInterpolation
from sen2water.msiresampling.downsampling import Downsampling
from sen2water.msiresampling.upsampling import Upsampling
from sen2water.msiresampling.tpinterpolation import TpInterpolation
from tests.msi_testdata import MsiTestData
import pytest

def test_upsampler():
    data = np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 8.0]])
    mode = 'bilinear'
    factor = 3
    image_shape = (2, 3)
    image_chunksize = (2, 3)
    block_id = (0, 0)
    result = Upsampling().upsample(
        data,
        mode=mode,
        factor=(factor, factor),
        src_image_shape=image_shape,
        src_image_chunksize=image_chunksize,
        block_id=block_id)
    print()
    print(result)
    if factor == 3:
        assert result[1,1] == data[0,0]
        assert result[4,4] == data[1,1]
        assert result[4,7] == data[1,2]
    print("done")

def test_downsampler():
    data = np.arange(36).reshape((6,6))
    mode = "mean"
    factor = 2
    result = Downsampling().downsample(data, mode=mode, factor=factor)
    print()
    print(result)
    assert result[1,2] == (16+17+22+23)/4

def test_interpolation():
    tp_data = np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 8.0], [3.0, 6.0, 12.0]])
    data = np.arange(166*166).reshape((166,166))
    result = TpInterpolation().interpolate_tp(
        data,
        tp_data=tp_data,
        resolution=60,
        tp_resolution=5000,
        image_shape=(166,166),
        image_chunksize=(166,166),
        block_id=(0,0))
    print()
    print(result)
    print("done")

def test_angles_interpolation():
    detector_index = MsiTestData.detector_index
    vz = MsiTestData.vza
    va = MsiTestData.vaa
    detectors = np.array([i for i in vz])
    vz5 = np.stack([vz[i] for i in vz])[:,0:5,0:5]
    va5 = np.stack([va[i] for i in va])[:,0:5,0:5]
    extended_vza, extended_vaa = AnglesInterpolation().expand_angles_per_detector(vz5, va5)

    result = AnglesInterpolation().interpolate_angles(
        detector_index,
        detectors=detectors,
        detector_angles=extended_vza,
        resolution=600,
        angles_resolution=5000,
        image_chunksize=(24,24),
        block_id=[0,0]
    )
    print()
    print(result)
    assert abs(result[1,1] - 3.2021043) < 0.001
    result = AnglesInterpolation().interpolate_angles(
        detector_index,
        detectors=detectors,
        detector_angles=extended_vaa,
        resolution=600,
        angles_resolution=5000,
        is_azimuth_angle=True,
        image_chunksize=(24,24),
        block_id=[0,0]
    )
    print()
    print(result)
    assert abs(result[1,1] - 142.8894) < 0.001