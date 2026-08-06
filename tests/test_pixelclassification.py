# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Olaf Danne (Brockmann Consult GmbH)"
__copyright__ = "Copyright 2023-2026, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import numpy as np
import pytest

from sen2water.msiidepix.pixelclassification import PixelClassification

dummy = (np.array([[1.0]]), np.array([[2.0]]))
dummy2 = (np.array([3.0]), np.array([2.0]))

refl_clear_land = (np.array([[0.1413]]), np.array([[0.1207]]), np.array([[0.1191]]), np.array([[0.1114]]),
                   np.array([[0.1559]]), np.array([[0.2464]]), np.array([[0.2963]]), np.array([[0.2941]]),
                   np.array([[0.3436]]), np.array([[0.0717]]), np.array([[0.0062]]), np.array([[0.2498]]),
                   np.array([[0.1356]]))

refl_clear_water = (np.array([[0.1562]]), np.array([[0.1304]]), np.array([[0.1174]]), np.array([[0.1032]]),
                   np.array([[0.1002]]), np.array([[0.0752]]), np.array([[0.0784]]), np.array([[0.0618]]),
                   np.array([[0.0612]]), np.array([[0.0145]]), np.array([[0.0021]]), np.array([[0.011]]),
                   np.array([[0.0149]]))

refl_cloud_ambiguous = (np.array([[0.164]]), np.array([[0.1426]]), np.array([[0.1383]]), np.array([[0.1012]]),
                   np.array([[0.1483]]), np.array([[0.3298]]), np.array([[0.4005]]), np.array([[0.3869]]),
                   np.array([[0.4387]]), np.array([[0.095]]), np.array([[0.0235]]), np.array([[0.215]]),
                   np.array([[0.1025]]))

refl_cirrus_ambiguous = (np.array([[0.164]]), np.array([[0.1426]]), np.array([[0.1383]]), np.array([[0.1012]]),
                   np.array([[0.1483]]), np.array([[0.3298]]), np.array([[0.4005]]), np.array([[0.3869]]),
                   np.array([[0.4387]]), np.array([[0.095]]), np.array([[0.0235]]), np.array([[0.215]]),
                   np.array([[0.1025]]))

refl_cloud_sure = (np.array([[0.8183]]), np.array([[0.7886]]), np.array([[0.753]]), np.array([[0.7907]]),
                   np.array([[0.7965]]), np.array([[0.8381]]), np.array([[0.884]]), np.array([[0.8293]]),
                   np.array([[0.8935]]), np.array([[0.3962]]), np.array([[0.0467]]), np.array([[0.3719]]),
                   np.array([[0.265]]))

refl_cirrus_sure = (np.array([[0.8183]]), np.array([[0.7886]]), np.array([[0.753]]), np.array([[0.7907]]),
                   np.array([[0.7965]]), np.array([[0.8381]]), np.array([[0.884]]), np.array([[0.8293]]),
                   np.array([[0.8935]]), np.array([[0.3962]]), np.array([[0.0467]]), np.array([[0.3719]]),
                   np.array([[0.265]]))

refl_bright = (np.array([[0.1853]]), np.array([[0.1552]]), np.array([[0.1491]]), np.array([[0.1311]]),
                   np.array([[0.1786]]), np.array([[0.3053]]), np.array([[0.3741]]), np.array([[0.3645]]),
                   np.array([[0.4304]]), np.array([[0.0867]]), np.array([[0.0054]]), np.array([[0.2773]]),
                   np.array([[0.142]]))

refl_white = (np.array([[0.4539]]), np.array([[0.4276]]), np.array([[0.3962]]), np.array([[0.4012]]),
                   np.array([[0.3964]]), np.array([[0.4153]]), np.array([[0.4372]]), np.array([[0.4101]]),
                   np.array([[0.4431]]), np.array([[0.1784]]), np.array([[0.0207]]), np.array([[0.2427]]),
                   np.array([[0.21]]))

refl_bright_white = (np.array([[0.8183]]), np.array([[0.7886]]), np.array([[0.753]]), np.array([[0.7907]]),
                   np.array([[0.7965]]), np.array([[0.8381]]), np.array([[0.884]]), np.array([[0.8293]]),
                   np.array([[0.8935]]), np.array([[0.3962]]), np.array([[0.0467]]), np.array([[0.3719]]),
                   np.array([[0.265]]))

refl_snow_ice = np.zeros(shape=13)  # TODO


def test_tc4_cirrus_value():
    assert PixelClassification().test_value(dummy) == 3.0
    assert PixelClassification().test_value(dummy2) == 5.0

    tc4_cirrus_value = PixelClassification().tc4_cirrus_value(refl_clear_land)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0538

    tc4_cirrus_value = PixelClassification().tc4_cirrus_value(refl_cloud_sure)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.3353

    tc4_cirrus_value = PixelClassification().tc4_cirrus_value(refl_cloud_ambiguous)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0953

    tc4_cirrus_value = PixelClassification().tc4_cirrus_value(refl_bright)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0715

    tc4_cirrus_value = PixelClassification().tc4_cirrus_value(refl_white)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.198

    print("done test_tc4_cirrus_value")

def test_tc4_value():
    tc4_cirrus_value = PixelClassification().tc4_value(refl_clear_land)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0476

    tc4_cirrus_value = PixelClassification().tc4_value(refl_cloud_sure)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.2886

    tc4_cirrus_value = PixelClassification().tc4_value(refl_cloud_ambiguous)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0718

    tc4_cirrus_value = PixelClassification().tc4_value(refl_bright)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0661

    tc4_cirrus_value = PixelClassification().tc4_value(refl_white)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.1773

    print("done test_tc4_cirrus_value")

def test_ndwi_value():
    ndwi_value = PixelClassification().ndwi_value(refl_clear_land)
    assert pytest.approx(ndwi_value[0][0], 1.E-3) == 0.1581

    ndwi_value = PixelClassification().ndwi_value(refl_cloud_sure)
    assert pytest.approx(ndwi_value[0][0], 1.E-3) == 0.4122

    ndwi_value = PixelClassification().ndwi_value(refl_cloud_ambiguous)
    assert pytest.approx(ndwi_value[0][0], 1.E-3) == 0.3422

    ndwi_value = PixelClassification().ndwi_value(refl_bright)
    assert pytest.approx(ndwi_value[0][0], 1.E-3) == 0.2163

    ndwi_value = PixelClassification().ndwi_value(refl_white)
    assert pytest.approx(ndwi_value[0][0], 1.E-3) == 0.2922

    print("done test_ndwi_value")

def test_b3_b11_value():
    b3_b11_value = PixelClassification().b3_b11_value(refl_clear_land)
    assert pytest.approx(b3_b11_value[0][0], 1.E-3) == 0.4768

    b3_b11_value = PixelClassification().b3_b11_value(refl_cloud_sure)
    assert pytest.approx(b3_b11_value[0][0], 1.E-3) == 2.0247

    b3_b11_value = PixelClassification().b3_b11_value(refl_cloud_ambiguous)
    assert pytest.approx(b3_b11_value[0][0], 1.E-3) == 0.6432

    b3_b11_value = PixelClassification().b3_b11_value(refl_bright)
    assert pytest.approx(b3_b11_value[0][0], 1.E-3) == 0.5377

    b3_b11_value = PixelClassification().b3_b11_value(refl_white)
    assert pytest.approx(b3_b11_value[0][0], 1.E-3) == 1.6324

    print("done test_b3_b11_value")

def test_vis_bright_value():
    vis_bright_value = PixelClassification().vis_bright_value(refl_clear_land)
    assert pytest.approx(vis_bright_value[0][0], 1.E-3) == 0.1170

    vis_bright_value = PixelClassification().vis_bright_value(refl_cloud_sure)
    assert pytest.approx(vis_bright_value[0][0], 1.E-3) == 0.7774

    vis_bright_value = PixelClassification().vis_bright_value(refl_cloud_ambiguous)
    assert pytest.approx(vis_bright_value[0][0], 1.E-3) == 0.1273

    vis_bright_value = PixelClassification().vis_bright_value(refl_bright)
    assert pytest.approx(vis_bright_value[0][0], 1.E-3) == 0.1451

    vis_bright_value = PixelClassification().vis_bright_value(refl_white)
    assert pytest.approx(vis_bright_value[0][0], 1.E-3) == 0.4083

    print("done test_vis_bright_value")

def test_is_cloud_sure():
    assert not PixelClassification().is_cloud_sure(refl_clear_land)[0][0]
    assert PixelClassification().is_cloud_sure(refl_cloud_sure)[0][0]
    assert not PixelClassification().is_cloud_sure(refl_cloud_ambiguous)[0][0]

    print("done test_is_cloud_sure")

def test_is_cloud_ambiguous():
    assert not PixelClassification().is_cloud_ambiguous(refl_clear_land)[0][0]
    assert not PixelClassification().is_cloud_ambiguous(refl_cloud_sure)[0][0]
    assert PixelClassification().is_cloud_ambiguous(refl_cloud_ambiguous)[0][0]

    print("done test_is_cloud_ambiguous")

def test_is_cloud():
    assert not PixelClassification().is_cloud(refl_clear_land)[0][0]
    assert PixelClassification().is_cloud(refl_cloud_sure)[0][0]
    assert PixelClassification().is_cloud(refl_cloud_ambiguous)[0][0]

    print("done test_is_cloud")

def test_is_clear_land():
    assert PixelClassification().is_clear_land(refl_clear_land)[0][0]
    assert not PixelClassification().is_clear_land(refl_clear_water)[0][0]
    assert not PixelClassification().is_clear_land(refl_cloud_sure)[0][0]
    assert not PixelClassification().is_clear_land(refl_cloud_ambiguous)[0][0]

    print("done test_is_clear_land")

def test_is_clear_water():
    assert not PixelClassification().is_clear_water(refl_clear_land)[0][0]
    assert PixelClassification().is_clear_water(refl_clear_water)[0][0]
    assert not PixelClassification().is_clear_water(refl_cloud_sure)[0][0]
    assert not PixelClassification().is_clear_water(refl_cloud_ambiguous)[0][0]

    print("done test_is_clear_water")

def test_is_bright():
    assert PixelClassification().is_bright(refl_bright)[0][0]

    print("done test_is_bright")

def test_is_white():
    assert PixelClassification().is_white(refl_white)[0][0]

    print("done test_is_white")

def test_is_bright_white():
    assert PixelClassification().is_bright_white(refl_bright_white)[0][0]

    print("done test_is_bright_white")

def test_is_snow_ice():
    # TODO
    assert True
    print("done test_is_snow_ice")

def test_is_cirrus():
    # TODO we need elevation
    assert True
    print("done test_is_cirrus")

def test_is_cirrus_ambiguous():
    # TODO we need elevation
    assert True
    print("done test_is_cirrus_ambiguous")

def test_spectral_slope():
    refl_slope_1 = np.array([[50.0]])
    refl_slope_2 = np.array([[100.0]])
    wvl_1 = 450.0
    wvl_2 = 460.0

    slope = PixelClassification().spectral_slope(refl_slope_1, refl_slope_2, wvl_1, wvl_2)
    assert pytest.approx(slope[0][0], 1.E-6) == 5.0

    refl_slope_1 = np.array([[500.0]])
    slope = PixelClassification().spectral_slope(refl_slope_1, refl_slope_2, wvl_1, wvl_2)
    assert pytest.approx(slope[0][0], 1.E-6) == -40.0

    refl_slope_1 = np.array([[50.0]])
    wvl_2 = 450.0
    slope = PixelClassification().spectral_slope(refl_slope_1, refl_slope_2, wvl_1, wvl_2)
    assert np.isinf(slope)

    print("done test_spectral_slope")