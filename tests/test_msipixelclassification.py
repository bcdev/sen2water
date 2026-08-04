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

from sen2water.msiidepix.msipixelclassification import MsiPixelClassification

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
    assert MsiPixelClassification().test_value(dummy) == 3.0
    assert MsiPixelClassification().test_value(dummy2) == 5.0

    tc4_cirrus_value = MsiPixelClassification().tc4_cirrus_value(refl_clear_land)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0538

    tc4_cirrus_value = MsiPixelClassification().tc4_cirrus_value(refl_cloud_sure)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.3353

    tc4_cirrus_value = MsiPixelClassification().tc4_cirrus_value(refl_cloud_ambiguous)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0953

    tc4_cirrus_value = MsiPixelClassification().tc4_cirrus_value(refl_bright)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0715

    tc4_cirrus_value = MsiPixelClassification().tc4_cirrus_value(refl_white)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.198

    print("done test_tc4_cirrus_value")

def test_tc4_value():
    tc4_cirrus_value = MsiPixelClassification().tc4_value(refl_clear_land)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0476

    tc4_cirrus_value = MsiPixelClassification().tc4_value(refl_cloud_sure)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.2886

    tc4_cirrus_value = MsiPixelClassification().tc4_value(refl_cloud_ambiguous)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0718

    tc4_cirrus_value = MsiPixelClassification().tc4_value(refl_bright)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.0661

    tc4_cirrus_value = MsiPixelClassification().tc4_value(refl_white)
    assert pytest.approx(tc4_cirrus_value[0][0], 1.E-3) == -0.1773

    print("done test_tc4_cirrus_value")
