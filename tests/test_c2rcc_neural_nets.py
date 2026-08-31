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


refl_snow_ice = np.zeros(shape=13)  # TODO


def test_nn_iop_rw():
    # TODO
    assert True
    print("done test_nn_iop_rw")

