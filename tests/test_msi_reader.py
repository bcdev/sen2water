# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import os
import xarray as xr
import pytest

@pytest.fixture
def INPUT_DIR():
    return "/windows/tmp/sen2water"

@pytest.fixture
def INPUT(INPUT_DIR):
    return os.path.join(
        INPUT_DIR,
        "S2A_MSIL1C_20230611T103631_N0509_R008_T32UME_20230611T142059.SAFE",
    )

# @pytest.mark.need_files
# @pytest.mark.integration
def test_msi_reader(INPUT):
    ds = xr.open_dataset(INPUT, engine="safe_msi_l1c")
    print(ds)
