"""Generic dask scheduler selection context"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import time
import dask
from sen2water.eoutils.eologging import logger


class Scheduler:

    def __init__(self, mode: str = "threads"):
        self._mode = mode

    def __enter__(self):
        if self._mode == "synchronous":
            dask.config.set(scheduler="synchronous")
        elif self._mode == "threads":
            dask.config.set(scheduler="threads")
        elif self._mode == "processes":
            from dask.distributed import Client, LocalCluster

            self._cluster = LocalCluster(
                processes=True, n_workers=4, threads_per_worker=1
            )
            self._client = Client(self._cluster)
            print(self._client.dashboard_link)
            time.sleep(10)
        else:
            raise ValueError(f"unknown mode {self._mode}")

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._mode == "processes":
            self._client.close()
            self._cluster.close()
