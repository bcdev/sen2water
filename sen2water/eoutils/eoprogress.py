# -*- coding: utf-8 -*-

"""Generic profiling context"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

class Progress:

    def __init__(self, with_progressbar: bool = False):
        self.with_progressbar = with_progressbar

    def __enter__(self):
        if self.with_progressbar:
            from dask.diagnostics import ProgressBar
            self.progress_bar = ProgressBar(dt=2.0)
            self.progress_bar.register()

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.with_progressbar:
            self.progress_bar.unregister()
