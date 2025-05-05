# -*- coding: utf-8 -*-

"""Generic logger"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

import logging

logging._levelToName = {
    logging.CRITICAL: "CRIT",
    logging.ERROR: "ERRO",
    logging.WARNING: "WARN",
    logging.INFO: "INFO",
    logging.DEBUG: "DEBU",
    logging.NOTSET: "NONE",
}

formatter = logging.Formatter(
    fmt="%(asctime)s.%(msecs)03d %(levelname)4s %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
handler = logging.StreamHandler()
handler.setFormatter(formatter)
handler.setLevel(logging.INFO)
logger = logging.getLogger("eo")
if not logger.hasHandlers():
    logger.addHandler(handler)
logger.setLevel(logging.INFO)
