#  Copyright (c) 2023 EUMETSAT (contract EUM/CO/19/4600002254/JCh)
#
#  Code authored 2023 by Brockmann Consult GmbH under EUMETSAT Contract
#  EUM/CO/19/4600002254/JCh "Sentinel-3 Synergy Cloud Mask Development".

from pathlib import Path

# noinspection PyPackageRequirements
import numpy as np

from .netreader import SnnReader


class Ann:
    """The class to encapsulate the TensorFlow neural network model."""

    _model_path: str | Path
    """The path to the neural network model file. """

    def __init__(self, model_path: str | Path):
        """Creates a new neural network instance.

        :param model_path: The path to the neural network model file.
        """
        self._model_path = model_path
        self._model = Ann._load_model(model_path)

    def __reduce__(self):
        """For serialization. Does not serialize the neural network model
        by intention."""
        return self.__class__, (self._model_path,)

    def predict(self, features: np.ndarray) -> np.ndarray:
        """Returns the predictions for given feature vectors. A feature
        vector comprises the top-of-atmosphere spectral reflectance for all
        21 OLCI channels.

        :param features: The feature vectors.
        :return: The predictions.
        """
        assert features.ndim == 2, f"Expected a two-dimensional array of features"
        assert (
            features.shape[1] == 21
        ), f"The second dimension is not the spectral dimension"
        return self._model(np.sqrt(features)).numpy().squeeze()

    @staticmethod
    def _load_model(path: str | Path):
        """Loads a neural network model from a file.

        :param path: The path to the neural network model (SNN format).
        :return: The neural network model.
        """
        model = SnnReader().read(path)
        model.compile()
        return model
