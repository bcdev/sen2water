#  Copyright (c) 2023 EUMETSAT (contract EUM/CO/19/4600002254/JCh)
#
#  Code authored 2023 by Brockmann Consult GmbH under EUMETSAT Contract
#  EUM/CO/19/4600002254/JCh "Sentinel-3 Synergy Cloud Mask Development".

from abc import ABC
from abc import abstractmethod

# noinspection PyPackageRequirements
import numpy as np

# noinspection PyPackageRequirements
import tensorflow as tf

# noinspection PyPackageRequirements
from keras.regularizers import Regularizer

from .layerbuilder import DenseBuilder
from .layerbuilder import InputBuilder
from .layerbuilder import LayerBuilder
from .layerbuilder import UnitRescaling

tfk = tf.keras


class NetBuilder(ABC):
    """Interface class to build neural networks."""

    @abstractmethod
    def build(self, name="model") -> tfk.Model:
        """Builds a new neural network.
        :param name: The name of the model.
        :return: the neural network built.
        """

    @abstractmethod
    def add(self, builder: LayerBuilder) -> object:
        """Adds a new layer to the neural network being build.
        :param builder: To build the new layer.
        :return: itself.
        """

    @abstractmethod
    def add_input(self, n: int) -> object:
        """Adds a new input layer to the neural network being build.
        :param n: The number of input neurons.
        :return: itself.
        """

    @abstractmethod
    def add_dense(
        self,
        n: int,
        activation: str | None = None,
        biases: str | np.ndarray | None = None,
        kernel: str | np.ndarray | None = None,
        kernel_regularizer: str | Regularizer | None = None,
    ) -> object:
        """Adds a new dense layer to the neural network being build.
        :param n: The number of input neurons.
        :param activation: The name of the activation function.
        :param biases: The biases initializer.
        :param kernel: The kernel initializer.
        :param kernel_regularizer: The kernel regularizer.
        :return: itself.
        """

    @abstractmethod
    def add_rescaling_unit(
        self,
        n: int,
        bounds: np.ndarray,
        signed: bool = False,
        invert: bool = False,
    ) -> object:
        """Adds a unit-rescaling layer to the neural network being build.
        :param n: The number of input neurons.
        :param bounds: The bounds. Must be an array of shape `(2, n)`.
        :param signed: Map to `[-1, 1]` rather than `[0, 1]`
        :param invert: Perform the inverse mapping.
        :return: itself.
        """

    @abstractmethod
    def get_layer_count(self) -> int:
        """Returns the number of neural network layers.
        :return: the number of neural network layers.
        """

    @abstractmethod
    def get_neuron_count(self, layer_index: int) -> int:
        """Returns the number of neurons for a layer of interest.

        :param layer_index: The layer index.
        :return: the number of neurons for the layer of interest.
        """


class SnnBuilder(NetBuilder):
    """Builds a new custom Stuttgart Neural Network (SNN) Simulator model.
    Helmut Schiller (deceased) used a streamlined SNN-based model for his
    feed-forward back-propagation (FFBP) neural network software. He
    developed a custom file format for storing his models, too, which is
    different from the original SNN format. We yet use the term SNN to refer
    to FFBP neural networks and file format.
    """

    _builders: list[LayerBuilder]

    def __init__(self):
        self._builders = []

    def build(self, name="sequential") -> tfk.Sequential:
        model = tf.keras.Sequential(name=name)
        for builder in self._builders:
            model.add(builder.build())
        return model

    def add(self, builder: LayerBuilder) -> object:
        self._builders.append(builder)
        return self

    def add_input(self, n: int) -> object:
        return self.add(InputBuilder(n))

    def add_dense(
        self,
        n: int,
        activation: str | None = None,
        biases: str | np.ndarray | None = None,
        kernel: str | np.ndarray | None = None,
        kernel_regularizer: str | Regularizer | None = None,
    ) -> object:
        if activation is not None:
            assert activation == "sigmoid", "Activation is not supported."
        m = self._builders[-1].get_n()
        b = DenseBuilder(n, m)
        b.activation("sigmoid")
        if biases is not None:
            b.biases(biases)
        if kernel is not None:
            if isinstance(kernel, np.ndarray):
                b.kernel(kernel.reshape(n, m).transpose())
            else:
                b.kernel(kernel)
        if kernel_regularizer is not None:
            b.kernel_regularizer(kernel_regularizer)
        b.name(f"{str(len(self._builders)).zfill(2)}_dense")
        return self.add(b)

    def add_rescaling_unit(
        self,
        n: int,
        bounds: np.ndarray,
        signed: bool = False,
        invert: bool = False,
    ) -> object:
        assert signed is False, "Signed unit normalisation is not supported"
        b = UnitRescaling.Builder(n, bounds[[0]], bounds[[1]])
        b.invert(invert)
        b.name(f"{str(len(self._builders)).rjust(2, '0')}_unit_rescaling")
        return self.add(b)

    def get_layer_count(self) -> int:
        return len(self._builders)

    def get_neuron_count(self, layer_index: int) -> int:
        return self._builders[layer_index].get_n()
