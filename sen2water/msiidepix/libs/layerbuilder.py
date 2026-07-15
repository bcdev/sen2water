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
from keras.layers import Rescaling

# noinspection PyPackageRequirements
from keras.regularizers import Regularizer

tfk = tf.keras


class LayerBuilder(ABC):
    """Interface class to build neural network layers."""

    @abstractmethod
    def get_n(self) -> int:
        """Returns the number of neurons within the layer to build."""

    @abstractmethod
    def build(self) -> tfk.layers.Layer:
        """Returns a new layer."""


class InputBuilder(LayerBuilder):
    """Class to build an input layer."""

    _n: int

    def __init__(self, n: int):
        """Creates a new instance of this class.

        :param n: The number of neurons in the input layer.
        """
        self._n = n

    def get_n(self) -> int:
        return self._n

    def build(self) -> tfk.layers.Layer:
        return tfk.layers.Input((self._n,), name="inputs")


class DenseBuilder(LayerBuilder):
    """Class to build a dense layer. The layer computes its output from its
    input by calling `output = f(matmul(input, kernel) + biases` where `f`
    is an activation function and `biases` and `kernel` are (trainable)
    parameters.

    Please note that the custom SNN model computes its output from its input
    by calling `output = activation(matmul( kernel.transpose(), input) +
    biases`. When reading a custom SSN model, the kernel read must be
    transposed when used for a dense layer. When writing a dense layer to an
    SSN model, the kernel must be transposed back again.
    """

    _n: int
    _m: int
    _activation: str | None
    _biases: str | tfk.initializers.Initializer | None
    _kernel: str | tfk.initializers.Initializer | None
    _kernel_regularizer: str | Regularizer | None
    _name: str | None

    def __init__(self, n: int, m: int):
        """Creates a new instance of this class.

        :param n: The number of neurons in the layer.
        :param m: The number of neurons in the preceding layer.
        """
        self._n = n
        self._m = m
        self._activation = None
        self._biases = None
        self._kernel = None
        self._kernel_regularizer = None
        self._name = None

    def activation(self, activation: str) -> object:
        """Defines the activation function.

        :param activation: The name of the activation function.
        :return: itself.
        """
        self._activation = activation
        return self

    def biases(self, b: str | np.ndarray) -> object:
        """To initialize the biases.

        The biases are an array of shape `(n,)` where `n` is the number of
        neurons in the layer being build.

        :param b: The name of an initializer class or the concrete biases.
        :return: itself.
        """
        if isinstance(b, np.ndarray):
            assert b.size == self._n
            self._biases = WeightInitializer(b)
        else:
            self._biases = b
        return self

    def kernel(self, k: str | np.ndarray) -> object:
        """To initialize the kernel.

        The kernel is an array of shape `(m, n)` where `n` is the number of
        neurons in the layer being built and `m` is the number of neurons in
        the preceding layer. If necessary, the array supplied as argument is
        reshaped to initialize the kernel. Array elements are not reordered
        unless transposition is requested explicitly.

        :param k: The name of an initializer class or the concrete kernel.
        :return: itself.
        """
        if isinstance(k, np.ndarray):
            assert k.size == self._n * self._m
            self._kernel = WeightInitializer(k)
        else:
            self._kernel = k
        return self

    def kernel_regularizer(self, r: str | Regularizer | None) -> object:
        self._kernel_regularizer = r
        return self

    def name(self, name: str) -> object:
        """Defines the name of the layer.

        :param name: The name of the layer.
        :return: itself.
        """
        self._name = name
        return self

    def get_n(self) -> int:
        return self._n

    def get_m(self) -> int:
        """Returns the number of neurons in the preceding layer."""
        return self._m

    def get_activation(self) -> str:
        """Returns the name of the activation function."""
        return self._activation

    def get_biases(self) -> {str, tfk.initializers.Initializer, None}:
        """Returns the biases initializer."""
        return self._biases

    def get_kernel(self) -> {str, tfk.initializers.Initializer, None}:
        """Returns the kernel initializer."""
        return self._kernel

    def get_kernel_regularizer(self) -> {str, Regularizer, None}:
        """Returns the kernel regularizer."""
        return self._kernel_regularizer

    def get_name(self) -> str:
        """Returns the name of the layer."""
        return self._name

    def build(self) -> tfk.layers.Layer:
        return tfk.layers.Dense(
            self._n,
            activation=self._activation,
            kernel_initializer=self._kernel,
            bias_initializer=self._biases,
            kernel_regularizer=self._kernel_regularizer,
            name=self._name,
        )


class UnitRescaling(tf.keras.layers.Layer):
    n: int
    a: np.ndarray
    b: np.ndarray
    signed: bool
    invert: bool
    rescaling: Rescaling

    def __init__(
        self,
        n: int,
        a: np.ndarray,
        b: np.ndarray,
        signed: bool = False,
        invert: bool = False,
        name: str | None = None,
    ):
        super(UnitRescaling, self).__init__(trainable=False, name=name)
        self.n = n
        self.a = a
        self.b = b
        self.signed = signed
        self.invert = invert
        if signed:
            scale = 0.5 * (b - a)
            shift = 0.5 * (a + b)
        else:
            scale = b - a
            shift = a
        if not invert:
            shift = self.deshift(scale, shift)
            scale = self.descale(scale)
        self.rescaling = tf.keras.layers.Rescaling(scale, shift)

    @staticmethod
    def descale(scale):
        return 1.0 / scale

    @staticmethod
    def deshift(scale, shift):
        return -shift / scale

    def build(self, input_shape):
        self.rescaling.build(input_shape)

    def call(self, inputs, **kwargs):
        return self.rescaling.call(inputs)

    def get_config(self):
        return {
            "n": self.n,
            "a": self.a,
            "b": self.b,
            "signed": self.signed,
            "invert": self.invert,
            "name": self.name,
        }

    class Builder(LayerBuilder):
        """Class to build a unit-rescaling layer."""

        _n: int
        _a: np.ndarray
        _b: np.ndarray
        _signed: bool
        _invert: bool
        _name: str | None

        def __init__(self, n: int, a: np.ndarray, b: np.ndarray):
            """Creates a new instance of this class.

            The layer built will shift and scale its input into a distribution
            on the unit interval `[0, 1]`. It accomplishes this by calling `(
            input - a) / (b - a)` at runtime where `a` and `b` are the minimum
            and maximum values of the input.

            :param n: The number of neurons in the layer.
            :param a: The `n` minimum values.
            :param b: The `n` maximum values.
            """
            assert a.size == n
            assert b.size == n
            self._n = n
            self._a = a
            self._b = b
            self._signed = False
            self._invert = False
            self._name = None

        def get_n(self) -> int:
            return self._n

        def invert(self, invert: bool) -> object:
            self._invert = invert
            return self

        def signed(self, signed: bool) -> object:
            self._signed = signed
            return self

        def name(self, name: str) -> object:
            """Defines the name of the layer.

            :param name: The name of the layer.
            :return: itself.
            """
            self._name = name
            return self

        def build(self) -> tfk.layers.Layer:
            return UnitRescaling(
                self._n,
                self._a,
                self._b,
                self._signed,
                self._invert,
                name=self._name,
            )


class WeightInitializer(tfk.initializers.Initializer):
    _w: np.ndarray

    def __init__(self, w: np.ndarray):
        self._w = w

    def __call__(self, shape, dtype=None, **kwargs):
        return tf.constant(self._w, dtype=dtype, shape=shape)

    def get_config(self):
        return {"w": self._w}
