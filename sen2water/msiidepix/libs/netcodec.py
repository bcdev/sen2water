#  Copyright (c) 2023 EUMETSAT (contract EUM/CO/19/4600002254/JCh)
#
#  Code authored 2023 by Brockmann Consult GmbH under EUMETSAT Contract
#  EUM/CO/19/4600002254/JCh "Sentinel-3 Synergy Cloud Mask Development".

from abc import ABC
from abc import abstractmethod
from typing import Any

# noinspection PyPackageRequirements
import numpy as np

# noinspection PyPackageRequirements
import tensorflow as tf

from .netbuilder import SnnBuilder

tfk = tf.keras


class NetCodec(ABC):

    @abstractmethod
    def decode(self, enc: dict[str, Any]) -> tfk.Model:
        """Decodes a dictionary into a model.

        :param enc: The encoding dictionary.
        :return: the model.
        """

    @abstractmethod
    def encode(self, model: tfk.Model) -> dict[str, Any]:
        """Encodes a model into a dictionary.

        :param model: The model.
        :return: the encoding dictionary.
        """


class SnnCodec(NetCodec):
    """Codec for custom Stuttgart Neural Network (SNN) Simulator models.
    Helmut Schiller (deceased) used a streamlined SNN-based model for his
    feed-forward back-propagation (FFBP) neural network software. He
    developed a custom file format for storing his models, too, which is
    different from the original SNN format. We yet use the term SNN to refer
    to FFBP neural networks and file format.
    """

    def decode(self, enc: dict[str, Any]) -> tfk.Sequential:
        builder = SnnBuilder()
        SnnCodec._decode_n_inputs(enc, builder)
        SnnCodec._decode_b_inputs(enc, builder)
        SnnCodec._decode_w_layers(enc, builder)
        SnnCodec._decode_b_output(enc, builder)
        return builder.build()

    def encode(self, model: tfk.Model) -> dict[str, Any]:
        enc: dict[str, Any] = {}
        SnnCodec._encode_b_inputs(model, enc)
        SnnCodec._encode_b_output(model, enc)
        SnnCodec._encode_w_layers(model, enc)
        return enc

    @staticmethod
    def _decode_b_inputs(enc, builder):
        n_inputs = enc["n_inputs"]
        b_inputs = enc["b_inputs"]
        assert len(b_inputs) == n_inputs
        builder.add_rescaling_unit(n_inputs, b_inputs.transpose())

    @staticmethod
    def _decode_b_output(enc, builder):
        n_output = enc["n_output"]
        b_output = enc["b_output"]
        assert len(b_output) == n_output
        builder.add_rescaling_unit(n_output, b_output.transpose(), invert=True)

    @staticmethod
    def _decode_n_inputs(enc, builder):
        n_inputs = enc["n_inputs"]
        builder.add_input(n_inputs)

    @staticmethod
    def _decode_w_layers(enc, builder):
        n_layers = enc["n_layers"]
        biases = enc["biases"]
        kernel = enc["kernel"]
        shapes = enc["shapes"]
        assert len(biases) == n_layers - 2
        assert len(kernel) == n_layers - 2
        assert len(shapes) == n_layers - 2
        for i in range(n_layers - 2):
            n, m = shapes[i]
            builder.add_dense(n, biases=biases[i], kernel=kernel[i])

    @staticmethod
    def _encode_b_inputs(model, enc):
        layers_config = model.get_config()["layers"]
        config = layers_config[0]
        assert config["class_name"] == "InputLayer"
        for config in layers_config:
            if config["class_name"] == "UnitRescaling":
                c = config["config"]
                assert c["signed"] is False
                assert c["invert"] is False
                n = c["n"]
                a = c["a"]
                b = c["b"]
                bounds = np.array([a, b]).transpose().reshape(n, 2)
                enc["n_inputs"] = n
                enc["b_inputs"] = bounds
                break
        assert config["class_name"] == "UnitRescaling"

    @staticmethod
    def _encode_b_output(model, enc):
        layers_config = model.get_config()["layers"]
        config = layers_config[-1]
        assert config["class_name"] == "UnitRescaling"
        config = config["config"]
        assert config["signed"] is False
        assert config["invert"] is True
        n = config["n"]
        a = config["a"]
        b = config["b"]
        bounds = np.array([a, b]).transpose().reshape(n, 2)
        enc["n_output"] = n
        enc["b_output"] = bounds

    @staticmethod
    def _encode_w_layers(model, enc):
        layers_config = model.get_config()["layers"]
        n_inputs = enc["n_inputs"]
        n_planes = 1
        for i, config in enumerate(layers_config):
            if config["class_name"] == "Dense":
                n_planes += 1
        biases = []
        kernel = []
        shapes = []
        planes = [n_planes, n_inputs]
        w = model.get_weights()
        for i in range(0, len(w), 2):
            shapes.append(w[i].shape)
            kernel.append(w[i].transpose().reshape(w[i].size))
            biases.append(w[i + 1])
            planes.append(w[i + 1].size)
        enc["biases"] = biases
        enc["kernel"] = kernel
        enc["shapes"] = shapes
        enc["planes"] = planes
