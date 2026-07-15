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

from .netcodec import SnnCodec

tfk = tf.keras


class NetReader(ABC):
    """The neural network reader interface class."""

    @abstractmethod
    def read(self, path: str) -> tfk.Model:
        """Reads a neural network from a file.

        :param path: The path to the neural network file.
        :return: The neural network read.
        """


class SnnReader(NetReader):
    """Reads a custom Stuttgart Neural Network (SNN) Simulator model from a
    file. Helmut Schiller (deceased) used a streamlined SNN-based model for
    his feed-forward back-propagation (FFBP) neural network software. He
    developed a custom file format for storing his models, too, which is
    different from the original SNN format. We yet use the term SNN to refer
    to FFBP neural networks and file format.
    """

    def read(self, path: str) -> tfk.Sequential:
        return SnnCodec().decode(SnnParser().parse(SnnReader._read(path)))

    @staticmethod
    def _read(path):
        with open(path) as f:
            lines = SnnReader._read_lines(SnnReader._skip_header(f))
        return lines

    @staticmethod
    def _read_lines(f):
        lines = []
        line = f.readline()
        while line:
            line = line.strip()
            if line:
                lines.append(line)
            line = f.readline()
        return lines

    @staticmethod
    def _skip_header(f):
        line = f.readline()
        while line:
            line = line.strip()
            if line.startswith("#"):
                break
            line = f.readline()
        return f


class SnnParser:
    """Parses an SNN model from text lines."""

    _config: dict[str, Any]

    def __init__(self, config: dict[str, Any] = None):
        """Creates a new parser with a given pre-configuration.

        :param config: The pre-configuration.
        """
        if config is None:
            self._config = {}
        else:
            self._config = config

    def parse(self, lines: list[str]) -> dict:
        """Parses an SNN model from text lines into an encoding dictionary.

        :param lines: The text lines
        :return: The encoding dictionary.
        """
        enc: dict[str, Any] = self._config.copy()
        SnnParser._parse_n_inputs(lines, enc)
        SnnParser._parse_b_inputs(lines, enc)
        SnnParser._parse_n_output(lines, enc)
        SnnParser._parse_b_output(lines, enc)
        SnnParser._parse_n_layers(lines, enc)
        SnnParser._parse_biases(lines, enc)
        SnnParser._parse_kernel(lines, enc)
        return enc

    @staticmethod
    def _parse_kernel(lines, enc):
        kernel = []
        shapes = []
        n_layers = enc["n_layers"]
        biases = enc["biases"]
        for i in range(n_layers - 2):
            popped = SnnParser._parse_mark(lines, "wgt ")
            splits = popped.split()
            assert len(splits) == 3
            m = int(splits[1])
            n = int(splits[2])
            assert n == len(biases[i])
            w = SnnParser._parse_values(lines, n * m)
            kernel.append(w)
            shapes.append((n, m))
        enc["kernel"] = kernel
        enc["shapes"] = shapes

    @staticmethod
    def _parse_biases(lines, enc):
        biases = []
        n_layers = enc["n_layers"]
        for i in range(n_layers - 2):
            popped = SnnParser._parse_mark(lines, "bias ")
            splits = popped.split()
            assert len(splits) == 2
            n = int(splits[1])
            b = SnnParser._parse_values(lines, n)
            biases.append(b)
        enc["biases"] = biases

    @staticmethod
    def _parse_values(lines, n):
        x = np.zeros(n, np.double)
        for i in range(len(x)):
            popped = lines.pop(0)
            x[i] = np.double(popped)
        return x

    @staticmethod
    def _parse_n_layers(lines, enc):
        n_inputs = enc["n_inputs"]
        n_output = enc["n_output"]
        # noinspection PyUnusedLocal
        popped = SnnParser._parse_mark(lines, "$")
        popped = SnnParser._parse_mark(lines, "#planes=")
        splits = popped.split()
        n_splits = len(splits)
        planes = np.zeros(n_splits, int)
        for i in range(n_splits):
            planes[i] = int(splits[i])
        assert planes[0] == n_splits - 1
        assert planes[1] == n_inputs
        assert planes[planes[0]] == n_output
        n_layers = planes[0] + 1
        enc["n_layers"] = n_layers

    @staticmethod
    def _parse_mark(lines, mark):
        popped = lines.pop(0)
        assert popped.startswith(mark)
        return popped.replace(mark, "")

    @staticmethod
    def _parse_b_output(lines, enc):
        n_output = enc["n_output"]
        b_output = np.zeros((n_output, 2), np.double)
        SnnParser._parse_bounds(lines, b_output)
        enc["b_output"] = b_output

    @staticmethod
    def _parse_n_output(lines, enc):
        popped = lines.pop(0)
        n_output = int(popped)
        enc["n_output"] = n_output

    @staticmethod
    def _parse_b_inputs(lines, enc):
        n_inputs = enc["n_inputs"]
        b_inputs = np.zeros((n_inputs, 2), np.double)
        SnnParser._parse_bounds(lines, b_inputs)
        enc["b_inputs"] = b_inputs

    @staticmethod
    def _parse_bounds(lines, bounds):
        for b in bounds:
            popped = lines.pop(0)
            splits = popped.split()
            assert len(splits) == 2
            b[0] = np.double(splits[0])
            b[1] = np.double(splits[1])

    @staticmethod
    def _parse_n_inputs(lines, enc):
        popped = lines.pop(0)
        n_inputs = int(popped)
        enc["n_inputs"] = n_inputs
