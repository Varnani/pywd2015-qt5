"""Validated, column-oriented observation input shared by all loaders."""
import warnings
import numpy


def read_columns(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        rows = numpy.loadtxt(path, ndmin=2)
    if rows.shape[0] == 0 or rows.shape[1] < 2:
        raise ValueError("Observation file requires at least one row and two columns: " + str(path))
    if not numpy.isfinite(rows).all():
        raise ValueError("Observation file contains NaN or infinity: " + str(path))
    return rows.T


def select_columns(data, x_index, y_index, weight_index=None):
    indices = (x_index, y_index) if weight_index is None else (x_index, y_index, weight_index)
    if any(not isinstance(i, int) or i < 0 or i >= len(data) for i in indices):
        raise ValueError("Selected observation column does not exist.")
    weights = numpy.ones(data.shape[1]) if weight_index is None else data[weight_index]
    if (weights < 0).any():
        raise ValueError("Observation weights must be nonnegative.")
    return data[x_index], data[y_index], weights


def read_timings(path, constant_weight=False):
    data = read_columns(path)
    if len(data) not in (2, 3):
        raise ValueError("Eclipse timings require two or three columns.")
    if not numpy.isin(data[1], (1, 2)).all():
        raise ValueError("Eclipse type must be 1 or 2.")
    return numpy.array(select_columns(data, 0, 1, None if constant_weight or len(data) == 2 else 2))


def validate_interval(start, stop, step):
    import math
    if not all(math.isfinite(value) for value in (start, stop, step)) or step <= 0 or stop < start:
        raise ValueError("Require finite ordered boundaries and a positive sampling step.")
    return math.floor((stop - start) / step) + 1
