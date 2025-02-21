import pytest
import numpy as np
from app.utils.dcor import dcov
from app.utils.dcor import dvar
from app.utils.dcor import cent_dist

def test_dcov():
    X = np.array([[1, 2], [3, 4], [5, 6]])
    Y = np.array([[1, 4], [2, 5], [3, 6]])

    result = dcov(X, Y)

    expected_result = np.sqrt(np.multiply(X, Y).sum()) / X.shape[0]

    assert np.isclose(result, expected_result), f"Expected {expected_result}, but got {result}"


def test_dvar():
    X = np.array([[1, 2], [3, 4], [5, 6]])

    result = dvar(X)

    expected_result = np.sqrt(np.sum(X ** 2) / X.shape[0] ** 2)

    assert np.isclose(result, expected_result), f"Expected {expected_result}, but got {result}"

def test_cent_dist():
    X = np.array([[1, 2], [3, 4], [5, 6]])

    result = cent_dist(X)

    M = np.array([[0.0, 2.82842712, 5.65685425],
                  [2.82842712, 0.0, 2.82842712],
                  [5.65685425, 2.82842712, 0.0]])
    rmean = M.mean(axis=1)
    cmean = M.mean(axis=0)
    gmean = rmean.mean()
    R = np.tile(rmean, (M.shape[0], 1)).transpose()
    C = np.tile(cmean, (M.shape[1], 1))
    G = np.tile(gmean, M.shape)
    expected_result = M - R - C + G

    assert np.allclose(result, expected_result), f"Expected {expected_result}, but got {result}"
