import numpy as np


def test_euler():
    assert np.allclose(np.exp(1j * np.pi) + 1, 0)
