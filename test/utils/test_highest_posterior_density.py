import pytest
import numpy as np
from app.utils.highest_posterior_density import hpd, calc_min_interval

def test_hpd_univariate():
    x = np.random.normal(0, 1, 1000)
    alpha = 0.05

    result = hpd(x, alpha)

    assert isinstance(result, np.ndarray)
    assert result.shape == (2,)
    assert result[0] < result[1]

def test_calc_min_interval():
    x = np.sort(np.random.normal(0, 1, 1000))
    alpha = 0.05

    result = calc_min_interval(x, alpha)

    assert isinstance(result, tuple)
    assert len(result) == 2
    assert result[0] < result[1]

if __name__ == "__main__":
    pytest.main()
