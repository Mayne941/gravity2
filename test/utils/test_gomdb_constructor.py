import pytest
import numpy as np
from app.utils.gomdb_constructor import GOMDB_Constructor

def test_GOMDB_Constructor():
    TaxoGroupingList = np.array(["A", "B", "A", "C"])
    PPHMMLocationTable = np.array([
        [1, 0, 1],
        [0, 1, 0],
        [1, 1, 1],
        [0, 0, 1]
    ])
    GOMIDList = ["A", "B", "C"]

    result = GOMDB_Constructor(TaxoGroupingList, PPHMMLocationTable, GOMIDList)

    expected_result = {
        "A": np.array([
            [1, 0, 1],
            [1, 1, 1]
        ]),
        "B": np.array([
            [0, 1, 0]
        ]),
        "C": np.array([
            [0, 0, 1]
        ])
    }

    for key in expected_result:
        assert np.array_equal(result[key], expected_result[key]), f"Expected {expected_result[key]} for key {key}, but got {result[key]}"

if __name__ == "__main__":
    pytest.main()
