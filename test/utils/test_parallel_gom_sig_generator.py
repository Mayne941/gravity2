import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from app.utils.parallel_gom_sig_generator import GOMSignatureTable_Constructor
from app.utils.parallel_gom_sig_generator import generate_gom_sigs

@patch("app.utils.parallel_gom_sig_generator.Pool")
@patch("app.utils.parallel_gom_sig_generator.tqdm")
@patch("app.utils.parallel_gom_sig_generator.progress_msg")
def test_GOMSignatureTable_Constructor(mock_progress_msg, mock_tqdm, mock_pool):
    PPHMMLocationTable = np.array([
        [1, 0, 1],
        [0, 1, 0],
        [1, 1, 1]
    ])
    GOMDB = {
        "GOM1": np.array([
            [1, 0, 1],
            [0, 1, 0],
            [1, 1, 1]
        ]),
        "GOM2": np.array([
            [0, 1, 0],
            [1, 0, 1],
            [0, 0, 0]
        ])
    }
    GOMIDList = ["GOM1", "GOM2"]
    payload = {"N_CPUs": 2}
    bootstrap = 0

    mock_pool.return_value.__enter__.return_value.apply_async.side_effect = [
        MagicMock(get=MagicMock(return_value=[0.5, 0.6, 0.7])),
        MagicMock(get=MagicMock(return_value=[0.2, 0.3, 0.4]))
    ]
    mock_tqdm.return_value.__enter__.return_value = MagicMock()

    result = GOMSignatureTable_Constructor(PPHMMLocationTable, GOMDB, GOMIDList, payload, bootstrap)

    expected_result = np.array([
        [0.5, 0.2],
        [0.6, 0.3],
        [0.7, 0.4]
    ])

    assert np.allclose(result, expected_result), f"Expected {expected_result}, but got {result}"
    mock_progress_msg.assert_called_once()
    mock_pool.assert_called_once()
    mock_tqdm.assert_called_once()

def test_generate_gom_sigs():
    GOM = "GOM1"
    PPHMMLocationTable = np.array([
        [1, 0, 1],
        [0, 1, 0],
        [1, 1, 1]
    ])
    GOMDB = {
        "GOM1": np.array([
            [1, 0, 1],
            [0, 1, 0],
            [1, 1, 1]
        ]),
        "GOM2": np.array([
            [0, 1, 0],
            [1, 0, 1],
            [0, 0, 0]
        ])
    }

    result = generate_gom_sigs(GOM, PPHMMLocationTable, GOMDB)

    expected_result = [0.5, 0.5, 1.0]

    assert isinstance(result, list)
    assert len(result) == len(PPHMMLocationTable)
    assert all(isinstance(x, float) for x in result)

if __name__ == "__main__":
    pytest.main()
