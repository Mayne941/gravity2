import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.gom_signature_table_constructor import generate_gom_sigs

@patch("app.utils.gom_signature_table_constructor.Pool")
@patch("app.utils.gom_signature_table_constructor.tqdm")
@patch("app.utils.gom_signature_table_constructor.progress_msg")
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
    bootstrap = 0

    mock_pool.return_value.__enter__.return_value.apply_async.side_effect = [
        MagicMock(get=MagicMock(return_value=[0.5, 0.6, 0.7])),
        MagicMock(get=MagicMock(return_value=[0.2, 0.3, 0.4]))
    ]
    mock_tqdm.return_value.__enter__.return_value = MagicMock()

    result = GOMSignatureTable_Constructor(PPHMMLocationTable, GOMDB, GOMIDList, bootstrap)

    expected_result = np.array([
        [0.5, 0.2],
        [0.6, 0.3],
        [0.7, 0.4]
    ])

    assert np.allclose(result, expected_result)
    mock_progress_msg.assert_called_once()
    mock_pool.assert_called_once()
    mock_tqdm.assert_called_once()

def test_generate_gom_sigs():
    i = "GOM1"
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

    result = generate_gom_sigs(i, PPHMMLocationTable, GOMDB)

    expected_result = [1.0, 1.0, 0.0]

    assert np.allclose(result, expected_result)

if __name__ == "__main__":
    pytest.main()
