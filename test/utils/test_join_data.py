import pytest
from unittest import mock
import pandas as pd
from app.utils.join_data import join_input

# FILE: app/utils/test_join_data.py


@mock.patch("app.utils.join_data.pd.read_csv")
@mock.patch("app.utils.join_data.pd.DataFrame.to_csv")
def test_join_input(mock_to_csv, mock_read_csv):
    # Mock data
    mock_df1 = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    mock_df2 = pd.DataFrame({"A": [5, 6], "B": [7, 8]})
    mock_read_csv.side_effect = [mock_df1, mock_df2]

    payload = {
        "vmr_1": "path/to/vmr_1.csv",
        "vmr_2": "path/to/vmr_2.csv",
        "vmr_joined": "path/to/vmr_joined.csv"
    }

    result = join_input(payload)

    # Assertions
    mock_read_csv.assert_any_call("path/to/vmr_1.csv")
    mock_read_csv.assert_any_call("path/to/vmr_2.csv")
    assert mock_to_csv.called
    assert result == "CSV join completed successfully. Output saved to path/to/vmr_joined.csv"
