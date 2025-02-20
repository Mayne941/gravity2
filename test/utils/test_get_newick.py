import pytest
from unittest.mock import MagicMock
from app.utils.get_newick import GetNewick

def test_GetNewick():
    node = MagicMock()
    node.is_leaf.side_effect = [False, True, True]
    node.get_left.return_value = node
    node.get_right.return_value = node
    node.dist = 0.5
    node.id = 0

    newick = ""
    parentdist = 1.0
    leaf_names = ["A", "B"]

    result = GetNewick(node, newick, parentdist, leaf_names)

    expected_result = "(A:0.000000,A:0.000000);"

    assert result == expected_result

if __name__ == "__main__":
    pytest.main()
