import pytest
import numpy as np
from app.utils.dist_mat_to_tree import DistMat2Tree

def test_DistMat2Tree():
    DistMat = np.array([
        [0.0, 0.2, 0.4],
        [0.2, 0.0, 0.6],
        [0.4, 0.6, 0.0]
    ])
    LeafList = ["A", "B", "C"]
    Dendrogram_LinkageMethod = "single"
    do_logscale = False

    result = DistMat2Tree(DistMat, LeafList, Dendrogram_LinkageMethod, do_logscale)

    expected_result = "((B:0.200000,A:0.200000):0.200000,C:0.400000);"

    assert result == expected_result
