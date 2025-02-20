import pytest
import numpy as np
from ete3 import Tree

from app.utils.classification_utils import PairwiseSimilarityScore_Cutoff_Dict_Constructor
from app.utils.classification_utils import TaxonomicAssignmentProposerAndEvaluator

def test_PairwiseSimilarityScore_Cutoff_Dict_Constructor():
    SimMat = np.array([
        [1.0, 0.8, 0.3],
        [0.8, 1.0, 0.4],
        [0.3, 0.4, 1.0]
    ])
    TaxoGroupingList = ["A", "A", "B"]
    N_PairwiseSimilarityScores = 2

    result = PairwiseSimilarityScore_Cutoff_Dict_Constructor(SimMat, TaxoGroupingList, N_PairwiseSimilarityScores)

    assert "A" in result
    assert "B" in result
    assert "PairwiseSimilarityScore_IntraClass_List" in result["A"]
    assert "PairwiseSimilarityScore_InterClass_List" in result["A"]
    assert "CutOff" in result["A"]
    assert "PairwiseSimilarityScore_IntraClass_List" in result["B"]
    assert "PairwiseSimilarityScore_InterClass_List" in result["B"]
    assert "CutOff" in result["B"]

    assert len(result["A"]["PairwiseSimilarityScore_IntraClass_List"]) <= N_PairwiseSimilarityScores
    assert len(result["A"]["PairwiseSimilarityScore_InterClass_List"]) <= N_PairwiseSimilarityScores
    assert len(result["B"]["PairwiseSimilarityScore_IntraClass_List"]) <= N_PairwiseSimilarityScores
    assert len(result["B"]["PairwiseSimilarityScore_InterClass_List"]) <= N_PairwiseSimilarityScores

    assert isinstance(result["A"]["CutOff"], float)
    assert isinstance(result["B"]["CutOff"], float)


def test_TaxonomicAssignmentProposerAndEvaluator():
    SimMat_UcfVirusesVSRefViruses = np.array([
        [0.9, 0.2, 0.1],
        [0.3, 0.8, 0.4],
        [0.2, 0.3, 0.7]
    ])
    TaxoGrouping = np.array(["A", "B", "C"])
    VirusDendrogram = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
    TaxoLabelList_RefVirus = ["A", "B", "C"]
    TaxoLabelList_UcfVirus = ["D", "E", "F"]
    PairwiseSimilarityScore_Cutoff_Dict = {
        "A": {"CutOff": 0.5},
        "B": {"CutOff": 0.5},
        "C": {"CutOff": 0.5}
    }

    MaxSimScoreList, TaxoOfMaxSimScoreList, TaxoAssignmentList, PhyloStatList = TaxonomicAssignmentProposerAndEvaluator(
        SimMat_UcfVirusesVSRefViruses, TaxoGrouping, VirusDendrogram, TaxoLabelList_RefVirus, TaxoLabelList_UcfVirus, PairwiseSimilarityScore_Cutoff_Dict)

    assert len(MaxSimScoreList) == 3
    assert len(TaxoOfMaxSimScoreList) == 3
    assert len(TaxoAssignmentList) == 1
    assert len(PhyloStatList) == 1

    assert TaxoAssignmentList[0] == "A"
    assert PhyloStatList[0] == "4"  # D is embedded within a clade of A
