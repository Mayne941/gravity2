import pytest
from app.utils.generate_fnames import generate_file_names

def test_generate_file_names():
    payload = {
        "ExperimentName": "test_experiment",
        "GenomeDescTableFile_FirstPass": "first_pass_desc_table.txt",
        "GenomeDescTableFile_SecondPass": "second_pass_desc_table.txt",
        "Bootstrap": True,
        "ExpDir_Pl1": "/path/to/pl1"
    }
    ExpDir = "/path/to/expdir"
    Pl2 = True

    fnames = generate_file_names(payload, ExpDir, Pl2)

    output_prefix = "/output"

    # Assertions for main directories
    assert fnames["ExpDir"] == ExpDir
    assert fnames["OutputDir"] == f"{ExpDir}{output_prefix}"

    # Assertions for Read Genome Desc Table
    assert fnames["ReadGenomeDescTablePickle"] == f"{fnames['OutputDir']}/ReadGenomeDescTable.p"

    # Assertions for PPHMMDB Construction
    assert fnames['MashDir'] == f"{fnames['ExpDir']}/Mash"
    assert fnames['MashQueryFile'] == f"{fnames['MashDir']}/Query.fasta"
    assert fnames['MashSubjectFile'] == f"{fnames['MashDir']}/Subjects.fasta"
    assert fnames['MashOutputFile'] == f"{fnames['MashDir']}/MashOutput.txt"
    assert fnames['MashSimFile'] == f"{fnames['MashDir']}/BitScoreMat.txt"
    assert fnames['MashProtClusterFile'] == f"{fnames['MashDir']}/ProtClusters.txt"
    assert fnames['ClustersDir'] == f"{fnames['MashDir']}/Clusters"
    assert fnames['HMMERDir'] == f"{fnames['ExpDir']}/HMMER"
    assert fnames['HMMER_PPHMMDir'] == f"{fnames['HMMERDir']}/HMMER_PPHMMs"
    assert fnames['HMMER_PPHMMDbDir'] == f"{fnames['HMMERDir']}/HMMER_PPHMMDb"
    assert fnames['HMMER_PPHMMDb'] == f"{fnames['HMMER_PPHMMDbDir']}/HMMER_PPHMMDb"
    assert fnames['HHsuiteDir'] == f"{fnames['ExpDir']}/HHsuite"
    assert fnames['HHsuite_PPHMMDir'] == f"{fnames['HHsuiteDir']}/HHsuite_PPHMMs"
    assert fnames['HHsuite_PPHMMDBDir'] == f"{fnames['HHsuiteDir']}/HHsuite_PPHMMDB"
    assert fnames['HHsuite_PPHMMDB'] == f"{fnames['HHsuite_PPHMMDBDir']}/HHsuite_PPHMMDB"
    assert fnames['RefSeqFile'] == f"{fnames['OutputDir']}/ref_seqs.fasta"
    assert fnames["PphmmdbPickle"] == f'{fnames["OutputDir"]}/PPHMMDBConstruction.p'

    # Assertions for Ref Virus Annotator
    assert fnames['RefAnnotatorPickle'] == f'{fnames["OutputDir"]}/RefVirusAnnotator.p'

    # Assertions for PL1 Graphs
    assert fnames['HeatmapWithDendrogramFile'] == f"{fnames['OutputDir']}/GRAViTy_heatmap.pdf"
    assert fnames['VirusGroupingFile'] == f"{fnames['OutputDir']}/VirusGrouping.txt"
    assert fnames['VirusDendrogramDistFile'] == f"{fnames['OutputDir']}/dendrogram_dist.nwk"
    assert fnames['Heatmap_DendrogramFile'] == f"{fnames['OutputDir']}/dendrogram_for_heatmap.nwk"
    assert fnames['BootstrappedDendrogramFile'] == f"{fnames['OutputDir']}/dendrogram_bootstrapped.nwk"

    # Assertions for shared PL graphs
    assert fnames["PphmmAndGomSigs"] == f"{fnames['OutputDir']}/PPHMMandGOMsignatures.csv"
    assert fnames["SharedPphmmHm"] == f"{fnames['OutputDir']}/shared_pphmms.png"
    assert fnames["NormPphmmRatioHm"] == f"{fnames['OutputDir']}/shared_norm_pphmm_ratio.png"
    assert fnames["PphmmLocs"] == f"{fnames['OutputDir']}/PPHMMLocs.csv"
    assert fnames["PphmmLocDists"] == f"{fnames['OutputDir']}/pphmm_location_distances.png"
    assert fnames["PphmmLocPairwise"] == f"{fnames['OutputDir']}/pphmm_location_distances_pairwise.png"

    # Assertions for MI scorer
    assert fnames['MutualInformationScoreDir'] == fnames['OutputDir']+"/MutualInformationScore"
    assert fnames['MiScorePickle'] == f'{fnames["OutputDir"]}/MutualInformationCalculator.p'
    assert fnames['MutualInformationScoreFile'] == f"{fnames['MutualInformationScoreDir']}/MIScore.csv"

    # Assertions for PL2 components
    assert fnames["ExpDir_Pl1"] == payload['ExpDir_Pl1']
    assert fnames["Pl1OutputDir"] == f'{fnames["ExpDir_Pl1"]}{output_prefix}'
    assert fnames['Pl1ReadDescTablePickle'] == f'{fnames["Pl1OutputDir"]}/ReadGenomeDescTable.p'
    assert fnames['Pl1RefAnnotatorPickle'] == f'{fnames["Pl1OutputDir"]}/RefVirusAnnotator.p'
    assert fnames['HMMERDir_RefVirus'] == f"{fnames['ExpDir_Pl1']}/HMMER"
    assert fnames['HMMER_PPHMMDB_RefVirus'] == f"{fnames['HMMERDir_RefVirus']}/HMMER_PPHMMDb/HMMER_PPHMMDb"
    assert fnames['UcfAnnotatorPickle'] == f'{fnames["OutputDir"]}/UcfVirusAnnotator.p'
    assert fnames['HMMERDir_UcfVirus'] == f"{fnames['ExpDir']}/HMMER"
    assert fnames['HMMER_PPHMMDBDir_UcfVirus'] == f"{fnames['HMMERDir_UcfVirus']}/HMMER_PPHMMDb"
    assert fnames['HMMER_PPHMMDB_UcfVirus'] == f"{fnames['HMMER_PPHMMDBDir_UcfVirus']}/HMMER_PPHMMDb"
    assert fnames['VirusDendrogramFile'] == f"{fnames['OutputDir']}/dendrogram.nwk"
    assert fnames['DendrogramDist'] == f"{fnames['OutputDir']}/dendrogram_dist.nwk"
    assert fnames['BootstrappedDendrogramFile'] == f"{fnames['OutputDir']}/dendrogram_bootstrapped.nwk"
    assert fnames['HeatmapWithDendrogramFile'] == f"{fnames['OutputDir']}/GRAViTy_heatmap.pdf"
    assert fnames['ClassificationResultFile'] == f"{fnames['OutputDir']}/ClassificationResults.txt"
    assert fnames['VirusGroupingFile'] == f"{fnames['OutputDir']}/VirusGrouping.txt"
    assert fnames['VirusClassifierPickle'] == f'{fnames["OutputDir"]}/VirusClassification.p'

if __name__ == "__main__":
    pytest.main()
