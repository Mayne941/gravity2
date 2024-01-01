import random, string

def generate_file_names(payload, ExpDir):
    '''Generate dictionary of file and folder names'''
    fnames = {}
    '''Main'''
    fnames["ExpDir"] = ExpDir
    fnames["OutputDir"] = f'{ExpDir}/output'

    '''Read Genome Desc Table'''
    fnames["ReadGenomeDescTablePickle"] = f'{fnames["OutputDir"]}/ReadGenomeDescTable.p'
    # If PL2 run (mandatory UCF file key is in payload)... make these too. << New fn

    '''PPHMMDB Cosntruction'''
    fnames = generate_pphmmdb_fnames(fnames)

    '''Ref Virus Annotator'''
    fnames = generate_ref_annotator_fnames(fnames)

    '''PL1 Graphs & MI'''
    fnames = generate_pl1_graph_fnames(fnames,payload)
    fnames = generate_mi_scorer_fnames(fnames)
    return fnames

def generate_pphmmdb_fnames(fnames):
    '''Generate file and folder names for PPHMMDB Constructor'''
    '''Mash Dirs'''
    fnames['MashDir'] = f"{fnames['ExpDir']}/Mash"
    fnames['MashQueryFile'] = f"{fnames['MashDir']}/Query.fasta"
    fnames['MashSubjectFile'] = f"{fnames['MashDir']}/Subjects.fasta"
    fnames['MashOutputFile'] = f"{fnames['MashDir']}/MashOutput.txt"
    fnames['MashSimFile'] = f"{fnames['MashDir']}/BitScoreMat.txt"
    fnames['MashProtClusterFile'] = f"{fnames['MashDir']}/ProtClusters.txt"
    fnames['ClustersDir'] = f"{fnames['MashDir']}/Clusters"
    '''HMMER Dirs'''
    fnames['HMMERDir'] = f"{fnames['ExpDir']}/HMMER"
    fnames['HMMER_PPHMMDir'] = f"{fnames['HMMERDir']}/HMMER_PPHMMs"
    fnames['HMMER_PPHMMDbDir'] = f"{fnames['HMMERDir']}/HMMER_PPHMMDb"
    fnames['HMMER_PPHMMDb'] = f"{fnames['HMMER_PPHMMDbDir']}/HMMER_PPHMMDb"
    '''HHSuite Dirs'''
    fnames['HHsuiteDir'] = f"{fnames['ExpDir']}/HHsuite"
    fnames['HHsuite_PPHMMDir'] = f"{fnames['HHsuiteDir']}/HHsuite_PPHMMs"
    fnames['HHsuite_PPHMMDBDir'] = f"{fnames['HHsuiteDir']}/HHsuite_PPHMMDB"
    fnames['HHsuite_PPHMMDB'] = f"{fnames['HHsuite_PPHMMDBDir']}/HHsuite_PPHMMDB"
    '''Misc'''
    fnames['RefSeqFile'] = f"{fnames['OutputDir']}/ref_seqs.fasta"
    fnames["PphmmdbPickle"] = f'{fnames["OutputDir"]}/PPHMMDBConstruction.p'
    # If PL2 run (mandatory UCF file key is in payload)... make these too.
    return fnames

def generate_ref_annotator_fnames(fnames):
    '''Generate file and folder names for Ref Virus Annotator'''
    '''HMMER Dirs'''
    fnames['HMMER_hmmscanDir'] = f"{fnames['HMMERDir']}/hmmscan_{''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))}"
    fnames['PPHMMDB_Summary'] = f"{fnames['HMMER_PPHMMDb']}_Summary.txt"
    fnames['PPHMMQueryFile'] = f"{fnames['HMMER_hmmscanDir']}/QProtSeqs.fasta"
    fnames['PPHMMScanOutFile'] = f"{fnames['HMMER_hmmscanDir']}/PPHMMScanOut.txt"
    '''Misc'''
    fnames['RefAnnotatorPickle'] = f'{fnames["OutputDir"]}/RefVirusAnnotator.p'
    return fnames

def generate_pl1_graph_fnames(fnames, payload):
    '''Generate filenames for graphical and tree output'''
    '''Mandatory files'''
    fnames['HeatmapWithDendrogramFile'] = f"{fnames['OutputDir']}/GRAViTy_heatmap.pdf"
    fnames['VirusGroupingFile'] = f"{fnames['OutputDir']}/virus_grouping.txt"
    '''Dendro files'''
    if payload['Bootstrap'] == True:
        fnames['VirusDendrogramDistFile'] = f"{fnames['OutputDir']}/dendrogram_dist.nwk"
        fnames['Heatmap_DendrogramFile'] = f"{fnames['OutputDir']}/dendrogram_for_heatmap.nwk"
        fnames['BootstrappedDendrogramFile'] = f"{fnames['OutputDir']}/dendrogram_bootstrapped.nwk"
    else:
        fnames['Heatmap_DendrogramFile'] = f"{fnames['OutputDir']}/dendrogram.nwk"
    return fnames

def generate_mi_scorer_fnames(fnames):
    fnames['MutualInformationScoreDir'] = fnames['OutputDir']+"/MutualInformationScore"
    fnames['MiScorePickle'] = f'{fnames["OutputDir"]}/MutualInformationCalculator.p'
    fnames['MutualInformationScoreFile'] = f"{fnames['MutualInformationScoreDir']}/MIScore.csv"
    return fnames
