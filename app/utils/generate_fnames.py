import random, string

def generate_file_names(payload, ExpDir, Pl2=False):
    '''Generate dictionary of file and folder names. If called from PL1, ignore
    PL2 args; if called from PL2, all vars prior to if statement pertain to PL2,
    then we generate dummy PL1 vars from the extra PL2 arguments.'''
    fnames = {}
    output_prefix = "/output"

    '''Main'''
    fnames["ExpDir"] = ExpDir
    fnames["OutputDir"] = f'{ExpDir}{output_prefix}'

    '''Read Genome Desc Table'''
    fnames["ReadGenomeDescTablePickle"] = f'{fnames["OutputDir"]}/ReadGenomeDescTable.p'

    '''PPHMMDB Cosntruction'''
    fnames = generate_pphmmdb_fnames(fnames)

    '''Ref Virus Annotator'''
    fnames = generate_ref_annotator_fnames(fnames)

    '''PL1 Graphs & MI'''
    fnames = generate_pl1_graph_fnames(fnames,payload)
    fnames = generate_mi_scorer_fnames(fnames)
    fnames = generate_shared_pl_graph_fnames(fnames)
    '''PL2 components'''
    if Pl2:
        '''Main'''
        fnames = generate_dummy_pl1_fnames(fnames, payload, output_prefix)
        '''UCF Annotator'''
        fnames = generate_ucf_annotator_fnames(fnames)
        '''Virus Clasifier'''
        fnames = generate_virus_classifier_fnames(fnames)

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
    return fnames

def generate_ref_annotator_fnames(fnames):
    '''Generate file and folder names for Ref Virus Annotator'''
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

def generate_shared_pl_graph_fnames(fnames):
    '''Generate fnames for graphs comparing PPHMM/GOM content, for both pipelines'''
    fnames["PphmmAndGomSigs"] = f"{fnames['OutputDir']}/PPHMMandGOMsignatures.csv"
    fnames["SharedPphmmHm"] = f"{fnames['OutputDir']}/shared_pphmms.png"
    fnames["NormPphmmRatioHm"] = f"{fnames['OutputDir']}/shared_norm_pphmm_ratio.png"
    fnames["PphmmLocs"] = f"{fnames['OutputDir']}/PPHMMLocs.csv"
    fnames["PphmmLocDists"] = f"{fnames['OutputDir']}/pphmm_location_distances.png"
    fnames["PphmmLocPairwise"] = f"{fnames['OutputDir']}/pphmm_location_distances_pairwise.png"
    return fnames

def generate_mi_scorer_fnames(fnames):
    '''Generate filenames for mutual information calculator'''
    fnames['MutualInformationScoreDir'] = fnames['OutputDir']+"/MutualInformationScore"
    fnames['MiScorePickle'] = f'{fnames["OutputDir"]}/MutualInformationCalculator.p'
    fnames['MutualInformationScoreFile'] = f"{fnames['MutualInformationScoreDir']}/MIScore.csv"
    return fnames

def generate_dummy_pl1_fnames(fnames, payload, output_pre):
    '''Generate filenames for PL1 reference data, when running PL2'''
    fnames["ExpDir_Pl1"] = payload['ExpDir_Pl1']
    fnames["Pl1OutputDir"] = f'{fnames["ExpDir_Pl1"]}{output_pre}'
    fnames['Pl1ReadDescTablePickle'] = f'{fnames["Pl1OutputDir"]}/ReadGenomeDescTable.p'
    fnames['Pl1RefAnnotatorPickle'] = f'{fnames["Pl1OutputDir"]}/RefVirusAnnotator.p'
    return fnames

def generate_ucf_annotator_fnames(fnames):
    '''Generate filenames for unclassified virus annotator'''
    fnames['HMMERDir_RefVirus'] = f"{fnames['ExpDir_Pl1']}/HMMER"
    fnames['HMMER_PPHMMDB_RefVirus'] = f"{fnames['HMMERDir_RefVirus']}/HMMER_PPHMMDb/HMMER_PPHMMDb"
    fnames['UcfAnnotatorPickle'] = f'{fnames["OutputDir"]}/UcfVirusAnnotator.p'
    return fnames

def generate_virus_classifier_fnames(fnames):
    '''Generate filenames for virus classifier'''
    fnames['HMMERDir_UcfVirus'] = f"{fnames['ExpDir']}/HMMER"
    fnames['HMMER_PPHMMDBDir_UcfVirus'] = f"{fnames['HMMERDir_UcfVirus']}/HMMER_PPHMMDb"
    fnames['HMMER_PPHMMDB_UcfVirus'] = f"{fnames['HMMER_PPHMMDBDir_UcfVirus']}/HMMER_PPHMMDb"
    fnames['VirusDendrogramFile'] = f"{fnames['OutputDir']}/dendrogram.nwk"
    fnames['DendrogramDist'] = f"{fnames['OutputDir']}/dendrogram_dist.nwk"
    fnames['BootstrappedDendrogramFile'] = f"{fnames['OutputDir']}/dendrogram_bootstrapped.nwk"
    fnames['HeatmapWithDendrogramFile'] = f"{fnames['OutputDir']}/GRAViTy_heatmap.pdf"
    fnames['ClassificationResultFile'] = f"{fnames['OutputDir']}/ClassificationResults.txt"
    fnames['VirusGroupingFile'] = f"{fnames['OutputDir']}/VirusGrouping.txt"
    fnames['VirusClassifierPickle'] = f'{fnames["OutputDir"]}/VirusClassification.p'
    return fnames
