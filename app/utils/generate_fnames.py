
def generate_file_names(payload, ExpDir):
    '''Generate dictionary of file and folder names'''
    fnames = {}
    '''Main'''
    fnames["ExpDir"] = ExpDir
    fnames["OutputDir"] = f'{ExpDir}/output'

    '''Read Genome Desc Table'''
    fnames["ReadGenomeDescTablePickle"] = f'{fnames["OutputDir"]}/ReadGenomeDescTable.p'

    '''PPHMMDB Cosntruction'''
    fnames = generate_pphmmdb_fnames(fnames)

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
