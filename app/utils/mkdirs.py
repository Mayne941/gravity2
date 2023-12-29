import os
import random
import string
from app.utils.shell_cmds import shell

def mkdir_pphmmdbc(fnames):
    '''Return all directories for db storage'''
    '''Blast dirs'''
    if os.path.exists(fnames['MashDir']):
        '''Clear previous results'''
        shell(f"rm -rf {fnames['MashDir']}")
    os.makedirs(fnames['MashDir'])
    os.makedirs(fnames['ClustersDir'])

    '''HMMER dirs (+ delete existing HMMER libraries)'''
    if os.path.exists(fnames['HMMERDir']):
        '''Clear previous results'''
        shell(f"rm -rf {fnames['HMMERDir']}")
    os.makedirs(fnames['HMMERDir'])
    os.makedirs(fnames['HMMER_PPHMMDir'])
    os.makedirs(fnames['HMMER_PPHMMDbDir'])

    '''HHsuite dirs, regardless of if it's enabled'''
    if os.path.exists(fnames['HHsuiteDir']):
        '''Clear previous results'''
        shell(f"rm -rf {fnames['HHsuiteDir']}")
    os.makedirs(fnames['HHsuiteDir'])
    os.makedirs(fnames['HHsuite_PPHMMDir'])
    os.makedirs(fnames['HHsuite_PPHMMDBDir'])

def mkdir_ref_annotator(ExpDir, RemoveSingletonPPHMMs, PPHMMSorting):
    ''' Return all dirs for db storage and retrieval'''
    HMMERDir = f"{ExpDir}/HMMER"
    HMMER_PPHMMDir = f"{HMMERDir}/HMMER_PPHMMs"
    HMMER_PPHMMDbDir = f"{HMMERDir}/HMMER_PPHMMDb"
    HMMER_PPHMMDb = f"{HMMER_PPHMMDbDir}/HMMER_PPHMMDb"

    if RemoveSingletonPPHMMs == True:
        ClustersDir = f"{ExpDir}/BLAST/Clusters"
    else:
        ClustersDir = ""

    if PPHMMSorting == True:
        ClustersDir = f"{ExpDir}/BLAST/Clusters"
        HHsuiteDir = f"{ExpDir}/HHsuite"
        if os.path.exists(HHsuiteDir):
            shell(f"rm -rf {HHsuiteDir}")
            os.makedirs(HHsuiteDir)
        HHsuite_PPHMMDir = f"{HHsuiteDir}/HHsuite_PPHMMs"
        os.makedirs(HHsuite_PPHMMDir)
        HHsuite_PPHMMDBDir = f"{HHsuiteDir}/HHsuite_PPHMMDB"
        os.makedirs(HHsuite_PPHMMDBDir)
        HHsuite_PPHMMDB = f"{HHsuite_PPHMMDBDir}/HHsuite_PPHMMDB"
    else:
        HHsuiteDir, HHsuite_PPHMMDir, HHsuite_PPHMMDB = "", "", ""

    HMMER_hmmscanDir = HMMERDir+"/hmmscan_" + \
        ''.join(random.choice(string.ascii_uppercase + string.digits)
                for _ in range(10))
    os.makedirs(HMMER_hmmscanDir)

    return HMMER_hmmscanDir, ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDb, HHsuite_PPHMMDB, HHsuite_PPHMMDir, HHsuite_PPHMMDB, HHsuiteDir
