import os
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

def mkdir_ref_annotator(fnames, PPHMMSorting):
    ''' Return all dirs for db storage and retrieval'''
    '''Hmmer '''
    os.makedirs(fnames['HMMER_hmmscanDir'])

    '''HHSuite for optional sorting fns'''
    if PPHMMSorting == True:
        if os.path.exists(fnames['HHsuiteDir']):
            shell(f"rm -rf {fnames['HHsuiteDir']}")
            os.makedirs(fnames['HHsuiteDir'])
        os.makedirs(fnames['HHsuite_PPHMMDir'])
        os.makedirs(fnames['HHsuite_PPHMMDBDir'])
    else:
        pass
