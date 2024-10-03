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
    '''HHSuite for optional sorting fns'''
    if PPHMMSorting == True:
        if os.path.exists(fnames['HHsuiteDir']):
            shell(f"rm -rf {fnames['HHsuiteDir']}")
            os.makedirs(fnames['HHsuiteDir'])
        os.makedirs(fnames['HHsuite_PPHMMDir'])
        os.makedirs(fnames['HHsuite_PPHMMDBDir'])
    else:
        pass

def mkdir_pl1_graphs(fnames, payload):
    '''Clean dirs for pl1 graph generator if already exists'''
    if payload['Bootstrap'] == True:
        '''Dendrogram distribution'''
        if os.path.isfile(fnames['VirusDendrogramDistFile']):
            os.remove(fnames['VirusDendrogramDistFile'])

        '''Bootstrapped Virus dendrogram'''
        if os.path.isfile(fnames['Heatmap_DendrogramFile']):
            os.remove(fnames['Heatmap_DendrogramFile'])

def mkdir_mi_scorer(fnames):
    '''Make MI Score dir'''
    if not os.path.exists(fnames['MutualInformationScoreDir']):
        os.makedirs(fnames['MutualInformationScoreDir'])

def mkdir_virus_classifier(fnames):
    if not os.path.exists(fnames['HMMER_PPHMMDB_UcfVirus']):
        raise ValueError(
            "Can't find HMMER PPHMM database of unclassified viruses")

    # if os.path.isfile(fnames['VirusDendrogramDistFile']):
    #     os.remove(fnames['VirusDendrogramDistFile'])

    '''Define path to bootstrapped dendrogram file'''
    if os.path.isfile(fnames['BootstrappedDendrogramFile']):
        os.remove(fnames['BootstrappedDendrogramFile'])
