import pickle 

def retrieve_genome_vars(VariableShelveDir, IncludeIncompleteGenomes) -> dict:
    '''Read pickle file to dict containing VMR data'''
    if IncludeIncompleteGenomes == True:
        VariableShelveFile = f"{VariableShelveDir}/ReadGenomeDescTable.AllGenomes.p"
    elif IncludeIncompleteGenomes == False:
        VariableShelveFile = f"{VariableShelveDir}/ReadGenomeDescTable.CompleteGenomes.p"

    return pickle.load(open(VariableShelveFile, "rb"))

def retrieve_ref_virus_vars(VariableShelveDir, IncludeIncompleteGenomes) -> dict:
    '''Read pickle file to dict containing annotated ref viruses'''
    if IncludeIncompleteGenomes == True:
        VariableShelveFile = f"{VariableShelveDir}/RefVirusAnnotator.AllGenomes.p"
    elif IncludeIncompleteGenomes == False:
        VariableShelveFile = f"{VariableShelveDir}/RefVirusAnnotator.CompleteGenomes.p"

    return pickle.load(open(VariableShelveFile, "rb"))

def retrieve_pphmmdb_construction(VariableShelveDir) -> dict:
    '''Read pickle file to dict containing PPHMMDB'''
    return pickle.load(open(f"{VariableShelveDir}/PPHMMDBConstruction.p", "rb"))
