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

def retrieve_ucf_annots(fpath, IncludeIncompleteGenomes_UcfVirus, IncludeIncompleteGenomes_RefVirus) -> dict:
    if IncludeIncompleteGenomes_UcfVirus == True:
        if IncludeIncompleteGenomes_RefVirus == True:
            VariableShelveFile_UcfVirus = f"{fpath}/UcfVirusAnnotator.AllUcfGenomes.AllRefGenomes.p"

        elif IncludeIncompleteGenomes_RefVirus == False:
            VariableShelveFile_UcfVirus = f"{fpath}/UcfVirusAnnotator.AllUcfGenomes.CompleteRefGenomes.p"

    elif IncludeIncompleteGenomes_UcfVirus == False:
        if IncludeIncompleteGenomes_RefVirus == True:
            VariableShelveFile_UcfVirus = f"{fpath}/UcfVirusAnnotator.CompleteUcfGenomes.AllRefGenomes.p"

        elif IncludeIncompleteGenomes_RefVirus == False:
            VariableShelveFile_UcfVirus = f"{fpath}/UcfVirusAnnotator.CompleteUcfGenomes.CompleteRefGenomes.p"

    return pickle.load(open(f"{VariableShelveFile_UcfVirus}", "rb"))
