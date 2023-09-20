import pickle

from app.utils.error_handlers import raise_gravity_error

def retrieve_genome_vars(VariableShelveDir, IncludeIncompleteGenomes) -> dict:
    '''Read pickle file to dict containing VMR data'''
    if IncludeIncompleteGenomes == True:
        VariableShelveFile = f"{VariableShelveDir}/ReadGenomeDescTable.AllGenomes.p"
    elif IncludeIncompleteGenomes == False:
        VariableShelveFile = f"{VariableShelveDir}/ReadGenomeDescTable.CompleteGenomes.p"

    f = pickle.load(open(VariableShelveFile, "rb"))
    if len(f["VirusNameList"]) == 0:
        raise_gravity_error("The list of processed genomes I retrieved is empty!\nThis usually results from errors reading the VMR or genomes in your VMR being incomplete (try enabling AnnotateIncompleteGenomes).")

    return f

def retrieve_ref_virus_vars(VariableShelveDir, IncludeIncompleteGenomes) -> dict:
    '''Read pickle file to dict containing annotated ref viruses'''
    if IncludeIncompleteGenomes == True:
        VariableShelveFile = f"{VariableShelveDir}/RefVirusAnnotator.AllGenomes.p"
    elif IncludeIncompleteGenomes == False:
        VariableShelveFile = f"{VariableShelveDir}/RefVirusAnnotator.CompleteGenomes.p"

    f = pickle.load(open(VariableShelveFile, "rb"))
    if f["ClusterIDList"].shape[0] == 0:
        raise_gravity_error("The list of processed genomes I retrieved is empty!\nThis usually results from errors reading the VMR or genomes in your VMR being incomplete (try enabling AnnotateIncompleteGenomes).")

    return f

def retrieve_pphmmdb_construction(VariableShelveDir) -> dict:
    '''Read pickle file to dict containing PPHMMDB''' # RM < TODO Add import test
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

    f = pickle.load(open(VariableShelveFile_UcfVirus, "rb"))
    if f["ClusterIDList"].shape[0] == 0:
        raise_gravity_error("The list of processed genomes I retrieved is empty!\nThis usually results from errors reading the VMR or genomes in your VMR being incomplete (try enabling AnnotateIncompleteGenomes).")


    return f
