import pickle 

def retrieve_variables(VariableShelveDir, IncludeIncompleteGenomes):
    '''Read pickle file containing VMR data'''
    if IncludeIncompleteGenomes == True:
        VariableShelveFile = VariableShelveDir + "/ReadGenomeDescTable.AllGenomes.p"
    elif IncludeIncompleteGenomes == False:
        VariableShelveFile = VariableShelveDir + "/ReadGenomeDescTable.CompleteGenomes.p"

    return pickle.load(open(VariableShelveFile, "rb"))
