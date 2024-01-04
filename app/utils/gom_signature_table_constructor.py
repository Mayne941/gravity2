import numpy as np
from .dcor import dcor
from alive_progress import alive_it


def GOMSignatureTable_Constructor(PPHMMLocationTable, GOMDB, GOMIDList, bootstrap=0):
    '''Generate organisational model signature table, for annotations, graphing and description functions'''
    if bootstrap != 0:
        print(f"- (Re-)Constructing GOM Signature Table, bootstrap iteration: {bootstrap}")
    N_Viruses = len(PPHMMLocationTable)
    GOMSignatureTable = np.empty((N_Viruses, 0))
    for GOM in alive_it(GOMIDList):
        GOMSignatureList = []
        for PPHMMLocation in PPHMMLocationTable:
            RelevantPPHMMIndices = np.where(list(map(any, list(
                zip(list(map(any, GOMDB[GOM].transpose() != 0)), PPHMMLocation != 0)))))[0]
            GOMSignatureList.append(dcor(
                GOMDB[GOM][:, RelevantPPHMMIndices].T, PPHMMLocation[RelevantPPHMMIndices].reshape(-1, 1)))

        GOMSignatureTable = np.column_stack(
            (GOMSignatureTable, GOMSignatureList))

    return GOMSignatureTable
