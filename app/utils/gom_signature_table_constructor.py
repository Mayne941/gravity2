# import numpy as np
# from .dcor import dcor
# from alive_progress import alive_it


# def GOMSignatureTable_Constructor(PPHMMLocationTable, GOMDB, GOMIDList, bootstrap=0):
#     '''Generate organisational model signature table, for annotations, graphing and description functions'''
#     if bootstrap != 0:
#         print(f"- (Re-)Constructing GOM Signature Table, bootstrap iteration: {bootstrap+1}")
#     N_Viruses = len(PPHMMLocationTable)
#     GOMSignatureTable = np.empty((N_Viruses, 0))
#     for GOM in alive_it(GOMIDList): ###############
#         GOMSignatureList = []
#         for PPHMMLocation in PPHMMLocationTable:
#             RelevantPPHMMIndices = np.where(list(map(any, list(
#                 zip(list(map(any, GOMDB[GOM].transpose() != 0)), PPHMMLocation != 0)))))[0]
#             GOMSignatureList.append(dcor(
#                 GOMDB[GOM][:, RelevantPPHMMIndices].T, PPHMMLocation[RelevantPPHMMIndices].reshape(-1, 1)))

#         GOMSignatureTable = np.column_stack(
#             (GOMSignatureTable, GOMSignatureList))

#     return GOMSignatureTable

import os
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
from alive_progress import alive_it

from app.utils.dcor import dcor
from app.utils.stdout_utils import progress_msg

def GOMSignatureTable_Constructor(PPHMMLocationTable, GOMDB, GOMIDList, bootstrap=0):
    '''Generate organisational model signature table, for annotations, graphing and description functions'''
    if bootstrap != 0:
        print(f"- (Re-)Constructing GOM Signature Table, bootstrap iteration: {bootstrap+1}")
    N_Viruses = len(PPHMMLocationTable)
    GOMSignatureTable = np.empty((N_Viruses, 0))

    progress_msg(f"-  Spinning up {os.cpu_count()-1} workers to generate GOM signatures. This may take a while...")
    pool = Pool(os.cpu_count()-1)
    with pool as p, tqdm(total=len(GOMIDList)) as pbar:
        res = [p.apply_async(
            generate_gom_sigs, args=(i,PPHMMLocationTable, GOMDB), callback=lambda _: pbar.update(1)) for i in GOMIDList]
        results = [r.get() for r in res]

    for result in results:
        GOMSignatureTable = np.column_stack(
            (GOMSignatureTable, result))

    return GOMSignatureTable

def generate_gom_sigs(GOM, PPHMMLocationTable, GOMDB):
    GOMSignatureList = []
    for PPHMMLocation in PPHMMLocationTable:
        RelevantPPHMMIndices = np.where(list(map(any, list(
            zip(list(map(any, GOMDB[GOM].transpose() != 0)), PPHMMLocation != 0)))))[0]
        GOMSignatureList.append(dcor(
            GOMDB[GOM][:, RelevantPPHMMIndices].T, PPHMMLocation[RelevantPPHMMIndices].reshape(-1, 1)))
    return GOMSignatureList
