import numpy as np

from app.utils.stdout_utils import clean_stdout
from app.utils.dcor import dcor

def SimilarityMat_Constructor(PPHMMSignatureTable, GOMSignatureTable, PPHMMLocationTable, SimilarityMeasurementScheme="PG", p=1.0):
    '''Construct similarity matrix according to specified scheme and p value'''
    if SimilarityMeasurementScheme not in ["P", "G", "L", "PG", "PL"]:
        return "'SimilarityMeasurementScheme' should be one of the following: 'P', 'G', 'L', 'PG', 'PL'."

    N_Viruses = PPHMMSignatureTable.shape[0]
    p = float(p)

    if "P" in SimilarityMeasurementScheme:
        PPHMMSignature_GJMat = np.zeros((N_Viruses, N_Viruses))
    if "G" in SimilarityMeasurementScheme:
        GOMSignature_GJMat = np.zeros((N_Viruses, N_Viruses))
    if "L" in SimilarityMeasurementScheme:
        PPHMMLocation_dCorMat = np.zeros((N_Viruses, N_Viruses))

    try:
        for i in range(N_Viruses):
            for j in range(i, N_Viruses):
                if "P" in SimilarityMeasurementScheme:
                    PPHMMSignature_i = PPHMMSignatureTable[i]
                    PPHMMSignature_j = PPHMMSignatureTable[j]
                    # TODO test ternary expression fixes zero div error. Break out to accessory fn? ("save divide")
                    PPHMMSignature_GJMat[i, j] = np.sum(np.minimum(
                        PPHMMSignature_i, PPHMMSignature_j))/np.sum(np.maximum(PPHMMSignature_i, PPHMMSignature_j)) if not np.sum(np.maximum(PPHMMSignature_i, PPHMMSignature_j)) == 0 else 0
                    PPHMMSignature_GJMat[j, i] = PPHMMSignature_GJMat[i, j]


                if "G" in SimilarityMeasurementScheme:
                    GOMSignature_i = GOMSignatureTable[i]
                    GOMSignature_j = GOMSignatureTable[j]
                    GOMSignature_GJMat[i, j] = np.sum(np.minimum( # TODO ditto save div
                        GOMSignature_i, GOMSignature_j))/np.sum(np.maximum(GOMSignature_i, GOMSignature_j))
                    GOMSignature_GJMat[j, i] = GOMSignature_GJMat[i, j]


                if "L" in SimilarityMeasurementScheme:
                    PPHMMLocation_i = PPHMMLocationTable[i]
                    PPHMMLocation_j = PPHMMLocationTable[j]
                    PresentPPHMM_IndexList = np.column_stack(
                        (PPHMMLocation_i != 0, PPHMMLocation_j != 0)).any(axis=1)
                    PPHMMLocation_dCorMat[i, j] = dcor(PPHMMLocation_i[PresentPPHMM_IndexList].reshape(
                        -1, 1), PPHMMLocation_j[PresentPPHMM_IndexList].reshape(-1, 1))
                    PPHMMLocation_dCorMat[j, i] = PPHMMLocation_dCorMat[i, j]

    except ZeroDivisionError as ex:
        raise SystemExit(f"{ex}\n"
                         f"This usually means that there is no similarity between any of the sequences in your database, so we can't generate any stats or figures.\n"
                         f"Check your input VMR: are there any entries, and do you expect them to be related? Check your input Genbank file: are there any sequences in it?")
    clean_stdout()

    if "P" in SimilarityMeasurementScheme:
        PPHMMSignature_GJMat[np.where(np.isnan(PPHMMSignature_GJMat))] = 0
        PPHMMSignature_GJMat[PPHMMSignature_GJMat < 0] = 0
    if "G" in SimilarityMeasurementScheme:
        GOMSignature_GJMat[np.where(np.isnan(GOMSignature_GJMat))] = 0
        GOMSignature_GJMat[GOMSignature_GJMat < 0] = 0
    if "L" in SimilarityMeasurementScheme:
        PPHMMLocation_dCorMat[np.where(np.isnan(PPHMMLocation_dCorMat))] = 0
        PPHMMLocation_dCorMat[PPHMMLocation_dCorMat < 0] = 0

    if SimilarityMeasurementScheme == "P":
        SimilarityMat = PPHMMSignature_GJMat
    elif SimilarityMeasurementScheme == "G":
        SimilarityMat = GOMSignature_GJMat
    elif SimilarityMeasurementScheme == "L":
        SimilarityMat = PPHMMLocation_dCorMat
    elif SimilarityMeasurementScheme == "PG":
        SimilarityMat = (PPHMMSignature_GJMat*GOMSignature_GJMat)**0.5
    elif SimilarityMeasurementScheme == "PL":
        SimilarityMat = PPHMMSignature_GJMat*PPHMMLocation_dCorMat
    else:
        return "'SimilarityMeasurementScheme' should be one of the following: 'P', 'G', 'L', 'PG', 'PL'."

    SimilarityMat[SimilarityMat < 0] = 0
    return SimilarityMat**p
