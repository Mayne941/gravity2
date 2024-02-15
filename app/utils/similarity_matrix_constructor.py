import numpy as np

from app.utils.stdout_utils import clean_stdout
from app.utils.dcor import dcor

def SimilarityMat_Constructor(PPHMMSignatureTable, GOMSignatureTable, PPHMMLocationTable, pphmm_neighbourhood_weight, pphmm_signature_score_threshold, SimilarityMeasurementScheme="PG", p=1.0, fnames=False):
    '''Construct similarity matrix according to specified scheme and p value'''
    if SimilarityMeasurementScheme not in ["P", "G", "L", "PG", "PL"]:
        return "'SimilarityMeasurementScheme' should be one of the following: 'P', 'G', 'L', 'PG', 'PL'."
    N_Viruses = PPHMMSignatureTable.shape[0]
    p = float(p)

    PPHMMSignature_GJMat = np.zeros((N_Viruses, N_Viruses))
    GOMSignature_GJMat = np.zeros((N_Viruses, N_Viruses))
    PPHMMLocation_dCorMat = np.zeros((N_Viruses, N_Viruses))

    import pandas as pd
    loc_df = pd.read_csv(fnames["PphmmLocs"], index_col=False)
    trim_df = loc_df.iloc[:,1:]
    trim_df = trim_df.apply(pd.to_numeric)
    '''Replace non-hits with NaN so as to not mess up mean calculations'''
    trim_df = trim_df.replace(0, np.nan)
    means = trim_df.mean().to_list()
    means_indices = np.argsort(means)

    # TODO SWITCH DISTANCE TO SIG TABLE AS PPHMM DIST MAX MAY != SIG MAX
    sig_df = pd.read_csv(fnames["PphmmAndGomSigs"], index_col=False)
    sig_trim_df = sig_df.iloc[:,2:]
    sig_trim_df = sig_trim_df.apply(pd.to_numeric)
    sig_trim_df = sig_trim_df.replace(0, np.nan)

    all_dists = []
    for row in np.array(sig_trim_df): # WAS trim_df
        dists = []
        for i in range(row.shape[0]):
            dists.append(round(np.abs(row[i] - means[i]), 4) if not row[i] == np.nan else 0)
        all_dists.append(dists)
    '''Rearrange matrix X to PPHMM loc order and Y to match GRAViTy heatmap (i.e. calculated tree)'''
    distance_arr = np.nan_to_num(np.array(all_dists).T[means_indices].T, 0)

    AMP_SINGLETONS = True
    neighbourhoods = {}
    for row_idx, row in enumerate(distance_arr):
        neighbourhood = 0
        for prof_idx, prof in enumerate(row):
            if prof_idx == 0:
                '''Skip first'''
                continue

            if prof != 0:
                '''If current profile is a hit...'''
                if AMP_SINGLETONS:
                    neighbourhood += 1
                else:
                    if row[prof_idx-1] != 0:
                        neighbourhood += 1
                if prof_idx == len(row) - 1:
                    '''If last entry and a hit'''
                    if not row_idx in neighbourhoods:
                        neighbourhoods[row_idx] = {}
                    # RM TODO < If neigh == prof_idx (i.e. all), row[prof_idx-neighbourhood-1] == -1 which breaks everything TODO CHECK LOGIC
                    # neighbourhoods[row_idx].update({prof_idx-neighbourhood: [neighbourhood+1, row[prof_idx-neighbourhood-1:prof_idx]]})
                    neighbourhoods[row_idx].update({prof_idx-neighbourhood: [neighbourhood+1, row[prof_idx-neighbourhood:prof_idx]]})

            else:
                '''If current not a hit...'''
                if neighbourhood == 0:
                    '''But no prior hit, continue'''
                    continue
                else:
                    '''If prior was a hit, sub dict: {start idx: [len, vals]}'''
                    if not row_idx in neighbourhoods:
                        neighbourhoods[row_idx] = {}
                    neighbourhoods[row_idx].update({prof_idx-neighbourhood: [neighbourhood+1, row[prof_idx-neighbourhood-1:prof_idx]]})
                    neighbourhood = 0

    MIN_NEIGH_SIZE = 5
    neigh_neg_weights = {}
    neigh_pos_weights = {}
    for row_key, row in neighbourhoods.items():
        for nei_key, nei in row.items():
            if nei[0] < MIN_NEIGH_SIZE:
                continue
            else:
                if row_key not in neigh_neg_weights.keys():
                    neigh_neg_weights[row_key] = []
                idxs = [i+nei_key for i in range(len(nei[1]))]
                '''Select largest parameter of neigh, select largest and non-largest'''
                if not row_key in neigh_pos_weights.keys():
                    neigh_pos_weights[row_key] = []
                neigh_pos_weights[row_key] = neigh_pos_weights[row_key] + [np.argmax(nei[1]) + nei_key]
                del idxs[np.argmax(nei[1])]                 # Delete largest number's idx
                neigh_neg_weights[row_key] = neigh_neg_weights[row_key] + idxs

    '''Reverse sort'''
    reversed = np.argsort(means_indices)
    unsort_neigh_neg_weights = {}
    for row_idx, row in neigh_neg_weights.items():
        '''Convert neg weight indices back to original'''
        unsort_neigh_neg_weights[row_idx] = []
        for prof in row:
            unsort_neigh_neg_weights[row_idx].append(np.where(reversed == prof)[0][0])
    unsort_neigh_pos_weights = {}
    for row_idx, row in neigh_pos_weights.items():
        '''Convert pos weight indices back to original'''
        unsort_neigh_pos_weights[row_idx] = []
        for prof in row:
            unsort_neigh_pos_weights[row_idx].append(np.where(reversed == prof)[0][0])

    try:
        for i in range(N_Viruses):
            for j in range(i, N_Viruses):
                if "P" in SimilarityMeasurementScheme:
                    PPHMMSignature_i = PPHMMSignatureTable[i]
                    PPHMMSignature_j = PPHMMSignatureTable[j]

                    '''Apply shared neighbourhood weights'''
                    NEG_WEIGHT = 1 - pphmm_neighbourhood_weight
                    POS_WEIGHT = 1 + pphmm_neighbourhood_weight

                    # NEG_WEIGHT = 1 - 0.95
                    # POS_WEIGHT = 1 + 0.2
                    try:
                        PPHMMSignature_i[unsort_neigh_neg_weights[i]] = PPHMMSignature_i[unsort_neigh_neg_weights[i]] * NEG_WEIGHT
                        PPHMMSignature_i[unsort_neigh_pos_weights[i]] = PPHMMSignature_i[unsort_neigh_pos_weights[i]] * POS_WEIGHT
                    except KeyError:
                        """No neighbourhood for this genome"""
                        ...
                    try:
                        PPHMMSignature_j[unsort_neigh_neg_weights[j]] = PPHMMSignature_j[unsort_neigh_neg_weights[j]] * NEG_WEIGHT
                        PPHMMSignature_j[unsort_neigh_pos_weights[j]] = PPHMMSignature_j[unsort_neigh_pos_weights[j]] * POS_WEIGHT
                    except KeyError:
                        """No neighbourhood for this genome"""
                        ...

                    '''THRESHOLD'''
                    PPHMMSignature_i[PPHMMSignature_i < pphmm_signature_score_threshold] = 0
                    PPHMMSignature_j[PPHMMSignature_j < pphmm_signature_score_threshold] = 0

                    # '''Apply absence weighting'''
                    # ABS_WEIGHT = 0.999
                    # n_matching = 0
                    # i_missing = 0
                    # j_missing = 0
                    # for idx, val in enumerate(PPHMMSignature_i):
                    #     if val != 0 and PPHMMSignature_j[idx] != 0:
                    #         n_matching += 1
                    #     elif val != 0 and PPHMMSignature_j[idx] == 0:
                    #         j_missing += 1
                    #     elif val == 0 and PPHMMSignature_j[idx] != 0:
                    #         i_missing += 1
                    # if sum(PPHMMSignature_i) > sum(PPHMMSignature_j):
                    #     small_arr, small_arr_missing = PPHMMSignature_j, j_missing
                    #     large_arr, large_arr_missing = PPHMMSignature_i, i_missing
                    # else:
                    #     small_arr, small_arr_missing = PPHMMSignature_i, i_missing
                    #     large_arr, large_arr_missing = PPHMMSignature_j, j_missing
                    # small_arr_hits, large_arr_hits = np.count_nonzero(small_arr), np.count_nonzero(large_arr)

                    # # PPHMMSignature_GJMat[i, j] = np.sum(small_arr) - ((small_arr_missing * np.median(small_arr[small_arr > 0])) * ABS_WEIGHT) / np.sum(
                    # #         large_arr) - ((large_arr_missing * np.median(large_arr[large_arr > 0])) * ABS_WEIGHT) if not np.sum(large_arr) == 0 else 0
                    # PPHMMSignature_GJMat[i, j] = np.sum(small_arr) * (1 - (small_arr_missing/small_arr_hits) * ABS_WEIGHT) / np.sum(
                    #         large_arr) * (1 - (large_arr_missing/large_arr_hits)* ABS_WEIGHT) if not np.sum(large_arr) == 0 else 0

                    PPHMMSignature_GJMat[i, j] = np.sum(np.minimum(
                        PPHMMSignature_i, PPHMMSignature_j))/np.sum(np.maximum(PPHMMSignature_i, PPHMMSignature_j)) if not np.sum(np.maximum(PPHMMSignature_i, PPHMMSignature_j)) == 0 else 0
                    PPHMMSignature_GJMat[j, i] = PPHMMSignature_GJMat[i, j]


                if "G" in SimilarityMeasurementScheme:
                    GOMSignature_i = GOMSignatureTable[i]
                    GOMSignature_j = GOMSignatureTable[j]
                    GOMSignature_GJMat[i, j] = np.sum(np.minimum(
                        GOMSignature_i, GOMSignature_j))/np.sum(np.maximum(GOMSignature_i, GOMSignature_j)) if not np.sum(np.maximum(GOMSignature_i, GOMSignature_j)) == 0 else 0
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

    sim_settings = {
            "P": PPHMMSignature_GJMat,
            "G": GOMSignature_GJMat,
            "L": PPHMMLocation_dCorMat
        }
    for k,v in sim_settings.items():
        '''Clear NaNs and negative numbers from GJ matrices for each scheme'''
        if k in SimilarityMeasurementScheme:
            v[np.where(np.isnan(v))] = 0
            v[v < 0] = 0

    if SimilarityMeasurementScheme == "P": # TODO map to a dict
        SimilarityMat = PPHMMSignature_GJMat
    elif SimilarityMeasurementScheme == "G":
        SimilarityMat = GOMSignature_GJMat
    elif SimilarityMeasurementScheme == "L":
        SimilarityMat = PPHMMLocation_dCorMat
    elif SimilarityMeasurementScheme == "PG":
        SimilarityMat = (PPHMMSignature_GJMat*GOMSignature_GJMat)**0.5
    elif SimilarityMeasurementScheme == "PL":
        SimilarityMat = PPHMMSignature_GJMat*PPHMMLocation_dCorMat

    SimilarityMat[SimilarityMat < 0] = 0

    return SimilarityMat**p
