import plotly.express as px
from sklearn.metrics.pairwise import pairwise_distances
import pandas as pd
import numpy as np

def get_sig_data(fnames, label_order) -> np.array:
    '''Load PPHMM/GOM sigs from csv, remove non-numerical cols, cast to array'''
    df = pd.read_csv(fnames['PphmmAndGomSigs'], index_col=0)
    df = df.iloc[:,1:]
    df = df.apply(pd.to_numeric)
    return np.asarray(df.reindex(label_order))

def build_hm(data, fnames, labels, map_type="null") -> None:
    '''Build and write heatmap'''
    if len(labels[0]) < 25:
        h, w = 1000, 1000
    else:
        h, w = 1500, 1500
    if map_type == "shared":
        tit = "Shared N PPHMMs"
        out_fname = fnames['SharedPphmmHm']
    elif map_type == "norm":
        tit = "Shared Normalised PPHMM Ratio (pairwise: n common PPHMMs / n PPHMMs)"
        out_fname = fnames["NormPphmmRatioHm"]
    elif map_type == "loc":
        tit = "Average normalised absolute distance (nt) between shared PPHMMs, Bray-Curtis distance (sum(|U_i - V_i|)/sum(|U_i + V_i|)) "
        out_fname = fnames["PphmmLocPairwise"]
    else:
        tit = "Distance of PPHMMs (nt) from their mean position, x=PPHMMs, y=genomes. PPHMMs sorted by estimated genome position (<-5', 3'->)"
        out_fname = fnames["PphmmLocDists"]
    fig = px.imshow(data, x=labels[0],y=labels[1],title=tit)
    fig.layout.height = h
    fig.layout.width  = w if not map_type == "null" else w*1.5
    fig.write_image(out_fname)

def shared_pphmm_ratio(label_order, fnames, labels) -> None:
    '''Build matrix of PPHMM co-occurrence, plot heatmap'''
    arr = get_sig_data(fnames, label_order)

    shared_pphmms = []
    for row in arr:
        row_res = []
        for i in range(arr.shape[0]):
            cur_row = row
            com_row = arr[i]
            n_hits  = 0
            for j in range(row.shape[0]):
                if cur_row[j] != 0 and com_row[j] != 0:
                    n_hits += 1
            row_res.append(n_hits)
        shared_pphmms.append(row_res)

    build_hm(shared_pphmms, fnames, labels, map_type="shared")
    return shared_pphmms

def shared_norm_pphmm_ratio(label_order, fnames, labels) -> None:
    '''Build matrix of normalised ratio of shared PPHMMs, plot heatmap'''
    arr = get_sig_data(fnames, label_order)

    def _f(a,b, max_possible):
        '''Count shared PPHMM n, return normalised save division'''
        hits = len([a[idx] for idx in range(a.shape[0]) if not a[idx] == 0 and not b[idx] == 0])
        return min([hits, max_possible]) / max([hits,max_possible])

    '''Construct matrix'''
    l = arr.shape[0]
    shared_pphmms_ratio = np.zeros((l,l))
    for idx_a, a in enumerate(arr):
        for idx_b, b in enumerate(arr):
            max_possible = min([np.count_nonzero(a), np.count_nonzero(b)])
            shared_pphmms_ratio[idx_b,idx_a] = _f(a,b, max_possible)

    build_hm(shared_pphmms_ratio, fnames, labels, "norm")
    return shared_pphmms_ratio

def get_dist_data(fnames, label_order):
    loc_df = pd.read_csv(fnames["PphmmLocs"], index_col=False)
    trim_df = loc_df.iloc[:,1:]
    trim_df = trim_df.apply(pd.to_numeric)
    '''Replace non-hits with NaN so as to not mess up mean calculations'''
    trim_df = trim_df.replace(0, np.nan)
    means = trim_df.mean().to_list()
    means_indices = np.argsort(means) ## HOW IS THIS ORIGINAL GRAV ORDER? THIS IS JUST X AXIS

    all_dists = []
    for row in np.array(trim_df):
        dists = []
        for i in range(row.shape[0]):
            dists.append(round(np.abs(row[i] - means[i]), 4) if not row[i] == np.nan else np.nan) # nan for blank squares
        all_dists.append(dists)
    '''Rearrange matrix X to PPHMM loc order and Y to match GRAViTy heatmap (i.e. calculated tree)'''
    # distance_arr = np.array(all_dists).T[means_indices].T
    distance_arr = np.array(all_dists).T[means_indices].T[label_order]
    return distance_arr, means_indices

def pphmm_loc_distances(fnames, pphmm_names, label_order, labels):
    '''Reindex output df by mean position of PPHMM (NOT distance from mean)'''
    distance_arr, means_indices = get_dist_data(fnames, label_order)
    reindexed_xlabels = np.array(pphmm_names)[means_indices]
    build_hm(distance_arr, fnames, [reindexed_xlabels, labels[1]])
    return distance_arr, reindexed_xlabels

def pphmm_loc_diffs_pairwise(fnames, label_order, labels):
    '''Do pairwise distances of PPHMM location from mean'''
    distance_arr, _ = get_dist_data(fnames, label_order)
    out_arr = pairwise_distances(np.nan_to_num(distance_arr, nan=0), metric="braycurtis")
    build_hm(out_arr, fnames, labels, map_type="loc")
    return out_arr
