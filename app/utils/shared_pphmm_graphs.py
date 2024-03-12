import plotly.express as px
from sklearn.metrics.pairwise import pairwise_distances
import pandas as pd
import numpy as np
from Bio import Phylo
import re
from copy import copy
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from app.utils.stdout_utils import progress_msg
from app.utils.heatmap_params import get_hmap_params, construct_hmap_lines, construct_wide_hmap_lines
from app.utils.make_heatmap_labels import make_labels, split_labels
from app.utils.error_handlers import raise_gravity_error

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

def build_discrete_hm(data, fnames, labels, map_type="null") -> None:
    '''Build and write heatmap'''
    if len(labels[0]) < 25:
        h, w = 1000, 1000
    else:
        h, w = 1500, 1500
    tit=""
    out_fname=f"{fnames['OutputDir']}/DISCRETE_CMAP_DISTANCES.png"
    fig = px.imshow(data, color_continuous_scale=["red","green","blue", "orange"], x=labels[0],y=labels[1],title=tit)
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

    # build_hm(shared_pphmms, fnames, labels, map_type="shared")
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

    # build_hm(shared_pphmms_ratio, fnames, labels, "norm")
    return shared_pphmms_ratio

def get_dist_data(fnames, label_order):
    loc_df = pd.read_csv(fnames["PphmmLocs"], index_col=False)
    trim_df = loc_df.iloc[:,1:]
    trim_df = trim_df.apply(pd.to_numeric)
    '''Replace non-hits with NaN so as to not mess up mean calculations'''
    trim_df = trim_df.replace(0, np.nan)
    # means = trim_df.mean().to_list()
    means = np.nanmedian(trim_df, axis=0)
    means_indices = np.argsort(means)

    all_dists = []
    for row in np.array(trim_df):
        dists = []
        for i in range(row.shape[0]):
            # dists.append(row[i] - means[i] if not row[i] == np.nan else np.nan) # nan for blank squares # TODO HOW IT WAS 02/02/24
            dists.append(row[i] if not row[i] == np.nan else np.nan) # nan for blank squares
        all_dists.append(dists)
    '''Rearrange matrix X to PPHMM loc order and Y to match GRAViTy heatmap (i.e. calculated tree)'''
    distance_arr = np.array(all_dists).T[means_indices].T[label_order]

    ROUNDING = 2
    cscale = np.round(np.array(means)[means_indices], ROUNDING)
    loc_dist_frm_mean_arr = np.copy(distance_arr)
    for row in loc_dist_frm_mean_arr:
        row_tot = np.nanmax(row) #new
        for idx, val in enumerate(row):
            if not np.isnan(val):
                row[idx] = round((cscale[idx] + val) / row_tot, ROUNDING)


    discrete_dist_arr = np.copy(distance_arr)
    for row in discrete_dist_arr:
        row_tot = np.nanmax(row)
        for idx, val in enumerate(row):
            if not np.isnan(val):
                # row[idx] = round(cscale[idx] + val, ROUNDING) # TODO HOW IT WAS 02/02/24
                row[idx] = round(val / row_tot, ROUNDING)

    return loc_dist_frm_mean_arr, means_indices, discrete_dist_arr

def pphmm_loc_distances(fnames, pphmm_names, label_order, labels):
    '''Reindex output df by mean position of PPHMM (NOT distance from mean)'''
    distance_arr, means_indices, discrete_dist_arr = get_dist_data(fnames, label_order)
    reindexed_xlabels = np.array(pphmm_names)[means_indices]
    # build_hm(distance_arr, fnames, [reindexed_xlabels, labels[1]])
    # build_discrete_hm(discrete_dist_arr, fnames, [reindexed_xlabels, labels[1]])

    return distance_arr, reindexed_xlabels, discrete_dist_arr

def pphmm_loc_diffs_pairwise(fnames, label_order, labels):
    '''Do pairwise distances of PPHMM location from mean'''
    distance_arr, _, _ = get_dist_data(fnames, label_order)
    out_arr = pairwise_distances(np.nan_to_num(distance_arr, nan=0), metric="braycurtis")
    # build_hm(out_arr, fnames, labels, map_type="loc")
    return out_arr

def supplementary_pphmm_heatmaps(pl1_ref_annotations, TaxoLabelList_AllVirus, N_RefViruses, TaxoLabelList_RefVirus,
                                TaxoAssignmentList, data, out_fname, fnames, payload, final_results, is_square=True, pl2_components=False):
    '''Construct GRAViTy heat map with dendrogram'''
    progress_msg(f"Constructing supplementary GRAViTy Heatmap: {out_fname}")
    '''Load the tree'''
    if payload["Bootstrap"]:
        VirusDendrogramFile = fnames['BootstrappedDendrogramFile']
    else:
        VirusDendrogramFile = fnames['Heatmap_DendrogramFile']

    VirusDendrogram = Phylo.read(VirusDendrogramFile, "newick")

    '''Determine virus order'''
    _ = VirusDendrogram.ladderize(reverse=True)
    OrderedTaxoLabelList = [
        Clade.name for Clade in VirusDendrogram.get_terminals()]
    VirusOrder = [TaxoLabelList_AllVirus.index(
        TaxoLabel) for TaxoLabel in OrderedTaxoLabelList]

    '''Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff'''
    N_InternalNodes = len(VirusDendrogram.get_nonterminals())
    for InternalNode_i in range(N_InternalNodes):
        try:
            '''Here we want to ensure there are no labels on the dendrogram where bootstrap confidence < cutoff (or not present if BS = False)'''
            if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence < payload['Heatmap_DendrogramSupport_Cutoff'] or np.isnan(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence):
                continue
            else:
                VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(
                    VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
        except:
            '''Exception as first entry will always be None'''
            continue

    if pl2_components:
        '''PL2: Colour terminal branches: reference virus's branch is blue, unclassified virus's branch is red'''
        N_Viruses = N_RefViruses + pl2_components[0]
        for Virus_i in range(N_Viruses):
            Taxolabel = VirusDendrogram.get_terminals()[Virus_i].name
            if Taxolabel in TaxoLabelList_RefVirus:
                VirusDendrogram.get_terminals()[Virus_i].color = "blue"
            elif Taxolabel in pl2_components[1]:
                VirusDendrogram.get_terminals()[Virus_i].color = "red"
    else:
        '''PL1'''
        N_Viruses = N_RefViruses
        for Virus_i in range(N_Viruses):
            Taxolabel = VirusDendrogram.get_terminals()[Virus_i].name

    '''Labels, label positions, and ticks'''
    ClassDendrogram = copy(VirusDendrogram)
    ClassDendrogram_label = Phylo.read(VirusDendrogramFile, "newick")
    TaxoGroupingList_AllVirus = pl1_ref_annotations["TaxoGroupingList"].tolist(
    ) + TaxoAssignmentList

    '''Construct taxo labels in ordered list for heatmap axis labels'''
    _, LineList_major = make_labels(ClassDendrogram, zip(
        TaxoLabelList_AllVirus, TaxoGroupingList_AllVirus))

    _, LineList_minor = make_labels(ClassDendrogram_label, zip(
        TaxoLabelList_AllVirus, OrderedTaxoLabelList))
    ClassLabelList_x, ClassLabelList_y = split_labels(OrderedTaxoLabelList)

    '''Plot configuration'''
    if not is_square:
        n_pphmms = data.shape[1]
        cmap='rainbow'
    else:
        n_pphmms = 99
        cmap="magma"
    hmap_params, fig, ax_Dendrogram, ax_Heatmap = get_hmap_params(len(ClassLabelList_x), n_pphmms=n_pphmms, is_square=is_square)

    try:
        Phylo.draw(VirusDendrogram, label_func=lambda x: "",
                do_show=False,  axes=ax_Dendrogram)
    except Exception as ex:
        raise raise_gravity_error(f"ERROR: {ex}\nThis will usually occur when there's been an error with bootstrapping. Try swapping bootstrap method or disabling, then try again.")

    VirusDendrogramDepth = max(
        [v for k, v in VirusDendrogram.depths().items()])
    ax_Dendrogram		.set_xlim(
        [(VirusDendrogramDepth - 1), VirusDendrogramDepth])
    ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
    ax_Dendrogram		.set_axis_off()

    '''Draw heatmap elements on axes'''
    vmax = np.max(np.nan_to_num(data, nan=0))
    Heatmap_Graphic = ax_Heatmap.imshow(
        data, cmap=cmap, aspect='auto', vmin=0, vmax=vmax, interpolation='nearest', filternorm=False) ## TODO test for slow speed

    '''Draw grouping major & minor lines'''
    if is_square:
        ax_Heatmap = construct_hmap_lines(ax_Heatmap, len(ClassLabelList_y), LineList_major, LineList_minor,
                                        hmap_params, ClassLabelList_x, ClassLabelList_y,
                                        TickLocList = np.array(
                                            list(map(np.mean, list(zip(LineList_minor[0:-1], LineList_minor[1:])))))
                                            )
    else:
        ax_Heatmap = construct_wide_hmap_lines(ax_Heatmap, data.shape[1], LineList_major, LineList_minor,
                                        hmap_params, final_results["PphmLabels"], ClassLabelList_y, np.array(
                                            list(map(np.mean, list(zip(LineList_minor[0:-1], LineList_minor[1:]))))),
                                            data
                                            )
        # n = 500 # TODO test, custom labels for v large dna graphs
        # [l.set_visible(False) for (i,l) in enumerate(ax_Heatmap.xaxis.get_ticklabels()) if i % n != 0]

    '''Selectively colour tick labels red if a UCF sample'''
    [i.set_color("red") for i in ax_Heatmap.get_xticklabels()
    if bool(re.match(r"Query", i.get_text()))]
    [i.set_color("red") for i in ax_Heatmap.get_yticklabels()
    if bool(re.match(r"Query", i.get_text()))]

    '''Reference virus colour bar'''
    ax_CBar = fig.add_axes(
        [hmap_params['ax_CBar_L'], hmap_params['ax_CBar_B'], hmap_params['ax_CBar_W'], hmap_params['ax_CBar_H']], frame_on=True, facecolor="white")
    CBar_Graphic = fig.colorbar(
        Heatmap_Graphic, cax=ax_CBar, orientation="horizontal", ticks=[round(i, -1) for i in np.linspace(0, vmax, 5)]) # DYNAMIC TICKS
    CBar_Graphic		.ax.set_xticklabels(
        ['0', '0.25', '0.50', '0.75', '1'], rotation=0, size=hmap_params['FontSize'])
    CBar_Graphic		.ax.set_xlabel('Distance', rotation=0, size=hmap_params['FontSize']+2)
    CBar_Graphic		.ax.tick_params(top=False,
                                bottom=True,
                                left=False,
                                right=False,
                                labeltop=False,
                                labelbottom=True,
                                labelleft=False,
                                labelright=False,
                                direction='out')

    plt.savefig(f"{out_fname}", format="pdf", bbox_inches = "tight", dpi=hmap_params['dpi']) ## FNAME
    plt.cla()
    plt.clf()
    plt.close()
    return VirusOrder, ClassLabelList_x
