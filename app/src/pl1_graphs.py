from app.utils.make_heatmap_labels import make_labels, split_labels
from app.utils.shell_cmds import shell
from app.utils.dist_mat_to_tree import DistMat2Tree
from app.utils.gomdb_constructor import GOMDB_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.taxo_label_constructor import TaxoLabel_Constructor
from app.utils.similarity_matrix_constructor import SimilarityMat_Constructor
from app.utils.virus_grouping_estimator import VirusGrouping_Estimator
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_pickle
from app.utils.stdout_utils import warning_msg, progress_msg
from app.utils.generate_fnames import generate_file_names
from app.utils.mkdirs import mkdir_pl1_graphs
from app.utils.heatmap_params import get_hmap_params, construct_hmap_lines
from app.utils.error_handlers import raise_gravity_error
from app.utils.shared_pphmm_graphs import supplementary_pphmm_heatmaps, shared_norm_pphmm_ratio, shared_pphmm_ratio, pphmm_loc_distances, pphmm_loc_diffs_pairwise

import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
import matplotlib
matplotlib.use('agg')


class GRAViTyDendrogramAndHeatmapConstruction:
    def __init__(self,
                 payload,
                 ExpDir,
                 ) -> None:
        '''Parameters'''
        self.payload = payload
        self.fnames = generate_file_names(payload, ExpDir)
        self.genomes = retrieve_genome_vars(self.fnames['ReadGenomeDescTablePickle'])
        self.ref_annotations = retrieve_pickle(self.fnames['RefAnnotatorPickle'])
        mkdir_pl1_graphs(self.fnames, self.payload)

    def dist_mat_constructor(self, PPHMMSignatureTable, GOMSignatureTable, PPHMMLocationTable, SPRSignatureTable):
        '''Estimate pairwise distances between viruses, return distance matrix (recip. sim matrix)'''
        SimMat = SimilarityMat_Constructor(PPHMMSignatureTable,
                                           GOMSignatureTable,
                                           PPHMMLocationTable,
                                           SPRSignatureTable,
                                           self.payload["PphmmNeighbourhoodWeight"],
                                           self.payload["PphmmSigScoreThreshold"],
                                           self.payload['SimilarityMeasurementScheme'],
                                           self.payload['p'],
                                           self.fnames
                                           )
        DistMat = 1 - SimMat
        DistMat[DistMat < 0] = 0
        return DistMat

    def dendrogram(self, DistMat, interim_Rscheme_matrix):
        '''Build dendrogram'''
        progress_msg("Constructing Dendrogram")
        '''Make TaxoLabelList'''
        TaxoLabelList = TaxoLabel_Constructor(SeqIDLists=self.genomes["SeqIDLists"],
                                              FamilyList=self.genomes["FamilyList"],
                                              GenusList=self.genomes["GenusList"],
                                              VirusNameList=self.genomes["VirusNameList"]
                                              )
        '''Make dendrogram'''
        VirusDendrogram = DistMat2Tree(DistMat=DistMat,
                                       LeafList=TaxoLabelList,
                                       Dendrogram_LinkageMethod=self.payload['Dendrogram_LinkageMethod'],
                                       do_logscale=self.payload["LogscaleDendrogram"]
                                       )
        with open(self.fnames['Heatmap_DendrogramFile'], "w") as VirusDendrogram_txt:
            VirusDendrogram_txt.write(VirusDendrogram)

        if self.payload['Bootstrap'] == True:
            '''Compute bootstrap support'''
            progress_msg("Computing Dendrogram Bootstrap Support")
            N_PPHMMs = self.ref_annotations["PPHMMSignatureTable"].shape[1]
            for Bootstrap_i in range(0, self.payload['N_Bootstrap']):
                print(f"Round {Bootstrap_i+1}/{self.payload['N_Bootstrap']}")
                '''Construct bootstrapped PPHMMSignatureTable and PPHMMLocationTable'''
                PPHMM_IndexList = np.random.choice(
                    list(range(N_PPHMMs)), N_PPHMMs, replace=True)
                BootstrappedPPHMMSignatureTable = self.ref_annotations[
                    "PPHMMSignatureTable"][:, PPHMM_IndexList]
                BootstrappedPPHMMLocationTable = self.ref_annotations[
                    "PPHMMSignatureTable"][:, PPHMM_IndexList]
                BootstrappedGOMSignatureTable = None

                if "G" in self.payload['SimilarityMeasurementScheme']:
                    '''Construct bootstrapped GOMSignatureTable'''
                    BootstrappedGOMDB = GOMDB_Constructor(TaxoGroupingList=self.genomes["TaxoGroupingList"],
                                                          PPHMMLocationTable=BootstrappedPPHMMLocationTable,
                                                          GOMIDList=self.ref_annotations["GOMIDList"]
                                                          )
                    BootstrappedGOMSignatureTable = GOMSignatureTable_Constructor(PPHMMLocationTable=BootstrappedPPHMMLocationTable,
                                                                                  GOMDB=BootstrappedGOMDB,
                                                                                  GOMIDList=self.ref_annotations["GOMIDList"],
                                                                                  bootstrap=Bootstrap_i
                                                                                  )

                '''Construct a dendrogram from the bootstrapped data'''
                BootstrappedDistMat = self.dist_mat_constructor(BootstrappedPPHMMSignatureTable,
                                                               BootstrappedGOMSignatureTable,
                                                               BootstrappedPPHMMLocationTable,
                                                               interim_Rscheme_matrix
                                                               )

                '''Generate bs dist matrix tree'''
                BootstrappedVirusDendrogram = DistMat2Tree(DistMat=BootstrappedDistMat,
                                                           LeafList=TaxoLabelList,
                                                           Dendrogram_LinkageMethod=self.payload['Dendrogram_LinkageMethod'],
                                                           do_logscale=self.payload["LogscaleDendrogram"]
                                                           )
                with open(self.fnames['VirusDendrogramDistFile'], "a") as VirusDendrogramDist_txt:
                    VirusDendrogramDist_txt.write(
                        BootstrappedVirusDendrogram+"\n")

            if os.path.isfile(self.fnames['BootstrappedDendrogramFile']):
                '''Bootstrap gets upset if trying to overwrite existing file'''
                os.remove(self.fnames['BootstrappedDendrogramFile'])

            '''Create bootstrapped dendrogram. N.b. No error handler as successful output goes to stderr'''
            if self.payload['Bootstrap_method'] == "booster":
                shell(f"./booster_linux64 -i {self.fnames['Heatmap_DendrogramFile']} -b {self.fnames['VirusDendrogramDistFile']} -o {self.fnames['BootstrappedDendrogramFile']} -@ {self.payload['N_CPUs']}")

            elif self.payload['Bootstrap_method'] == "sumtrees":
                shell(f"sumtrees.py --decimals=2 --no-annotations --preserve-underscores --force-rooted --output-tree-format=newick --output-tree-filepath={self.fnames['BootstrappedDendrogramFile']} --target={self.fnames['Heatmap_DendrogramFile']} {self.fnames['VirusDendrogramDistFile']}")
            else:
                warning_msg(
                    f"WARNING: 'Bootstrap_method' can either be 'booster' or 'sumtrees', but you gave me {self.payload['Bootstrap_method']}.")

    def heatmap_with_dendro(self, DistMat):
        '''6/7: (OPT) Construct GRAViTy heat map with dendrogram, save as pdf.
        Needs to be independent of stage 4 as user may provide pre-computer dendrogram file'''
        progress_msg("Generating GRAViTy Heatmap")
        N_Viruses = len(DistMat)

        '''Load tree'''
        if self.payload["Bootstrap"]:
            VirusDendrogram = Phylo.read(self.fnames['BootstrappedDendrogramFile'], "newick")
        else:
            VirusDendrogram = Phylo.read(self.fnames['Heatmap_DendrogramFile'], "newick")
        '''Determine virus order'''
        _ = VirusDendrogram.ladderize(reverse=True)

        '''Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff'''
        N_InternalNodes = len(VirusDendrogram.get_nonterminals())
        for InternalNode_i in range(N_InternalNodes):
            try:
                if np.isnan(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence):
                    '''Some bootstrap programs output NaN instead of zero, which breaks the parser.'''
                    VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = None

                if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence < self.payload['Heatmap_DendrogramSupport_Cutoff']:
                    VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = None

                else:
                    VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(
                        VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
            except:
                continue

        TaxoLabelList = TaxoLabel_Constructor(SeqIDLists=self.genomes["SeqIDLists"],
                                              FamilyList=self.genomes["FamilyList"],
                                              GenusList=self.genomes["GenusList"],
                                              VirusNameList=self.genomes["VirusNameList"]
                                              )

        for Virus_i in range(N_Viruses):
            Taxolabel = VirusDendrogram.get_terminals()[Virus_i].name
            if not "Query" in str(Taxolabel):
                VirusDendrogram.get_terminals()[Virus_i].color = "blue"
            else:
                VirusDendrogram.get_terminals()[Virus_i].color = "red"

        OrderedTaxoLabelList = [
            Clade.name for Clade in VirusDendrogram.get_terminals()]
        VirusOrder = [TaxoLabelList.index(TaxoLabel)
                      for TaxoLabel in OrderedTaxoLabelList]

        '''Re-order the distance matrix'''
        OrderedDistMat = DistMat[VirusOrder][:, VirusOrder]

        '''Labels, label positions, and ticks'''
        ClassDendrogram_grp = VirusDendrogram
        if self.payload["Bootstrap"]:
            ClassDendrogram_label = Phylo.read(  # Needs to be separately read in as Bio entities are linked
                self.fnames['BootstrappedDendrogramFile'], "newick")
        else:
            ClassDendrogram_label = Phylo.read(self.fnames['Heatmap_DendrogramFile'], "newick")

        '''Construct taxo labels in ordered list for heatmap axis labels'''
        _, LineList_major = make_labels(ClassDendrogram_grp, zip(
            TaxoLabelList, self.genomes["TaxoGroupingList"]))
        _, LineList_minor = make_labels(ClassDendrogram_label, zip(
            TaxoLabelList, self.genomes["VirusNameList"]))
        ClassLabelList_x, ClassLabelList_y = split_labels(OrderedTaxoLabelList)

        '''Plot configuration'''
        hmap_params, fig, ax_Dendrogram, ax_Heatmap = get_hmap_params(len(ClassLabelList_x),virus_dendrogram=VirusDendrogram, do_logscale=self.payload["LogscaleDendrogram"])

        try:
            Phylo.draw(VirusDendrogram, label_func=lambda x: "",
                    do_show=False,  axes=ax_Dendrogram)
            plt.rcParams['lines.linewidth'] = 1.0
        except Exception as ex:
            raise raise_gravity_error(f"ERROR: {ex}\nThis will usually occur when there's been an error with bootstrapping. Try disabling this feature and trying again.")
        dendro_x_min = 1

        ## ??
        # VirusDendrogramDepth = max(
        #         [v for k, v in VirusDendrogram.depths().items()])
        # if self.payload["LogscaleDendrogram"]:
        #     VirusDendrogramDepth = np.log10(VirusDendrogramDepth)
        #     dendro_x_min = 10**dendro_x_min
        # ax_Dendrogram		.set_xlim(
        #     [(VirusDendrogramDepth - dendro_x_min), VirusDendrogramDepth])
        # ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
        # ax_Dendrogram		.set_axis_off()
        ## ??

        if self.payload["LogscaleDendrogram"]:
            dendro_x_min = 10**dendro_x_min
        VirusDendrogramDepth = max(
            [v for k, v in VirusDendrogram.depths().items()])
        ax_Dendrogram		.set_xlim(
            [(VirusDendrogramDepth-dendro_x_min),VirusDendrogramDepth])
        ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
        ax_Dendrogram		.set_axis_off()

        '''Draw heatmap elements on axes'''
        Heatmap_Graphic = ax_Heatmap.imshow(
            OrderedDistMat, cmap='magma', aspect='auto', vmin=0, vmax=1, interpolation='none')

        '''Draw grouping major & minor lines'''
        ax_Heatmap = construct_hmap_lines(ax_Heatmap, len(ClassLabelList_y), LineList_major, LineList_minor,
                                          hmap_params, ClassLabelList_x, ClassLabelList_y,
                                          TickLocList = np.arange(0, len(ClassLabelList_y)))

        [i.set_color("red") for i in ax_Heatmap.get_xticklabels()
         if bool(re.search(r"Query", i.get_text().replace(" ","")))]
        [i.set_color("red") for i in ax_Heatmap.get_yticklabels()
         if bool(re.search(r"Query", i.get_text()))]

        '''Heatmap colourbars'''
        ax_CBar = fig.add_axes(
            [hmap_params['ax_CBar_L'], hmap_params['ax_CBar_B'], hmap_params['ax_CBar_W'], hmap_params['ax_CBar_H']], frame_on=True, facecolor="white")
        CBar_Graphic = fig.colorbar(
            Heatmap_Graphic, cax=ax_CBar, orientation="horizontal", ticks=[0, 0.25, 0.50, 0.75, 1])
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

        '''Save fig'''
        plt.savefig(self.fnames['HeatmapWithDendrogramFile'], format="pdf", bbox_inches = "tight", dpi=hmap_params['dpi'])
        plt.clf()
        plt.cla()
        return VirusOrder, TaxoLabelList, ClassLabelList_x

    def virus_grouping(self, DistMat):
        '''Group viruses via Thiels-U and other metrics; save as txt.'''
        progress_msg("Calculating virus groupings")
        VirusGroupingList, OptDistance_Cutoff, CorrelationScore, Theils_u_TaxoGroupingListGivenPred, \
            Theils_u_PredGivenTaxoGroupingList = VirusGrouping_Estimator(
                DistMat, self.payload['Dendrogram_LinkageMethod'], self.genomes["TaxoGroupingList"])
        np.savetxt(fname=self.fnames['VirusGroupingFile'],
                   X=np.column_stack((list(map(", ".join, self.genomes["SeqIDLists"])),
                                      self.genomes["FamilyList"],
                                      self.genomes["GenusList"],
                                      self.genomes["VirusNameList"],
                                      self.genomes["TaxoGroupingList"],
                                      VirusGroupingList,
                                      )),
                   fmt='%s',
                   delimiter="\t",
                   header="Sequence identifier\tFamily\tGenus\tVirus name\tClass\tGrouping")

        with open(self.fnames['VirusGroupingFile'], "a") as VirusGrouping_txt:
            VirusGrouping_txt.write(
                f"\n"
                f"Distance cut off: {OptDistance_Cutoff}\n"
                f"Theil's uncertainty correlation for the reference assignments given the predicted grouping U(Ref|Pred): {Theils_u_TaxoGroupingListGivenPred}\n"
                f"Theil's uncertainty correlation for the predicted grouping given the reference assignments U(Pred|Ref): {Theils_u_PredGivenTaxoGroupingList}\n"
                f"Symmetrical Theil's uncertainty correlation between the reference assignments and the predicted grouping U(Ref, Pred): {CorrelationScore}\n"
                f"U(X|Y) == 1 means that knowing Y implies a perfect knowledge of X, but not vice-versa\n"
                f"U(X,Y) == 1 means that knowing Y implies a perfect knowledge of X and vice-versa\n"
            )

    def create_supplementary_pphmm_tables(self):
        '''Combine PPHMM/GOM sigs, save as CSV for shared PPHMM heatmaps'''
        pl1_ref_annotations = {**self.genomes, **self.ref_annotations}
        TaxoLabelList = TaxoLabel_Constructor(SeqIDLists=self.genomes["SeqIDLists"],
                                              FamilyList=self.genomes["FamilyList"],
                                              GenusList=self.genomes["GenusList"],
                                              VirusNameList=self.genomes["VirusNameList"]
                                              )
        sigs_data = np.column_stack((TaxoLabelList,
                                    pl1_ref_annotations["PPHMMSignatureTable"]))
        pphmm_names = [f"PPHMM|{ClusterDesc}" for ClusterDesc in pl1_ref_annotations["ClusterDescList"].astype("str")]
        pphmm_names = [pphmm_names[i] if not pphmm_names[i] == "PPHMM|~|" else f"PPHMM|UCF{i}|" for i in range(len(pphmm_names))]
        sigs_df = pd.DataFrame(sigs_data, index=None, columns=
                            ["Virus name"] + \
                            pphmm_names
                )
        sigs_df.to_csv(self.fnames["PphmmAndGomSigs"])
        '''Get PPHMM Locations, save as CSV for location heatmaps'''
        pphmm_names = [i for i in range(0,len(pl1_ref_annotations["NaivePPHMMLocationTable"][0]))] # TODO TEST
        locs_df = pd.DataFrame(np.column_stack((TaxoLabelList, pl1_ref_annotations["NaivePPHMMLocationTable"])), columns = ["Virus name"] + pphmm_names) # TODO FIND OUT WHY DIDNT WORK
        locs_df.to_csv(self.fnames["PphmmLocs"], index=False)
        interim_Rscheme_matrix = shared_norm_pphmm_ratio(TaxoLabelList, self.fnames, labels=[pphmm_names,  [int(i) for i in range(len(TaxoLabelList))]])
        return pphmm_names, interim_Rscheme_matrix

    def call_pphmm_sig_graphs(self, label_order, TaxoLabelList, x_labels, pphmm_names):
        pl1_ref_annotations = {**self.genomes, **self.ref_annotations}
        '''Call all heatmaps for PPHMM location and signatures'''
        pphmm_tables = {}
        y_label_acc_ids = np.array([i.split("_")[0] if not i[0:5] == "Query" else i.split("_")[-1] for i in TaxoLabelList])[label_order]
        pphmm_tables["SharedPphmmRatioMatrix"] = shared_pphmm_ratio(label_order, self.fnames, labels=[x_labels, y_label_acc_ids])
        pphmm_tables["SharedNormPphmmRatioMatrix"] = shared_norm_pphmm_ratio(label_order, self.fnames, labels=[x_labels, y_label_acc_ids])
        pphmm_tables["PphmmLocDistances"], pphmm_tables["PphmLabels"], pphmm_tables["PphmmPos"] = pphmm_loc_distances(self.fnames, pphmm_names, label_order, labels=[x_labels, y_label_acc_ids])
        pphmm_tables["PphmmLocMatrix"] = pphmm_loc_diffs_pairwise(self.fnames, label_order, labels=[x_labels, y_label_acc_ids])
        pphmm_heatmaps = {
            "Shared_norm_pphmm_ratio_matrix": [pphmm_tables["SharedNormPphmmRatioMatrix"], True],
            "Pphmm_locations": [pphmm_tables["PphmmPos"], False],
            "Shared_pphmm_ratio_matrix": [pphmm_tables["SharedPphmmRatioMatrix"], True],
            "Pphmm_loc_distances_pairwise": [pphmm_tables["PphmmLocMatrix"], True],
            "Pphmm_loc_distances": [pphmm_tables["PphmmLocDistances"], False]
        }
        for fname, data in pphmm_heatmaps.items():
            _, _ = supplementary_pphmm_heatmaps(pl1_ref_annotations, TaxoLabelList,
                                                     len(pl1_ref_annotations["SeqIDLists"]), TaxoLabelList, [],
                                                     data[0], f"{self.fnames['OutputDir']}/{fname}.pdf",
                                                     self.fnames, self.payload, pphmm_tables, is_square=data[1],
                                                     )

    def main(self):
        '''Generate GRAViTy dendrogram and heat map	'''
        section_header("Generate GRAViTy dendrogram and heat map")

        '''Estimate virus pairwise distances'''
        pphmm_names, interim_Rscheme_matrix = self.create_supplementary_pphmm_tables()
        DistMat = self.dist_mat_constructor(self.ref_annotations["PPHMMSignatureTable"],
                                           self.ref_annotations["GOMSignatureTable"],
                                           self.ref_annotations["PPHMMLocationTable"],
                                           interim_Rscheme_matrix)

        '''Do plotting'''
        self.dendrogram(DistMat, interim_Rscheme_matrix)
        label_order, taxo_labels, x_labels = self.heatmap_with_dendro(DistMat)

        if self.payload['VirusGrouping'] == True:
            '''(OPT) Group viruses'''
            self.virus_grouping(DistMat)

        '''Call PPHMM Heatmaps'''
        self.call_pphmm_sig_graphs(label_order, taxo_labels, x_labels, pphmm_names)
