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

import os
import numpy as np
import plotly.express as px
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

    def est_pairwise_dists(self):
        '''3/7: Estimate pairwise distances between viruses, return distance matrix (recip. sim matrix)'''
        SimMat = SimilarityMat_Constructor(PPHMMSignatureTable=self.ref_annotations["PPHMMSignatureTable"],
                                           GOMSignatureTable=self.ref_annotations["GOMSignatureTable"],
                                           PPHMMLocationTable=self.ref_annotations["PPHMMLocationTable"],
                                           SimilarityMeasurementScheme=self.payload['SimilarityMeasurementScheme'],
                                           p=self.payload['p']
                                           )
        DistMat = 1 - SimMat
        DistMat[DistMat < 0] = 0
        return DistMat

    def dendrogram(self, DistMat):
        '''4/7: Build dendrogram'''
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
                                       Dendrogram_LinkageMethod=self.payload['Dendrogram_LinkageMethod']
                                       )
        with open(self.fnames['Heatmap_DendrogramFile'], "w") as VirusDendrogram_txt:
            VirusDendrogram_txt.write(VirusDendrogram)

        if self.payload['Bootstrap'] == True:
            '''Compute bootstrap support'''
            progress_msg("Computing Dendrogram Bootstrap Support")
            N_PPHMMs = self.ref_annotations["PPHMMSignatureTable"].shape[1]
            for Bootstrap_i in range(0, self.payload['N_Bootstrap']):
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
                BootstrappedSimMat = SimilarityMat_Constructor(PPHMMSignatureTable=BootstrappedPPHMMSignatureTable,
                                                               GOMSignatureTable=BootstrappedGOMSignatureTable,
                                                               PPHMMLocationTable=BootstrappedPPHMMLocationTable,
                                                               SimilarityMeasurementScheme=self.payload['SimilarityMeasurementScheme'],
                                                               p=self.payload['p']
                                                               )
                BootstrappedDistMat = 1 - BootstrappedSimMat
                BootstrappedDistMat[BootstrappedDistMat < 0] = 0

                '''Generate bs dist matrix tree'''
                BootstrappedVirusDendrogram = DistMat2Tree(DistMat=BootstrappedDistMat,
                                                           LeafList=TaxoLabelList,
                                                           Dendrogram_LinkageMethod=self.payload['Dendrogram_LinkageMethod']
                                                           )
                with open(self.fnames['VirusDendrogramDistFile'], "a") as VirusDendrogramDist_txt:
                    VirusDendrogramDist_txt.write(
                        BootstrappedVirusDendrogram+"\n")

            if os.path.isfile(self.fnames['BootstrappedDendrogramFile']):
                '''Bootstrap gets really upset if trying to overwrite'''
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
        VirusDendrogram = Phylo.read(self.fnames['BootstrappedDendrogramFile'], "newick")

        '''Determine virus order'''
        TaxoLabelList = TaxoLabel_Constructor(SeqIDLists=self.genomes["SeqIDLists"],
                                              FamilyList=self.genomes["FamilyList"],
                                              GenusList=self.genomes["GenusList"],
                                              VirusNameList=self.genomes["VirusNameList"]
                                              )

        _ = VirusDendrogram.ladderize(reverse=True)
        OrderedTaxoLabelList = [
            Clade.name for Clade in VirusDendrogram.get_terminals()]
        VirusOrder = [TaxoLabelList.index(TaxoLabel)
                      for TaxoLabel in OrderedTaxoLabelList]

        '''Re-order the distance matrix'''
        OrderedDistMat = DistMat[VirusOrder][:, VirusOrder]

        '''Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff'''
        N_InternalNodes = len(VirusDendrogram.get_nonterminals())
        for InternalNode_i in range(N_InternalNodes):
            try:
                if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence < self.payload['Heatmap_DendrogramSupport_Cutoff'] or np.isnan(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence):
                    continue
                else:
                    VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(
                        VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
            except:
                continue

        '''Labels, label positions, and ticks'''
        ClassDendrogram_grp = VirusDendrogram
        ClassDendrogram_label = Phylo.read(  # Needs to be separately read in as Bio entities are linked
            self.fnames['BootstrappedDendrogramFile'], "newick")

        '''Construct taxo labels in ordered list for heatmap axis labels'''
        _, LineList_major = make_labels(ClassDendrogram_grp, zip(
            TaxoLabelList, self.genomes["TaxoGroupingList"]))
        _, LineList_minor = make_labels(ClassDendrogram_label, zip(
            TaxoLabelList, self.genomes["VirusNameList"]))
        ClassLabelList_x, ClassLabelList_y = split_labels(OrderedTaxoLabelList)

        '''Plot configuration'''
        Heatmap_width = float(12)
        Heatmap_height = Heatmap_width
        TaxoLable_space = 1.00

        CBar_Heatmap_gap = 0.05
        CBar_width = Heatmap_width
        CBar_height = 0.25
        CBarLable_space = 0.25

        Dendrogram_width = Heatmap_width/3
        Dendrogram_height = Heatmap_height
        Dendrogram_Heatmap_gap = 0.1

        ScaleBar_Dendrogram_gap = CBar_Heatmap_gap
        ScaleBar_width = Dendrogram_width
        ScaleBar_height = CBar_height
        ScaleBarLable_space = CBarLable_space

        Outer_margin = 0.5
        FontSize = 6
        dpi=600

        if len(ClassLabelList_x) >= 200:
            '''Reduce draw parameter size if large n'''
            FontSize = FontSize/2
            dpi=dpi*1.5

        elif len(ClassLabelList_x) >= 400:
            FontSize = FontSize/4
            dpi=dpi*3


        Fig_width = Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap + \
            Heatmap_width + TaxoLable_space + Outer_margin
        Fig_height = Outer_margin + CBarLable_space + CBar_height + \
            CBar_Heatmap_gap + Heatmap_height + TaxoLable_space + Outer_margin

        ax_Dendrogram_L = Outer_margin/Fig_width
        ax_Dendrogram_B = (Outer_margin + ScaleBarLable_space +
                           ScaleBar_height + ScaleBar_Dendrogram_gap)/Fig_height
        ax_Dendrogram_W = Dendrogram_width/Fig_width
        ax_Dendrogram_H = Dendrogram_height/Fig_height

        ax_ScaleBar_L = Outer_margin/Fig_width
        ax_ScaleBar_B = (Outer_margin + ScaleBarLable_space)/Fig_height
        ax_ScaleBar_W = ScaleBar_width/Fig_width
        ax_ScaleBar_H = ScaleBar_height/Fig_height

        ax_Heatmap_L = (Outer_margin + Dendrogram_width +
                        Dendrogram_Heatmap_gap)/Fig_width
        ax_Heatmap_B = (Outer_margin + CBarLable_space +
                        CBar_height + CBar_Heatmap_gap)/Fig_height
        ax_Heatmap_W = Heatmap_width/Fig_width
        ax_Heatmap_H = Heatmap_height/Fig_height

        ax_CBar_L = (Outer_margin + Dendrogram_width +
                     Dendrogram_Heatmap_gap)/Fig_width
        ax_CBar_B = (Outer_margin + CBarLable_space)/Fig_height
        ax_CBar_W = CBar_width/Fig_width
        ax_CBar_H = CBar_height/Fig_height

        '''Plot the heat map'''
        fig = plt.figure(figsize=(Fig_width, Fig_height), dpi=dpi)

        '''Dendrogram'''
        ax_Dendrogram = fig.add_axes(
            [ax_Dendrogram_L, ax_Dendrogram_B, ax_Dendrogram_W, ax_Dendrogram_H], frame_on=False, facecolor="white")

        try:
            Phylo.draw(VirusDendrogram, label_func=lambda x: "",
                    do_show=False,  axes=ax_Dendrogram)
        except Exception as ex:
            raise SystemExit(f"ERROR: {ex}\nThis will usually occur when there's been an error with bootstrapping. Try disabling this feature and trying again.")

        VirusDendrogramDepth = max(
            [v for k, v in VirusDendrogram.depths().items()])
        ax_Dendrogram		.set_xlim(
            [(VirusDendrogramDepth - 1), VirusDendrogramDepth])
        ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
        ax_Dendrogram		.set_axis_off()

        '''Dendrogram scale bar'''
        ax_ScaleBar = fig.add_axes(
            [ax_ScaleBar_L, ax_ScaleBar_B, ax_ScaleBar_W, ax_ScaleBar_H], frame_on=False, facecolor="white")
        ax_ScaleBar.plot([0, 1], [0, 0], 'k-')
        ScaleBarTicks = [0, 0.25, 0.5, 0.75, 1]
        for Tick in ScaleBarTicks:
            ax_ScaleBar.plot([Tick, Tick], [-0.05, 0.05], 'k-')

        ax_ScaleBar			.set_xlim([1, 0])
        ax_ScaleBar			.set_xticks(ScaleBarTicks)
        ax_ScaleBar			.set_xticklabels(
            list(map(str, ScaleBarTicks)), rotation=0, size=FontSize)
        ax_ScaleBar			.set_xlabel('Distance', rotation=0, size=FontSize+2)
        ax_ScaleBar			.xaxis.set_label_position('bottom')
        ax_ScaleBar			.tick_params(top=False,
                                   bottom=False,
                                   left=False,
                                   right=False,
                                   labeltop=False,
                                   labelbottom=True,
                                   labelleft=False,
                                   labelright=False,
                                   direction='out')

        '''Heatmap'''
        ax_Heatmap = fig.add_axes(
            [ax_Heatmap_L, ax_Heatmap_B, ax_Heatmap_W, ax_Heatmap_H], frame_on=True, facecolor="white")
        Heatmap_Graphic = ax_Heatmap.imshow(
            OrderedDistMat, cmap='magma', aspect='auto', vmin=0, vmax=1, interpolation='none')

        '''Draw grouping major & minor lines'''
        for l in LineList_major:
            ax_Heatmap.axvline(l, color='k', lw=0.4)
            ax_Heatmap.axhline(l, color='k', lw=0.4)

        for l in LineList_minor:
            ax_Heatmap.axvline(l, color='gray', lw=0.2)
            ax_Heatmap.axhline(l, color='gray', lw=0.2)

        '''Draw gridlines for individual samples'''
        TickLocList = np.arange(0, len(ClassLabelList_y))

        ax_Heatmap			.set_xticks(TickLocList)
        ax_Heatmap			.set_xticklabels(
            ClassLabelList_x, rotation=90, size=FontSize)

        ax_Heatmap			.set_yticks(TickLocList)
        ax_Heatmap			.set_yticklabels(
            ClassLabelList_y, rotation=0, size=FontSize)
        ax_Heatmap			.tick_params(top=True,
                                  bottom=False,
                                  left=False,
                                  right=True,
                                  labeltop=True,
                                  labelbottom=False,
                                  labelleft=False,
                                  labelright=True,
                                  direction='out')

        '''Heatmap colourbars'''
        ax_CBar = fig.add_axes(
            [ax_CBar_L, ax_CBar_B, ax_CBar_W, ax_CBar_H], frame_on=True, facecolor="white")
        CBar_Graphic = fig.colorbar(
            Heatmap_Graphic, cax=ax_CBar, orientation="horizontal", ticks=[0, 0.25, 0.50, 0.75, 1])
        CBar_Graphic		.ax.set_xticklabels(
            ['0', '0.25', '0.50', '0.75', '1'], rotation=0, size=FontSize)
        CBar_Graphic		.ax.set_xlabel('Distance', rotation=0, size=FontSize+2)
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
        plt.savefig(self.fnames['HeatmapWithDendrogramFile'], format="pdf", bbox_inches = "tight", dpi=dpi)
        return VirusOrder

    def virus_grouping(self, DistMat):
        '''7/7: (OPT) Group viruses via Thiels-U and other metrics; save as txt.'''
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

    def main(self):
        '''Generate GRAViTy dendrogram and heat map	'''
        section_header("Generate GRAViTy dendrogram and heat map")
        '''1/7: Make fpaths'''
        mkdir_pl1_graphs(self.fnames, self.payload)

        '''3/7: Estimate virus pairwise distances'''
        DistMat = self.est_pairwise_dists()

        '''X/7: Do plotting'''
        self.dendrogram(DistMat)
        label_order = self.heatmap_with_dendro(DistMat)

        if self.payload['VirusGrouping'] == True:
            '''7/7: (OPT) Group viruses'''
            self.virus_grouping(DistMat)

        ########################################################################################
        ################## New feats in development

        self.shared_pphmm_ratio(label_order)
        self.shared_norm_pphmm_ratio(label_order)


    def shared_pphmm_ratio(self, label_order):
        out_fname = f"{self.fnames['OutputDir']}/shared_pphmms.png" # TODO Parameterise
        df = pd.read_csv(f"{self.fnames['OutputDir']}/PPHMMandGOMsignatures.csv", index_col=0) # TODO Parameterise
        df = df.drop(columns=["Virus name", "Accession", "Order", "Family", "Subfamily", "Genus", "Class", df.columns[-1]])
        df = df.apply(pd.to_numeric)
        arr = np.asarray(df.reindex(label_order))

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

        fig = px.imshow(shared_pphmms, x=self.genomes["VirusNameList"][label_order],y=list([i[0] for i in self.genomes["SeqIDLists"][label_order]]),title="Shared N PPHMMs")
        fig.layout.height=1000; fig.layout.width=1000
        fig.write_image(out_fname)

    def shared_norm_pphmm_ratio(self, label_order):
        out_fname = f"{self.fnames['OutputDir']}/shared_norm_pphmm_ratio.png" # TODO Parameterise
        df = pd.read_csv(f"{self.fnames['OutputDir']}/PPHMMandGOMsignatures.csv", index_col=0) # TODO Parameterise
        df = df.drop(columns=["Virus name", "Accession", "Order", "Family", "Subfamily", "Genus", "Class", df.columns[-1]])
        df = df.apply(pd.to_numeric)
        arr = np.asarray(df.reindex(label_order))
        def func(a,b, max_possible):
            hits = len([a[idx] for idx in range(a.shape[0]) if not a[idx] == 0 and not b[idx] == 0])
            return min([hits, max_possible]) / max([hits,max_possible])

        l = arr.shape[0]
        shared_pphmm_ratio = np.zeros((l,l))
        for idx_a, a in enumerate(arr):
            for idx_b, b in enumerate(arr):
                max_possible = min([np.count_nonzero(a), np.count_nonzero(b)])
                shared_pphmm_ratio[idx_b,idx_a] = func(a,b, max_possible)

        fig = px.imshow(shared_pphmm_ratio, x=self.genomes["VirusNameList"][label_order],y=list([i[0] for i in self.genomes["SeqIDLists"][label_order]]),title="Shared Normalised PPHMM Ratio (pairwise: n common PPHMMs / n PPHMMs)")
        fig.layout.height=1000; fig.layout.width=1000
        fig.write_image(out_fname)

        # RM < TODO probably don't need these now it's all in ref annots? FROM MAIN LOOP
        # '''2/7: Retrieve additional sparse coordinate matrices, if present'''
        # if "PPHMMSignatureTable_coo" in self.ref_annotations.keys():
        #     self.ref_annotations["PPHMMSignatureTable_coo"] = self.ref_annotations["PPHMMSignatureTable_coo"].toarray(
        #     )
        # if "PPHMMLocationTable_coo" in self.ref_annotations.keys():
        #     self.ref_annotations["PPHMMLocationTable_coo"] = self.ref_annotations["PPHMMLocationTable_coo"].toarray(
        #     )
