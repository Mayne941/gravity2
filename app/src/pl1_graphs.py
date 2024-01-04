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
from app.utils.shared_pphmm_graphs import shared_norm_pphmm_ratio, shared_pphmm_ratio

import os
import numpy as np
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

    def est_pairwise_dists(self):
        '''Estimate pairwise distances between viruses, return distance matrix (recip. sim matrix)'''
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
                                       Dendrogram_LinkageMethod=self.payload['Dendrogram_LinkageMethod']
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
        _ = VirusDendrogram.ladderize(reverse=True)
        TaxoLabelList = TaxoLabel_Constructor(SeqIDLists=self.genomes["SeqIDLists"],
                                              FamilyList=self.genomes["FamilyList"],
                                              GenusList=self.genomes["GenusList"],
                                              VirusNameList=self.genomes["VirusNameList"]
                                              )

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
        hmap_params, fig, ax_Dendrogram, ax_Heatmap = get_hmap_params(len(ClassLabelList_x))

        try:
            Phylo.draw(VirusDendrogram, label_func=lambda x: "",
                    do_show=False,  axes=ax_Dendrogram)
        except Exception as ex:
            raise raise_gravity_error(f"ERROR: {ex}\nThis will usually occur when there's been an error with bootstrapping. Try disabling this feature and trying again.")

        VirusDendrogramDepth = max(
            [v for k, v in VirusDendrogram.depths().items()])
        ax_Dendrogram		.set_xlim(
            [(VirusDendrogramDepth - 1), VirusDendrogramDepth])
        ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
        ax_Dendrogram		.set_axis_off()

        '''Draw heatmap elements on axes'''
        Heatmap_Graphic = ax_Heatmap.imshow(
            OrderedDistMat, cmap='magma', aspect='auto', vmin=0, vmax=1, interpolation='none')

        '''Draw grouping major & minor lines'''
        ax_Heatmap = construct_hmap_lines(ax_Heatmap, LineList_major, LineList_minor,
                                          hmap_params, ClassLabelList_x, ClassLabelList_y,
                                          TickLocList = np.arange(0, len(ClassLabelList_y)))

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
        return VirusOrder

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

    def call_pphmm_sig_graphs(self, label_order):
        '''Call PPHMM Heatmaps'''
        labels = [self.genomes["VirusNameList"][label_order], list([i[0] for i in self.genomes["SeqIDLists"][label_order]])]
        shared_pphmm_ratio(label_order, self.fnames, labels)
        shared_norm_pphmm_ratio(label_order, self.fnames, labels)

    def main(self):
        '''Generate GRAViTy dendrogram and heat map	'''
        section_header("Generate GRAViTy dendrogram and heat map")

        '''Estimate virus pairwise distances'''
        DistMat = self.est_pairwise_dists()

        '''Do plotting'''
        self.dendrogram(DistMat)
        label_order = self.heatmap_with_dendro(DistMat)

        if self.payload['VirusGrouping'] == True:
            '''(OPT) Group viruses'''
            self.virus_grouping(DistMat)

        '''Call PPHMM Heatmaps'''
        self.call_pphmm_sig_graphs(label_order)
