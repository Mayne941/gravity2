from Bio import Phylo
from copy import copy
from scipy.cluster.hierarchy import linkage, fcluster
from collections import Counter
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
import sys
import random
import string
import operator
import pickle
import re

from app.utils.dist_mat_to_tree import DistMat2Tree
from app.utils.similarity_matrix_constructor import SimilarityMat_Constructor
from app.utils.taxo_label_constructor import TaxoLabel_Constructor
from app.utils.pphmm_signature_table_constructor import PPHMMSignatureTable_Constructor
from app.utils.gomdb_constructor import GOMDB_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.virus_grouping_estimator import VirusGrouping_Estimator
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_pickle
from app.utils.highest_posterior_density import hpd
from app.utils.classification_utils import PairwiseSimilarityScore_Cutoff_Dict_Constructor, TaxonomicAssignmentProposerAndEvaluator
from app.utils.make_heatmap_labels import make_labels, split_labels
from app.utils.generate_fnames import generate_file_names
from app.utils.error_handlers import error_handler_virus_classifier
from app.utils.stdout_utils import progress_msg
from app.utils.shell_cmds import shell
from app.utils.mkdirs import mkdir_virus_classifier
from app.utils.heatmap_params import get_blue_cmap, get_red_cmap, get_purple_cmap, get_hmap_params, construct_hmap_lines
from app.utils.error_handlers import raise_gravity_error

matplotlib.use('agg')
sys.setrecursionlimit(10000)


class VirusClassificationAndEvaluation:
    def __init__(self,
                 payload,
                 ExpDir,
                 ):
        self.payload = payload
        self.fnames = generate_file_names(payload, ExpDir, Pl2=True)
        self.genomes = retrieve_genome_vars(self.fnames['ReadGenomeDescTablePickle'])
        self.ucf_annots = retrieve_pickle(self.fnames['UcfAnnotatorPickle'])
        if "PPHMMSignatureTable_Dict_coo" in self.ucf_annots.keys(): # RM < TODO Optimal if could remove RefVirusGroup (_) from UCF Annotator to avoid this
            self.ucf_annots["PPHMMSignatureTable_Dict"] = [PPHMMSignatureTable_coo.toarray() for _, PPHMMSignatureTable_coo in self.ucf_annots["PPHMMSignatureTable_Dict_coo"].items()][0]
        if "PPHMMLocationTable_Dict_coo" in self.ucf_annots.keys():
            self.ucf_annots["PPHMMLocationTable_Dict"] = [PPHMMLocationTable_coo.toarray() for _, PPHMMLocationTable_coo in self.ucf_annots["PPHMMLocationTable_Dict_coo"].items()][0]
        self.N_UcfViruses = len(self.genomes["SeqIDLists"])
        self.TaxoLabelList_UcfVirus = [f"Query_{i+1}_{self.genomes['SeqIDLists'][i][0]}" for i in range(self.N_UcfViruses)]
        self.final_results = {}

    def use_ucf_virus_pphmms(self):
        '''4/8: Scan unclassified viruses against their PPHMM DB to create additional PPHMMSignatureTable, and PPHMMLocationTable'''
        progress_msg('Generating PPHMMSignatureTable and PPHMMLocationTable tables for Unclassified Viruses')
        PPHMMSignatureTable_UcfVirusVSUcfDB, PPHMMLocationTable_UcfVirusVSUcfDB = PPHMMSignatureTable_Constructor( # TODO Surplus now
                self.genomes,
                self.payload,
                self.fnames,
                GenomeSeqFile=self.payload['GenomeSeqFile_UcfVirus'],
                HMMER_PPHMMDB=self.fnames['HMMER_PPHMMDB_UcfVirus'],
                Pl2=True
            )


        '''Update unclassified viruses' PPHMMSignatureTables, and PPHMMLocationTables'''
        self.ucf_annots["PPHMMSignatureTable_Dict"] = np.hstack((self.ucf_annots["PPHMMSignatureTable_Dict"], PPHMMSignatureTable_UcfVirusVSUcfDB))
        self.ucf_annots["PPHMMLocationTable_Dict"] = np.hstack((self.ucf_annots["PPHMMLocationTable_Dict"], PPHMMLocationTable_UcfVirusVSUcfDB))

    def classify(self):
        '''6/8: Classify viruses'''
        progress_msg(f"Classifying unknown viruses")
        MaxSimScoreTable, TaxoOfMaxSimScoreTable, TaxoAssignmentTable, \
            PhyloStatTable = np.zeros((self.N_UcfViruses, 0)), np.zeros(
                (self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0))

        '''Load in genomes variables'''
        pl1_ref_annotations = {**retrieve_genome_vars(self.fnames['Pl1ReadDescTablePickle']),
                                **retrieve_pickle(self.fnames['Pl1RefAnnotatorPickle'])}
        N_RefViruses = len(pl1_ref_annotations["SeqIDLists"])
        if "PPHMMSignatureTable_coo" in pl1_ref_annotations.keys():
            pl1_ref_annotations["PPHMMSignatureTable"] = pl1_ref_annotations["PPHMMSignatureTable_coo"].toarray(
            )
        if "PPHMMLocationTable_coo" in pl1_ref_annotations.keys():
            pl1_ref_annotations["PPHMMLocationTable"] = pl1_ref_annotations["PPHMMLocationTable_coo"].toarray(
            )

        if self.payload['UseUcfVirusPPHMMs']:
            '''Scan reference viruses against the PPHMM database of unclassified viruses to generate additional PPHMMSignatureTable, and PPHMMLocationTable'''
            '''Make OR update HMMER_hmmscanDir'''
            PPHMMSignatureTable_UcfVirusVSUcfDB, PPHMMLocationTable_UcfVirusVSUcfDB = PPHMMSignatureTable_Constructor( # TODO Surplus now
                    retrieve_genome_vars(self.fnames['Pl1ReadDescTablePickle']),
                    self.payload,
                    self.fnames,
                    GenomeSeqFile=self.payload['GenomeSeqFiles_RefVirus'],
                    HMMER_PPHMMDB=self.fnames['HMMER_PPHMMDB_UcfVirus'],
                    Pl2=True
                )

            pl1_ref_annotations["PPHMMSignatureTable"] = np.hstack(
                (pl1_ref_annotations["PPHMMSignatureTable"], PPHMMSignatureTable_UcfVirusVSUcfDB))
            pl1_ref_annotations["PPHMMLocationTable"] = np.hstack(
                (pl1_ref_annotations["PPHMMLocationTable"],  PPHMMLocationTable_UcfVirusVSUcfDB))
            UpdatedGOMDB_RefVirus = GOMDB_Constructor(
                pl1_ref_annotations["TaxoGroupingList"], pl1_ref_annotations["PPHMMLocationTable"], pl1_ref_annotations["GOMIDList"])

            pl1_ref_annotations["GOMSignatureTable"] = GOMSignatureTable_Constructor(
                pl1_ref_annotations["PPHMMLocationTable"], UpdatedGOMDB_RefVirus, pl1_ref_annotations["GOMIDList"])

            '''Update unclassified viruses' GOMSignatureTable'''
            self.ucf_annots["GOMSignatureTable_Dict"] = GOMSignatureTable_Constructor(
                self.ucf_annots["PPHMMLocationTable_Dict"], UpdatedGOMDB_RefVirus, pl1_ref_annotations["GOMIDList"])

        '''Build the dendrogram, including all sequences'''
        '''Generate TaxoLabelList of reference viruses'''
        TaxoLabelList_RefVirus = TaxoLabel_Constructor(
            pl1_ref_annotations["SeqIDLists"], pl1_ref_annotations["FamilyList"], pl1_ref_annotations["GenusList"], pl1_ref_annotations["VirusNameList"])

        '''Compute pairwise distances'''
        PPHMMSignatureTable_AllVirus = np.vstack(
            (pl1_ref_annotations["PPHMMSignatureTable"], self.ucf_annots["PPHMMSignatureTable_Dict"]))
        GOMSignatureTable_AllVirus = np.vstack(
            (pl1_ref_annotations["GOMSignatureTable"],   self.ucf_annots["GOMSignatureTable_Dict"]))
        PPHMMLocationTable_AllVirus = np.vstack(
            (pl1_ref_annotations["PPHMMLocationTable"],  self.ucf_annots["PPHMMLocationTable_Dict"]))
        TaxoLabelList_AllVirus = TaxoLabelList_RefVirus + self.TaxoLabelList_UcfVirus
        SimMat = SimilarityMat_Constructor(PPHMMSignatureTable_AllVirus, GOMSignatureTable_AllVirus,
                                            PPHMMLocationTable_AllVirus, self.payload['SimilarityMeasurementScheme'], self.payload['p'])
        DistMat = 1 - SimMat
        DistMat[DistMat < 0] = 0

        '''Generate dendrogram'''
        VirusDendrogram = DistMat2Tree(
            DistMat, TaxoLabelList_AllVirus, self.payload['Dendrogram_LinkageMethod'])

        with open(self.fnames['VirusDendrogramFile'], "w") as VirusDendrogram_txt:
            VirusDendrogram_txt.write(VirusDendrogram)

        '''Compute similarity cut off for each taxonomic class'''
        self.final_results["PairwiseSimilarityScore_Cutoff_Dict"] = PairwiseSimilarityScore_Cutoff_Dict_Constructor(
            SimMat[:N_RefViruses][:, :N_RefViruses], pl1_ref_annotations["TaxoGroupingList"], self.payload['N_PairwiseSimilarityScores'])

        '''Propose a taxonomic class to each unclassified virus and evaluate the taxonomic assignments'''
        MaxSimScoreList, TaxoOfMaxSimScoreList, \
            TaxoAssignmentList, PhyloStatList = TaxonomicAssignmentProposerAndEvaluator(SimMat[-self.N_UcfViruses:][:, :N_RefViruses],
                                                                                        pl1_ref_annotations["TaxoGroupingList"],
                                                                                        VirusDendrogram,
                                                                                        TaxoLabelList_RefVirus,
                                                                                        self.TaxoLabelList_UcfVirus,
                                                                                        self.final_results["PairwiseSimilarityScore_Cutoff_Dict"])

        try:
            self.final_results["MaxSimScoreTable"] = np.column_stack(
                (MaxSimScoreTable, MaxSimScoreList))
            self.final_results["TaxoOfMaxSimScoreTable"] = np.column_stack(
                (TaxoOfMaxSimScoreTable, TaxoOfMaxSimScoreList))
            self.final_results["TaxoAssignmentTable"] = np.column_stack(
                (TaxoAssignmentTable, TaxoAssignmentList))
            self.final_results["PhyloStatTable"] = np.column_stack(
                (PhyloStatTable, PhyloStatList))

        except ValueError:
            '''If rows omitted then change shape of array to fit'''
            TaxoAssignmentTable = np.zeros((len(TaxoAssignmentList), 0))
            PhyloStatTable = np.zeros((len(TaxoAssignmentList), 0))
            MaxSimScoreTable = np.zeros((len(MaxSimScoreList), 0))
            TaxoOfMaxSimScoreTable = np.zeros(
                (len(TaxoOfMaxSimScoreList), 0))

            self.final_results["MaxSimScoreTable"] = np.column_stack(
                (MaxSimScoreTable, MaxSimScoreList))
            self.final_results["TaxoOfMaxSimScoreTable"] = np.column_stack(
                (TaxoOfMaxSimScoreTable, TaxoOfMaxSimScoreList))
            self.final_results["TaxoAssignmentTable"] = np.column_stack(
                (TaxoAssignmentTable, TaxoAssignmentList))
            self.final_results["PhyloStatTable"] = np.column_stack(
                (PhyloStatTable, PhyloStatList))

        if self.payload['VirusGrouping']:
            '''(OPT) Virus grouping'''
            self.do_virus_groupings(DistMat, N_RefViruses, pl1_ref_annotations)

        if self.payload['Bootstrap']:
            '''(OPT) Bootstrap stats'''
            self.construct_bootstrap_distributions(pl1_ref_annotations, PPHMMSignatureTable_AllVirus,
                                                    PPHMMLocationTable_AllVirus, N_RefViruses, TaxoLabelList_AllVirus, TaxoLabelList_RefVirus)

        '''Draw hm/dendro'''
        self.heatmap_with_dendrogram(pl1_ref_annotations, TaxoLabelList_AllVirus,
                                        DistMat, N_RefViruses, TaxoLabelList_RefVirus, TaxoAssignmentList)


    def do_virus_groupings(self, DistMat, N_RefViruses, pl1_ref_annotations) -> None:
        '''Accessory fn to CLASSIFY (6):  Group viruses with scipy fcluster/linkage, save results to txt'''
        self.final_results["VirusGrouping"] = {}

        '''Compute distance cutoff'''
        _, OptDistance_Cutoff, CorrelationScore, Theils_u_TaxoGroupingListGivenPred, \
            Theils_u_PredGivenTaxoGroupingList = VirusGrouping_Estimator(
                DistMat[:N_RefViruses][:, :N_RefViruses], self.payload['Dendrogram_LinkageMethod'], pl1_ref_annotations["TaxoGroupingList"])

        '''Virus grouping'''
        VirusGroupingList = fcluster(linkage(DistMat[np.triu_indices_from(
            DistMat, k=1)], self.payload['Dendrogram_LinkageMethod']), OptDistance_Cutoff, 'distance')
        self.final_results["VirusGrouping"]["VirusGroupingList"] = VirusGroupingList
        self.final_results["VirusGrouping"]["OptDistance_Cutoff"] = OptDistance_Cutoff
        self.final_results["VirusGrouping"]["CorrelationScore"] = CorrelationScore
        self.final_results["VirusGrouping"]["Theils_u_TaxoGroupingListGivenPred"] = Theils_u_TaxoGroupingListGivenPred
        self.final_results["VirusGrouping"]["Theils_u_PredGivenTaxoGroupingList"] = Theils_u_PredGivenTaxoGroupingList


        '''Save result to txt and CSV files'''
        output_data = np.column_stack((list(map(", ".join, pl1_ref_annotations["SeqIDLists"].tolist()+self.genomes["SeqIDLists"].tolist())), pl1_ref_annotations["FamilyList"].tolist( )+[""]*self.N_UcfViruses, pl1_ref_annotations["GenusList"].tolist( )+[""]*self.N_UcfViruses, pl1_ref_annotations["VirusNameList"].tolist( )+self.genomes["VirusNameList"].tolist(), pl1_ref_annotations["TaxoGroupingList"].tolist( )+[""]*self.N_UcfViruses, VirusGroupingList, ))
        output_header ="Sequence identifier\tFamily\tGenus\tVirus name\tClass\tGrouping"
        np.savetxt(fname=self.fnames['VirusGroupingFile'],
                   X=output_data,
                   fmt='%s',
                   delimiter="\t",
                   header=output_header)

        with open(self.fnames['VirusGroupingFile'], "a") as VirusGrouping_txt:
            VirusGrouping_txt.write("\n" +
                                    f"Distance cut off: {OptDistance_Cutoff}\n"
                                    f"Theil's uncertainty correlation for the reference assignments given the predicted grouping U(Ref|Pred): {Theils_u_TaxoGroupingListGivenPred}\n"
                                    f"Theil's uncertainty correlation for the predicted grouping given the reference assignments U(Pred|Ref): {Theils_u_PredGivenTaxoGroupingList}\n"
                                    f"Symmetrical Theil's uncertainty correlation between the reference assignments and the predicted grouping U(Ref, Pred): {CorrelationScore}\n"
                                    f"U(X|Y) == 1 means that knowing Y implies a perfect knowledge of X, but not vice-versa\n"
                                    f"U(X,Y) == 1 means that knowing Y implies a perfect knowledge of X and vice-versa\n"
                                    )
        df = pd.DataFrame(output_data)
        df.columns = output_header.split("\t")
        df.to_csv(self.fnames['VirusGroupingFile'].replace(".txt", ".csv"), index=False)

    def construct_bootstrap_distributions(self, pl1_ref_annotations, PPHMMSignatureTable_AllVirus, PPHMMLocationTable_AllVirus, N_RefViruses, TaxoLabelList_AllVirus, TaxoLabelList_RefVirus):
        '''Accessory fn to CLASSIFY (6). RM << DOCSTRING'''
        '''Construct result distributions by using the bootstrapping technique'''
        '''Define path to dendrogram distribution file'''
        progress_msg(f"Constructing Bootstrap Distributions")
        BootsrappedTaxoAssignmentTable, BootsrappedMaxSimScoreTable, BootsrappedTaxoOfMaxSimScoreTable, \
            BootsrappedPhyloStatTable = np.zeros((self.N_UcfViruses, 0)), np.zeros(
                (self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0))
        N_PPHMMs = PPHMMSignatureTable_AllVirus.shape[1]
        self.final_results["PairwiseSimilarityScore_CutoffDist_Dict"] = {Bootstrap_i: {} for Bootstrap_i in range(self.payload['N_Bootstrap'])}

        for Bootstrap_i in range(0, self.payload['N_Bootstrap']):
            '''Bootstrap the data'''
            print(f"Round {Bootstrap_i+1}")

            '''Construct bootstrapped PPHMMSignatureTable and PPHMMLocationTable'''
            PPHMM_IndexList = sorted(np.random.choice(
                list(range(N_PPHMMs)), N_PPHMMs, replace=True))
            BootstrappedPPHMMSignatureTable = PPHMMSignatureTable_AllVirus[:,
                                                                           PPHMM_IndexList]
            BootstrappedPPHMMLocationTable = PPHMMLocationTable_AllVirus[:,
                                                                         PPHMM_IndexList]
            BootstrappedGOMSignatureTable = None

            if "G" in self.payload['SimilarityMeasurementScheme']:
                '''Construct bootstrapped GOMSignatureTable'''
                BootstrappedGOMDB = GOMDB_Constructor(
                    pl1_ref_annotations["TaxoGroupingList"], BootstrappedPPHMMLocationTable[:N_RefViruses], pl1_ref_annotations["GOMIDList"])
                BootstrappedGOMSignatureTable = GOMSignatureTable_Constructor(
                    BootstrappedPPHMMLocationTable, BootstrappedGOMDB, pl1_ref_annotations["GOMIDList"])

            '''Construct a dendrogram from the bootstrapped data'''
            BootstrappedSimMat = SimilarityMat_Constructor(
                BootstrappedPPHMMSignatureTable, BootstrappedGOMSignatureTable, BootstrappedPPHMMLocationTable, self.payload['SimilarityMeasurementScheme'], self.payload['p'])
            BootstrappedDistMat = 1 - BootstrappedSimMat
            BootstrappedDistMat[BootstrappedDistMat < 0] = 0
            BootstrappedVirusDendrogram = DistMat2Tree(
                BootstrappedDistMat, TaxoLabelList_AllVirus, self.payload['Dendrogram_LinkageMethod'])

            with open(self.fnames['VirusDendrogramDistFile'], "a") as VirusDendrogramDist_txt:
                VirusDendrogramDist_txt.write(BootstrappedVirusDendrogram+"\n")

            '''Compute similarity cut off for each taxonomic class based on the bootstrapped data'''
            self.final_results["PairwiseSimilarityScore_CutoffDist_Dict"][Bootstrap_i] = PairwiseSimilarityScore_Cutoff_Dict_Constructor(BootstrappedSimMat[:N_RefViruses][:, :N_RefViruses],
                                                                                                                                       pl1_ref_annotations[
                                                                                                                                           "TaxoGroupingList"],
                                                                                                                                       self.payload['N_PairwiseSimilarityScores'])

            '''Propose a taxonomic class to each unclassified virus and evaluate the taxonomic assignments based on the bootstrapped data'''
            BootsrappedMaxSimScoreList, BootsrappedTaxoOfMaxSimScoreList, BootsrappedTaxoAssignmentList, \
                BootsrappedPhyloStatList = TaxonomicAssignmentProposerAndEvaluator(BootstrappedSimMat[-self.N_UcfViruses:][:, :N_RefViruses],
                                                                                   pl1_ref_annotations["TaxoGroupingList"],
                                                                                   BootstrappedVirusDendrogram,
                                                                                   TaxoLabelList_RefVirus,
                                                                                   self.TaxoLabelList_UcfVirus,
                                                                                   self.final_results["PairwiseSimilarityScore_CutoffDist_Dict"][Bootstrap_i])

            BootsrappedMaxSimScoreTable = np.column_stack(
                (BootsrappedMaxSimScoreTable, BootsrappedMaxSimScoreList))
            BootsrappedTaxoOfMaxSimScoreTable = np.column_stack(
                (BootsrappedTaxoOfMaxSimScoreTable, BootsrappedTaxoOfMaxSimScoreList))
            BootsrappedTaxoAssignmentTable = np.column_stack(
                (BootsrappedTaxoAssignmentTable, BootsrappedTaxoAssignmentList))
            BootsrappedPhyloStatTable = np.column_stack(
                (BootsrappedPhyloStatTable, BootsrappedPhyloStatList))

        '''Update globals with bootstrap results'''
        self.final_results["MaxSimScoreDistTable"] = np.array([i for i in BootsrappedMaxSimScoreTable.T])
        self.final_results["TaxoOfMaxSimScoreDistTable"] = np.array([i for i in BootsrappedTaxoOfMaxSimScoreTable.T])
        self.final_results["TaxoAssignmentDistTable"] = np.array([i for i in BootsrappedTaxoAssignmentTable.T])
        self.final_results["PhyloStatDistTable"] = np.array([i for i in BootsrappedPhyloStatTable.T])

        '''Create a bootstrapped dendrogram'''
        if self.payload['Bootstrap_method'] == "booster":
            '''Can't call error handler here as tool writes to stderr'''
            shell(f"./booster_linux64 -i {self.fnames['VirusDendrogramFile']} -b {self.fnames['VirusDendrogramDistFile']} -o {self.fnames['BootstrappedVirusDendrogramFile']} -@ {self.payload['N_CPUs']} ")

        elif self.payload['Bootstrap_method'] == "sumtrees":
            shell(f"sumtrees.py --decimals=2 --no-annotations --preserve-underscores --force-rooted --output-tree-format=newick --output-tree-filepath={self.fnames['BootstrappedVirusDendrogramFile']} --target={self.fnames['VirusDendrogramFile']} {self.fnames['VirusDendrogramDistFile']}")

        else:
            print("WARNING ** Bootstrap dendrogram not constructed: 'Bootstrap_method' can either be 'booster' or 'sumtrees'.")

    def heatmap_with_dendrogram(self, pl1_ref_annotations, TaxoLabelList_AllVirus, DistMat, N_RefViruses, TaxoLabelList_RefVirus, TaxoAssignmentList):
        '''Construct GRAViTy heat map with dendrogram'''
        progress_msg("Constructing GRAViTy Heatmap and Dendrogram")
        '''Load the tree'''
        if self.payload['Bootstrap']:
            VirusDendrogramFile = self.fnames['BootstrappedVirusDendrogramFile']
        else:
            VirusDendrogramFile = self.fnames['VirusDendrogramFile']

        VirusDendrogram = Phylo.read(VirusDendrogramFile, "newick")

        '''Determine virus order'''
        _ = VirusDendrogram.ladderize(reverse=True)
        OrderedTaxoLabelList = [
            Clade.name for Clade in VirusDendrogram.get_terminals()]
        VirusOrder = [TaxoLabelList_AllVirus.index(
            TaxoLabel) for TaxoLabel in OrderedTaxoLabelList]

        '''Re-order the distance matrix for heatmap squares'''
        OrderedDistMat = DistMat[VirusOrder][:, VirusOrder]

        '''Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff'''
        N_InternalNodes = len(VirusDendrogram.get_nonterminals())
        for InternalNode_i in range(N_InternalNodes):
            try:
                '''Here we want to ensure there are no labels on the dendrogram where bootstrap confidence < cutoff (or not present if BS = False)'''
                if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence < self.payload['Heatmap_DendrogramSupport_Cutoff'] or np.isnan(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence):
                    continue
                else:
                    VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(
                        VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
            except:
                '''Exception as first entry will always be None'''
                continue

        '''Colour terminal branches: reference virus's branch is blue, unclassified virus's branch is red'''
        N_Viruses = N_RefViruses + self.N_UcfViruses
        for Virus_i in range(N_Viruses):
            Taxolabel = VirusDendrogram.get_terminals()[Virus_i].name
            if Taxolabel in TaxoLabelList_RefVirus:
                VirusDendrogram.get_terminals()[Virus_i].color = "blue"
            elif Taxolabel in self.TaxoLabelList_UcfVirus:
                VirusDendrogram.get_terminals()[Virus_i].color = "red"

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

        '''Heat map colour indicators'''
        IndicatorMat_RefVirus = np.tile(
            [True]*N_RefViruses + [False]*self.N_UcfViruses, N_Viruses).reshape(N_Viruses, -1)
        IndicatorMat_RefVirus = IndicatorMat_RefVirus * IndicatorMat_RefVirus.T
        IndicatorMat_RefVirus = IndicatorMat_RefVirus[VirusOrder][:, VirusOrder]
        IndicatorMat_UcfVirus = np.tile(
            [False]*N_RefViruses + [True]*self.N_UcfViruses, N_Viruses).reshape(N_Viruses, -1)
        IndicatorMat_UcfVirus = IndicatorMat_UcfVirus * IndicatorMat_UcfVirus.T
        IndicatorMat_UcfVirus = IndicatorMat_UcfVirus[VirusOrder][:, VirusOrder]
        IndicatorMat_CrossGroup = ~IndicatorMat_RefVirus * ~IndicatorMat_UcfVirus

        '''Masked OrderedDistMat'''
        OrderedDistMat_RefVirus = np.ma.masked_where(
            ~IndicatorMat_RefVirus, OrderedDistMat)
        OrderedDistMat_UcfVirus = np.ma.masked_where(
            ~IndicatorMat_UcfVirus, OrderedDistMat)
        OrderedDistMat_CrossGroup = np.ma.masked_where(
            ~IndicatorMat_CrossGroup, OrderedDistMat)

        '''Colour map construction'''
        MyBlues = get_blue_cmap()
        MyReds = get_red_cmap()
        MyPurples = get_purple_cmap()

        '''Plot configuration'''
        hmap_params, fig, ax_Dendrogram, ax_Heatmap = get_hmap_params(len(ClassLabelList_x))

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
        ax_Heatmap		.imshow(OrderedDistMat_RefVirus, cmap=MyBlues,
                            aspect='auto', vmin=0, vmax=1, interpolation='none')
        ax_Heatmap		.imshow(OrderedDistMat_UcfVirus, cmap=MyReds,
                            aspect='auto', vmin=0, vmax=1, interpolation='none')
        ax_Heatmap		.imshow(OrderedDistMat_CrossGroup, cmap=MyPurples,
                            aspect='auto', vmin=0, vmax=1, interpolation='none')

        '''Draw grouping major & minor lines'''
        ax_Heatmap = construct_hmap_lines(ax_Heatmap, LineList_major, LineList_minor,
                                          hmap_params, ClassLabelList_x, ClassLabelList_y,
                                          TickLocList = np.array(
                                            list(map(np.mean, list(zip(LineList_minor[0:-1], LineList_minor[1:]))))))

        '''Selectively colour tick labels red if a UCF sample'''
        [i.set_color("red") for i in ax_Heatmap.get_xticklabels()
         if bool(re.match(r"Query", i.get_text()))]
        [i.set_color("red") for i in ax_Heatmap.get_yticklabels()
         if bool(re.match(r"Query", i.get_text()))]

        '''Reference virus colour bar'''
        ax_CBar_RefVirus = fig.add_axes(
            [hmap_params['ax_CBar_L'], hmap_params['ax_CBar_B'] + 2*hmap_params['ax_CBar_H']/3, hmap_params['ax_CBar_W'], hmap_params['ax_CBar_H']/3], frame_on=True, facecolor="white")
        ax_CBar_RefVirus	.imshow(np.linspace(0, 1, 1025).reshape(
            1, -1), cmap=MyBlues, aspect='auto', vmin=0, vmax=1, interpolation='none')
        ax_CBar_RefVirus	.set_yticks([0.0])
        ax_CBar_RefVirus	.set_yticklabels(
            ["Ref viruses"], rotation=0, size=hmap_params['FontSize'] + 2)
        ax_CBar_RefVirus	.tick_params(top=False,
                                      bottom=False,
                                      left=False,
                                      right=True,
                                      labeltop=False,
                                      labelbottom=False,
                                      labelleft=False,
                                      labelright=True,
                                      direction='out')

        '''Unclassified virus colour bar'''
        ax_CBar_UcfVirus = fig.add_axes(
            [hmap_params['ax_CBar_L'], hmap_params['ax_CBar_B'] + 1*hmap_params['ax_CBar_H']/3, hmap_params['ax_CBar_W'], hmap_params['ax_CBar_H']/3], frame_on=True, facecolor="white")
        ax_CBar_UcfVirus	.imshow(np.linspace(0, 1, 1025).reshape(
            1, -1), cmap=MyReds, aspect='auto', vmin=0, vmax=1, interpolation='none')
        ax_CBar_UcfVirus	.set_yticks([0.0])
        ax_CBar_UcfVirus	.set_yticklabels(
            ["Ucf viruses"], rotation=0, size=hmap_params['FontSize'] + 2)
        ax_CBar_UcfVirus	.tick_params(top=False,
                                      bottom=False,
                                      left=False,
                                      right=True,
                                      labeltop=False,
                                      labelbottom=False,
                                      labelleft=False,
                                      labelright=True,
                                      direction='out')

        '''Ref/Ucf merge colour bar'''
        ax_CBar_CrossGroup = fig.add_axes(
            [hmap_params['ax_CBar_L'], hmap_params['ax_CBar_B'] + 0*hmap_params['ax_CBar_H']/3, hmap_params['ax_CBar_W'], hmap_params['ax_CBar_H']/3], frame_on=True, facecolor="white")
        ax_CBar_CrossGroup	.imshow(np.linspace(0, 1, 1025).reshape(
            1, -1), cmap=MyPurples, aspect='auto', vmin=0, vmax=1, interpolation='none')
        ax_CBar_CrossGroup	.set_yticks([0.0])
        ax_CBar_CrossGroup	.set_yticklabels(
            ["Ref VS Ucf viruses"], rotation=0, size=hmap_params['FontSize'] + 2)
        ax_CBar_CrossGroup	.set_xticks(
            np.array([0, 0.25, 0.50, 0.75, 1])*(1025)-0.5)
        ax_CBar_CrossGroup	.set_xticklabels(
            ['0', '0.25', '0.50', '0.75', '1'], rotation=0, size=hmap_params['FontSize'])
        ax_CBar_CrossGroup	.set_xlabel(
            "Distance", rotation=0, size=hmap_params['FontSize'] + 2)
        ax_CBar_CrossGroup	.tick_params(top=False,
                                        bottom=True,
                                        left=False,
                                        right=True,
                                        labeltop=False,
                                        labelbottom=True,
                                        labelleft=False,
                                        labelright=True,
                                        direction='out')
        plt.savefig(self.fnames['HeatmapWithDendrogramFile'], format="pdf", bbox_inches = "tight", dpi=hmap_params['dpi'])

    def group(self):
        '''7/8 : Pool results from all classfiers, and finalise the taxonomic assignment (and virus grouping)'''
        progress_msg("Pooling results from classifiers, making final taxonomic groupings")
        if self.payload['VirusGrouping']:
            FinalisedTaxoAssignmentList, FinalisedVirusGroupingList, FinalisedVirusGroupingDict, FinalisedVirusGroupingIndex = [], [], {}, 1
            for UcfVirus_i in range(self.N_UcfViruses):
                MaxSimScore_i = np.argmax(
                    self.final_results["MaxSimScoreTable"][UcfVirus_i])
                FinalisedTaxoAssignment = self.final_results["TaxoAssignmentTable"][UcfVirus_i][MaxSimScore_i]

                if "Unclassified" not in FinalisedTaxoAssignment:
                    FinalisedTaxoAssignmentList.append(
                        f"{FinalisedTaxoAssignment}")
                    FinalisedVirusGroupingList.append(
                        f"{FinalisedTaxoAssignment}")
                else:
                    if self.final_results["MaxSimScoreTable"][UcfVirus_i][MaxSimScore_i] <= self.payload['DatabaseAssignmentSimilarityScore_Cutoff']:
                        FinalisedTaxoAssignmentList.append("Unclassified (NA)")
                        FinalisedVirusGroupingList.append("Unclassified (NA)")
                    else:
                        FinalisedTaxoAssignmentList.append(
                            f"Unclassified")
                        VirusGrouping = self.final_results["VirusGrouping"][
                            "VirusGroupingList"][-self.N_UcfViruses:][UcfVirus_i]
                        '''Novel taxa are labelled unassigned taxonomy units (UTUs)'''
                        if VirusGrouping not in FinalisedVirusGroupingDict.keys():
                            FinalisedVirusGroupingDict[
                                VirusGrouping] = f"UTU{FinalisedVirusGroupingIndex}"
                            FinalisedVirusGroupingIndex += 1

                        FinalisedVirusGroupingList.append(
                            f"{FinalisedVirusGroupingDict[VirusGrouping]}")

        else:
            FinalisedTaxoAssignmentList, FinalisedVirusGroupingList = [], [""]*self.N_UcfViruses
            for UcfVirus_i in range(self.N_UcfViruses):
                MaxSimScore_i = np.argmax(
                    self.final_results["MaxSimScoreTable"][UcfVirus_i])
                FinalisedTaxoAssignment = self.final_results["TaxoAssignmentTable"][UcfVirus_i][MaxSimScore_i]

                if "Unclassified" not in FinalisedTaxoAssignment:
                    FinalisedTaxoAssignmentList.append(
                        f"{FinalisedTaxoAssignment}")
                else:
                    if self.final_results["MaxSimScoreTable"][UcfVirus_i][MaxSimScore_i] <= self.payload['DatabaseAssignmentSimilarityScore_Cutoff']:
                        FinalisedTaxoAssignmentList.append("Unclassified (NA)")
                    else:
                        FinalisedTaxoAssignmentList.append(
                            f"Unclassified")

        RemainedUcfVirus_IndexList = np.where(list(map(
            all, self.final_results["MaxSimScoreTable"] <= self.payload['DatabaseAssignmentSimilarityScore_Cutoff'])))[0]
        self.final_results["FinalisedTaxoAssignmentList"] = FinalisedTaxoAssignmentList
        self.final_results["FinalisedVirusGroupingList"] = FinalisedVirusGroupingList
        self.final_results["SeqIDLists_RemainedUcfVirus"] = self.genomes["SeqIDLists"][RemainedUcfVirus_IndexList]
        self.final_results["VirusNameList_RemainedUcfVirus"] = self.genomes["VirusNameList"][RemainedUcfVirus_IndexList]

        if self.payload['Bootstrap']:
            TaxoOfMaxSimScoreRangeTable, MaxSimScoreRangeTable, PhyloStatRangeTable, \
                TaxoAssignmentRangeTable = np.zeros((self.N_UcfViruses, 0)), np.zeros(
                    (self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0))

            TaxoOfMaxSimScoreRangeList, MaxSimScoreRangeList, PhyloStatRangeList, TaxoAssignmentRangeList = [], [], [], []
            for UcfVirus_i in range(self.N_UcfViruses):
                TaxoOfMaxSimScoreDist_Counter = sorted(list(Counter(
                    self.final_results["TaxoOfMaxSimScoreDistTable"][:, UcfVirus_i]).items()), key=operator.itemgetter(1), reverse=True)
                TaxoOfMaxSimScoreRangeList.append(", ".join(
                    [f"{Taxo}: {round(Count/float(self.payload['N_Bootstrap']))}" for Taxo, Count in TaxoOfMaxSimScoreDist_Counter]))
                MaxSimScoreRangeList.append(", ".join(["%s: %s" % (Taxo, "-".join(map(str, np.around(hpd(x=self.final_results["MaxSimScoreDistTable"][np.where(
                    self.final_results["TaxoOfMaxSimScoreDistTable"][:, UcfVirus_i] == Taxo)[0], UcfVirus_i], alpha=0.05), 3)))) for Taxo in list(zip(*TaxoOfMaxSimScoreDist_Counter))[0]]))
                PhyloStatRangeList.append(", ".join(["%s: %s" % (Taxo, ",".join(["p(%s): %s" % (x[0], x[1]/float(self.payload['N_Bootstrap'])) for x in sorted(list(Counter(self.final_results["PhyloStatDistTable"][np.where(
                    self.final_results["TaxoOfMaxSimScoreDistTable"][:, UcfVirus_i] == Taxo)[0], UcfVirus_i]).items()), key=operator.itemgetter(1), reverse=True)])) for Taxo in list(zip(*TaxoOfMaxSimScoreDist_Counter))[0]]))
                TaxoAssignmentRangeList.append(", ".join(["%s: %.2f" % (Taxo, Count/float(self.payload['N_Bootstrap'])) for Taxo, Count in sorted(list(
                    Counter(self.final_results["TaxoAssignmentDistTable"][:, UcfVirus_i]).items()), key=operator.itemgetter(1), reverse=True)]))

            TaxoOfMaxSimScoreRangeTable = np.column_stack(
                (TaxoOfMaxSimScoreRangeTable, TaxoOfMaxSimScoreRangeList))
            MaxSimScoreRangeTable = np.column_stack(
                (MaxSimScoreRangeTable, MaxSimScoreRangeList))
            PhyloStatRangeTable = np.column_stack(
                (PhyloStatRangeTable, PhyloStatRangeList))
            TaxoAssignmentRangeTable = np.column_stack(
                (TaxoAssignmentRangeTable, TaxoAssignmentRangeList))
            TaxoOfMaxSimScoreRangeTable = np.column_stack(
                (TaxoOfMaxSimScoreRangeTable, TaxoOfMaxSimScoreRangeList))
            MaxSimScoreRangeTable = np.column_stack(
                (MaxSimScoreRangeTable, MaxSimScoreRangeList))
            PhyloStatRangeTable = np.column_stack(
                (PhyloStatRangeTable, PhyloStatRangeList))
            TaxoAssignmentRangeTable = np.column_stack(
                (TaxoAssignmentRangeTable, TaxoAssignmentRangeList))

            '''Save final results from iteratively updated tables'''
            self.final_results["TaxoOfMaxSimScoreRangeTable"] = TaxoOfMaxSimScoreRangeTable
            self.final_results["MaxSimScoreRangeTable"] = MaxSimScoreRangeTable
            self.final_results["PhyloStatRangeTable"] = PhyloStatRangeTable
            self.final_results["TaxoAssignmentRangeTable"] = TaxoAssignmentRangeTable

            FinalisedTaxoAssignmentRangeList = []
            for UcfVirus_i in range(self.N_UcfViruses):
                BootstrappedFinalisedTaxoAssignmentList = []
                for Bootstrap_i in range(0, self.payload['N_Bootstrap']):
                    MaxSimScore_i = np.argmax(
                        self.final_results["MaxSimScoreDistTable"][Bootstrap_i][UcfVirus_i])
                    FinalisedTaxoAssignment = self.final_results["TaxoAssignmentDistTable"][
                        Bootstrap_i][UcfVirus_i]
                    if "Unclassified" not in FinalisedTaxoAssignment:
                        BootstrappedFinalisedTaxoAssignmentList.append(
                            f"{FinalisedTaxoAssignment}")
                    else:
                        if self.final_results["MaxSimScoreDistTable"][Bootstrap_i][UcfVirus_i] <= self.payload['DatabaseAssignmentSimilarityScore_Cutoff']:
                            BootstrappedFinalisedTaxoAssignmentList.append(
                                "Unclassified (NA)")
                        else:
                            BootstrappedFinalisedTaxoAssignmentList.append(
                                f"Unclassified")
                FinalisedTaxoAssignmentRangeList.append(", ".join(["%s: %.2f" % (Taxo, Count/float(self.payload['N_Bootstrap'])) for Taxo, Count in sorted(
                    list(Counter(BootstrappedFinalisedTaxoAssignmentList).items()), key=operator.itemgetter(1), reverse=True)]))
            self.final_results["FinalisedTaxoAssignmentRangeList"] = FinalisedTaxoAssignmentRangeList

    def write(self):
        '''8/8: Write results to persistent storage'''
        progress_msg("Writing results to file")

        '''Prepare common fields of output CSV'''
        closest_taxa = pd.DataFrame()
        closest_taxa["isolate_name"] = self.genomes["VirusNameList"].copy()
        closest_taxa["run_designation"] = list(
            map(', '.join, self.genomes["SeqIDLists"]))

        '''Get Cutoff scores for each grouping, for printing in TXT output'''
        sim_cutoff_string = " ".join([f'Grouping: {i}, Cutoff: {self.final_results["PairwiseSimilarityScore_Cutoff_Dict"][i]["CutOff"]}\n' for i in self.final_results["PairwiseSimilarityScore_Cutoff_Dict"].keys()])

        '''Prepare output TXT and CSV depending on bootstrap method'''
        if self.payload['Bootstrap']:
            np.savetxt(fname=self.fnames['ClassificationResultFile'],
                       X=np.column_stack((list(map(', '.join, self.genomes["SeqIDLists"])),
                                          self.genomes["VirusNameList"],
                                          [", ".join(["%s (%s)" % (Taxo, Range) for Taxo, Range in zip(TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange[0], TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange[1])])
                                           for TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange in zip(self.final_results["TaxoOfMaxSimScoreTable"], self.final_results["TaxoOfMaxSimScoreRangeTable"])],
                                          [", ".join(["%s (%s)" % (Score, Range) for Score, Range in zip(MaxSimScore_MaxSimScoreRange[0], MaxSimScore_MaxSimScoreRange[1])]) for MaxSimScore_MaxSimScoreRange in zip(
                                              np.around(self.final_results["MaxSimScoreTable"], 3).astype("str"), self.final_results["MaxSimScoreRangeTable"])],
                                          [", ".join(["%s (%s)" % (Stat, Range) for Stat, Range in zip(PhyloStat_PhyloStatRange[0], PhyloStat_PhyloStatRange[1])])
                                           for PhyloStat_PhyloStatRange in zip(self.final_results["PhyloStatTable"], self.final_results["PhyloStatRangeTable"])],
                                          [", ".join(["%s (%s)" % (Taxo, Range) for Taxo, Range in zip(TaxoAssignment_TaxoAssignmentRange[0], TaxoAssignment_TaxoAssignmentRange[1])])
                                           for TaxoAssignment_TaxoAssignmentRange in zip(self.final_results["TaxoAssignmentTable"], self.final_results["TaxoAssignmentRangeTable"])],
                                          ["%s (%s)" % (FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange[0], FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange[1]) for FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange in zip(
                                              self.final_results["FinalisedTaxoAssignmentList"], self.final_results["FinalisedTaxoAssignmentRangeList"])],
                                          self.final_results["FinalisedVirusGroupingList"],
                                          )),
                       fmt='%s',
                       delimiter="\t",
                       header="Sequence identifier\tSequence description\tCandidate class (class of the best match reference virus)\tSimilarity score*\tSupport from dendrogram**\tEvaluated taxonomic assignment\tThe best taxonomic assignment\tProvisional virus taxonomy")

            '''Update grouping CSV'''
            closest_taxa["closest_taxon"] = [", ".join(["%s (%s)" % (Taxo, Range) for Taxo, Range in zip(TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange[0], TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange[1])])
                                             for TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange in zip(self.final_results["TaxoOfMaxSimScoreTable"], self.final_results["TaxoOfMaxSimScoreRangeTable"])]
            closest_taxa["score"] = [", ".join(["%s (%s)" % (Score, Range) for Score, Range in zip(MaxSimScore_MaxSimScoreRange[0], MaxSimScore_MaxSimScoreRange[1])]) for MaxSimScore_MaxSimScoreRange in zip(
                np.around(self.final_results["MaxSimScoreTable"], 3).astype("str"), self.final_results["MaxSimScoreRangeTable"])]
            closest_taxa["phylo_stat"] = [", ".join(["%s (%s)" % (Stat, Range) for Stat, Range in zip(PhyloStat_PhyloStatRange[0], PhyloStat_PhyloStatRange[1])])
                                          for PhyloStat_PhyloStatRange in zip(self.final_results["PhyloStatTable"], self.final_results["PhyloStatRangeTable"])]
            closest_taxa["taxo_assignment"] = ["%s (%s)" % (FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange[0], FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange[1]) for FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange in zip(
                self.final_results["FinalisedTaxoAssignmentList"], self.final_results["FinalisedTaxoAssignmentRangeList"])]
            closest_taxa["taxo_grouping"] = self.final_results["FinalisedVirusGroupingList"]

        else:
            np.savetxt(fname=self.fnames['ClassificationResultFile'],
                       X=np.column_stack((list(map(', '.join, self.genomes["SeqIDLists"])),
                                          self.genomes["VirusNameList"],
                                          list(
                                              map(', '.join, self.final_results["TaxoOfMaxSimScoreTable"])),
                                          list(map(', '.join, np.around(
                                              self.final_results["MaxSimScoreTable"], 3).astype("str"))),
                                          list(
                                              map(', '.join, self.final_results["PhyloStatTable"])),
                                          list(
                                              map(', '.join, self.final_results["TaxoAssignmentTable"])),
                                          self.final_results["FinalisedTaxoAssignmentList"],
                                          self.final_results["FinalisedVirusGroupingList"],
                                          )),
                       fmt='%s',
                       delimiter="\t",
                       header="Sequence identifier\tSequence description\tCandidate class (class of the best match reference virus)\tSimilarity score*\tSupport from dendrogram**\tEvaluated taxonomic assignment\tThe best taxonomic assignment\tProvisional virus taxonomy")

            '''Update grouping CSV'''
            closest_taxa["closest_taxon"] = list(
                map(', '.join, self.final_results["TaxoOfMaxSimScoreTable"]))
            closest_taxa["score"] = list(map(', '.join, np.around(
                self.final_results["MaxSimScoreTable"], 3).astype("str")))
            closest_taxa["phylo_stat"] = list(
                map(', '.join, self.final_results["PhyloStatTable"]))
            closest_taxa["taxo_assignment"] = self.final_results["FinalisedTaxoAssignmentList"]
            closest_taxa["taxo_grouping"] = self.final_results["FinalisedVirusGroupingList"]

        '''Save output CSV, TXT and Pickle files with classification results'''
        with open(self.fnames['ClassificationResultFile'], "a") as ClassificationResult_txt:
            ClassificationResult_txt.write(f"\n"
                                            f"*Similarity score cutoff\n"
                                            f"{sim_cutoff_string}\n"
                                            f"\n\n"
                                            f"**Support from dendrogram\n"
                                            f"NA:the sequence is not similar enough to any of the reference sequences to be assigned to any classes\n"
                                            f"1: the sequence is embedded within a clade of the candidate class\n"
                                            f"2: the sequence has a sister relationship with the candidate class and they are similar enough (the 1st criterion)\n"
                                            f"3: the sequence is 'sandwiched' between 2 branches of the candidate class\n"
                                            f"4: the sequence has a paraphyletic relationship with the candidate class (just inside)\n"
                                            f"5: the sequence has a paraphyletic relationship with the candidate class (just outside)\n"
                                            f"6: the candidate class is not supported by the dendrogram\n"
                                            )
        closest_taxa.to_csv(self.fnames['ClassificationResultFile'].replace(".txt", ".csv"))
        pickle.dump(self.final_results, open(self.fnames['VirusClassifierPickle'], "wb"))

    def main(self):
        '''Classify viruses and evaluate results'''
        section_header("Classify viruses and evaluate the results")
        error_handler_virus_classifier(self.payload)

        '''2/8: Make fpaths'''
        mkdir_virus_classifier(self.fnames)

        if self.payload['UseUcfVirusPPHMMs']:
            '''(OPT) 4/8: Scan unclassified viruses against PPHMMDB, create additional sig and loc table'''
            self.use_ucf_virus_pphmms()

        '''6/8: Classify viruses'''
        self.classify()

        '''7/8: Finalise taxonomic assignment and grouping'''
        self.group()

        '''8/8: Write results'''
        self.write()
