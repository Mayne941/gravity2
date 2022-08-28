from Bio import Phylo
from copy import copy
from ete3 import Tree
from sklearn.svm import SVC
from scipy.cluster.hierarchy import linkage, fcluster
from collections import Counter
import matplotlib
matplotlib.use('agg')
import numpy as np
from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import os, shelve, sys, random, string, subprocess, operator, pickle
sys.setrecursionlimit(10000)

#Local functions
from app.utils.ordered_set import OrderedSet
from app.utils.dist_mat_to_tree import DistMat2Tree
from app.utils.similarity_matrix_constructor import SimilarityMat_Constructor
from app.utils.taxo_label_constructor import TaxoLabel_Constructor
from app.utils.pphmm_signature_table_constructor import PPHMMSignatureTable_Constructor
from app.utils.gomdb_constructor import GOMDB_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.virus_grouping_estimator import VirusGrouping_Estimator
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_ucf_annots, retrieve_ref_virus_vars, retrieve_pphmmdb_construction
from app.utils.heatmap_with_dendrogram import heatmap_with_dendrogram
from app.utils.highest_posterior_density import hpd

def PairwiseSimilarityScore_Cutoff_Dict_Constructor (SimMat, TaxoGroupingList, N_PairwiseSimilarityScores):
	N_Viruses = len(SimMat)
	PairwiseSimilarityScore_Cutoff_Dict = {TaxoGrouping:{"PairwiseSimilarityScore_InterClass_List":[],"PairwiseSimilarityScore_IntraClass_List":[]} for TaxoGrouping in OrderedSet(TaxoGroupingList)}
	for Virus_i in range(N_Viruses):
		for Virus_j in range(Virus_i, N_Viruses):
			Class_Virus_i = TaxoGroupingList[Virus_i]
			Class_Virus_j = TaxoGroupingList[Virus_j]
			if Class_Virus_i == Class_Virus_j:
				PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_i]["PairwiseSimilarityScore_IntraClass_List"].append(SimMat[Virus_i][Virus_j])
			else:
				PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_i]["PairwiseSimilarityScore_InterClass_List"].append(SimMat[Virus_i][Virus_j])
				PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_j]["PairwiseSimilarityScore_InterClass_List"].append(SimMat[Virus_i][Virus_j])
	
	for TaxoGrouping in OrderedSet(TaxoGroupingList):
		PairwiseSimilarityScore_IntraClass_List = PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_IntraClass_List"]
		PairwiseSimilarityScore_InterClass_List = PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_InterClass_List"]
		
		if len(PairwiseSimilarityScore_IntraClass_List) > N_PairwiseSimilarityScores:
			PairwiseSimilarityScore_IntraClass_List = np.random.choice(PairwiseSimilarityScore_IntraClass_List, N_PairwiseSimilarityScores, replace = False).tolist()
			PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_IntraClass_List"] = PairwiseSimilarityScore_IntraClass_List
		
		if len(PairwiseSimilarityScore_InterClass_List) > N_PairwiseSimilarityScores:
			PairwiseSimilarityScore_InterClass_List = np.random.choice(PairwiseSimilarityScore_InterClass_List, N_PairwiseSimilarityScores, replace = False).tolist()
			PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_InterClass_List"] = PairwiseSimilarityScore_InterClass_List
		
		PairwiseSimilarityList	= np.array(PairwiseSimilarityScore_IntraClass_List+PairwiseSimilarityScore_InterClass_List).reshape(-1,1)
		Labels			= np.array(["Intra"]*len(PairwiseSimilarityScore_IntraClass_List)+["Inter"]*len(PairwiseSimilarityScore_InterClass_List))
		
		clf			= SVC(kernel = "linear", class_weight = "balanced")
		clf			.fit(PairwiseSimilarityList, Labels)
		
		PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["CutOff"] = -clf.intercept_[0]/clf.coef_[0][0]
	
	return PairwiseSimilarityScore_Cutoff_Dict

def TaxonomicAssignmentProposerAndEvaluator (SimMat_UcfVirusesVSRefViruses, TaxoGrouping, VirusDendrogram, TaxoLabelList_RefVirus, TaxoLabelList_UcfVirus, PairwiseSimilarityScore_Cutoff_Dict):
	print("\t\tPropose a taxonomic class to each unclassified viruses, using the 1-nearest neighbour algorithm")
	#-------------------------------------------------------------------------------
	MaxSimScoreList			= np.max(SimMat_UcfVirusesVSRefViruses, axis = 1)
	TaxoOfMaxSimScoreList		= TaxoGrouping[np.argmax(SimMat_UcfVirusesVSRefViruses, axis = 1)]
	
	print("\t\tEvaluate the taxonomic assignments")
	#-------------------------------------------------------------------------------
	#Rename the reference virus leaves in the dendrogram to class label
	#-------------------------------------------------------------------------------
	VirusDendrogram = Tree(VirusDendrogram)
	for LeafNode in VirusDendrogram.get_leaves():
		if LeafNode.name in TaxoLabelList_RefVirus:
			LeafNode.name = TaxoGrouping[np.where(np.array(TaxoLabelList_RefVirus) == LeafNode.name)[0]][0]
	
	TaxoAssignmentList	= []
	PhyloStatList		= []
	N_UcfViruses		= len(SimMat_UcfVirusesVSRefViruses)
	for UcfVirus_i in range(N_UcfViruses):
		CandidateTaxoAssignment = TaxoOfMaxSimScoreList[UcfVirus_i]
		MaxSimScore = MaxSimScoreList[UcfVirus_i]
		
		#1st criterion: see if the similarity score between the unclassified virus and the best match is greater than the similarity cutoff of the proposed class
		#-------------------------------------------------------------------------------
		if MaxSimScore > PairwiseSimilarityScore_Cutoff_Dict[CandidateTaxoAssignment]["CutOff"]:	#if pass the 1st criterion ...
			#Prune the dendrogram so that it contains only reference sequences and the virus of interest (VoI)
			#-------------------------------------------------------------------------------
			VirusDendrogram_tmp	= VirusDendrogram.copy()
			for q in np.delete(TaxoLabelList_UcfVirus, UcfVirus_i):
				VirusDendrogram_tmp.get_leaves_by_name(q)[0].delete()
			
			#Get the leaf names (taxonomic assignments) of the sister clade of the unclassified virus
			#-------------------------------------------------------------------------------
			VoI			= TaxoLabelList_UcfVirus[UcfVirus_i]
			VoISisterClade		= VirusDendrogram_tmp.search_nodes(name = VoI)[0].get_sisters()[0]
			VoISisterClade_Leafnames= OrderedSet(VoISisterClade.get_leaf_names())
			if VoISisterClade.is_leaf():
				VoISisterClade_Leafnames1 = OrderedSet(VoISisterClade.get_leaf_names())
				VoISisterClade_Leafnames2 = []
			else:
				VoISisterClade_Leafnames1 = OrderedSet(VoISisterClade.get_children()[0].get_leaf_names())
				VoISisterClade_Leafnames2 = OrderedSet(VoISisterClade.get_children()[1].get_leaf_names())
			
			#Get the leaf names (taxonomic assignments) of the immediate out group of the unclassified virus
			#-------------------------------------------------------------------------------
			ImmediateAncestorNode = VirusDendrogram_tmp.search_nodes(name = VoI)[0].up
			#if not ImmediateAncestorNode.is_root():	#that is, if the branch isn't the most outer branch itself...
			if len(ImmediateAncestorNode.get_sisters())!=0: #that is, if the branch isn't the most outer branch itself...
				ImmediateOutCladeLeafnames = OrderedSet(ImmediateAncestorNode.get_sisters()[0].get_leaf_names())
			else:
				ImmediateOutCladeLeafnames = []
			
			#2nd criterion: see if the candidate taxonomic assignment is supported by the dendrogram...
			#-------------------------------------------------------------------------------
			if len(VoISisterClade_Leafnames) == 1\
			and [CandidateTaxoAssignment] == VoISisterClade_Leafnames:		#If the sister clade consists entirely of the candidate class...
				if VoISisterClade_Leafnames == ImmediateOutCladeLeafnames:	#... and the out group is the same, then assign VoI to the candidate class
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("1") #Embedded within a clade of a single class
				else:							#... but the out group isn't...
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("2") #Having a sister relationship with the proposed class and they are similar enough (the 1st criterion)
			elif len(ImmediateOutCladeLeafnames) == 1\
			and [CandidateTaxoAssignment] == ImmediateOutCladeLeafnames:	#If the immediate outgroup consists entirely of the candidate class...
				if(VoISisterClade_Leafnames1 == [CandidateTaxoAssignment] or\
				VoISisterClade_Leafnames2 == [CandidateTaxoAssignment]):		#... and the VoI is sandwished between 2 branches of the candidate class... 
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("3") #Sandwished between 2 branches of the candidate class
				else:							#... the candidate class is accepted on the ground that it has a paraphyletic relationship with the candidate class (just inside)
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("4") #Having a paraphyletic relationship with the candidate class (just inside)
			elif (VoISisterClade_Leafnames1 == [CandidateTaxoAssignment] or\
				VoISisterClade_Leafnames2 == [CandidateTaxoAssignment]):	#If one of the two branches in the sister clade consists entirely of the candidate class...
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("5") #Having a paraphyletic relationship with the candidate class (just outside)
			else:
				TaxoAssignmentList.append("Unclassified")
				PhyloStatList.append("6") #The candidate class is not supported by the dendrogram
		else: #The unclassified virus isn't similar enough to the members of the candidate class
			TaxoAssignmentList.append("Unclassified")
			PhyloStatList.append("NA")
		
		#Progress bar
		sys.stdout.write("\033[K" + "Taxonomic assignment evaluation: [%-20s] %d/%d viruses" % ('='*int(float(UcfVirus_i+1)/N_UcfViruses*20), UcfVirus_i+1, N_UcfViruses) + "\r")
		sys.stdout.flush()
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	
	return MaxSimScoreList, TaxoOfMaxSimScoreList, TaxoAssignmentList, PhyloStatList

class VirusClassificationAndEvaluation:
	def __init__(self,
				ShelveDir_UcfVirus,
				ShelveDirs_RefVirus,
				IncludeIncompleteGenomes_UcfVirus	= True,
				IncludeIncompleteGenomes_RefVirus	= False,
				UseUcfVirusPPHMMs				= True,
				GenomeSeqFile_UcfVirus			= None,
				GenomeSeqFiles_RefVirus			= None,
				SeqLength_Cutoff				= 0,
				HMMER_N_CPUs					= 20,
				HMMER_C_EValue_Cutoff			= 1E-3,
				HMMER_HitScore_Cutoff			= 0,
				SimilarityMeasurementScheme		= "PG",
				p								= 1,
				Dendrogram_LinkageMethod		= "average", #'single', 'complete', 'average', 'weighted'
				DatabaseAssignmentSimilarityScore_Cutoff= 0.01,
				N_PairwiseSimilarityScores		= 10000,
				Heatmap_WithDendrogram			= True,
				Heatmap_DendrogramSupport_Cutoff= 0.75,
				Bootstrap						= True,
				N_Bootstrap						= 10,
				Bootstrap_method				= "booster", #"sumtrees", "booster"
				Bootstrap_N_CPUs				= 20,
				VirusGrouping					= True,
			):
		'''Fpaths'''
		self.ShelveDir_UcfVirus  = ShelveDir_UcfVirus
		self.ShelveDirs_RefVirus = ShelveDirs_RefVirus.split(", ")
		self.VariableShelveDir_UcfVirus			= f"{ShelveDir_UcfVirus}/Shelves"
		self.HMMERDir_UcfVirus = self.HMMER_PPHMMDB_UcfVirus = "placeholder"
		'''Parameters'''
		self.IncludeIncompleteGenomes_UcfVirus	= IncludeIncompleteGenomes_UcfVirus
		self.IncludeIncompleteGenomes_RefVirus	= IncludeIncompleteGenomes_RefVirus
		self.UseUcfVirusPPHMMs				= UseUcfVirusPPHMMs
		self.GenomeSeqFile_UcfVirus			= GenomeSeqFile_UcfVirus
		self.GenomeSeqFiles_RefVirus		= GenomeSeqFiles_RefVirus
		self.SeqLength_Cutoff				= SeqLength_Cutoff
		self.HMMER_N_CPUs					= HMMER_N_CPUs
		self.HMMER_C_EValue_Cutoff			= HMMER_C_EValue_Cutoff
		self.HMMER_HitScore_Cutoff			= HMMER_HitScore_Cutoff
		self.SimilarityMeasurementScheme	= SimilarityMeasurementScheme
		self.p								= p	
		self.Dendrogram_LinkageMethod		= Dendrogram_LinkageMethod
		self.DatabaseAssignmentSimilarityScore_Cutoff= DatabaseAssignmentSimilarityScore_Cutoff
		self.N_PairwiseSimilarityScores		= N_PairwiseSimilarityScores
		self.Heatmap_WithDendrogram			= Heatmap_WithDendrogram
		self.Heatmap_DendrogramSupport_Cutoff= Heatmap_DendrogramSupport_Cutoff
		self.Bootstrap						= Bootstrap
		self.N_Bootstrap					= N_Bootstrap
		self.Bootstrap_method				= Bootstrap_method
		self.Bootstrap_N_CPUs				= Bootstrap_N_CPUs
		self.VirusGrouping					= VirusGrouping
		'''Placeholder arrays'''
		self.ucf_genomes, self.PairwiseSimilarityScore_Cutoff_Dict, self.VirusGroupingDict = {}, {}, {}
		self.TaxoLabelList_UcfVirus, self.N_UcfViruses, self.N_RefVirusGroups = [], [], []
		self.MaxSimScoreTable, self.TaxoOfMaxSimScoreTable, self.TaxoAssignmentTable, self.PhyloStatTable = np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0))
		self.PairwiseSimilarityScore_CutoffDist_Dict, self.MaxSimScoreDistTable, \
			self.TaxoOfMaxSimScoreDistTable, self.TaxoAssignmentDistTable, self.PhyloStatDistTable = None, None, None, None, None
		
	def input_test(self):
		'''1/ : Guide user to specify input PPHMM file'''
		if self.GenomeSeqFiles_RefVirus == None or self.GenomeSeqFile_UcfVirus == None:
			raise ValueError("'UseUcfVirusPPHMMs' == True. 'GenomeSeqFiles_RefVirus' and 'GenomeSeqFile_UcfVirus' cannot be None")
		
		self.GenomeSeqFiles_RefVirus = self.GenomeSeqFiles_RefVirus.split(", ")
		if len(self.ShelveDirs_RefVirus) != len(self.GenomeSeqFiles_RefVirus):
			raise ValueError("Length of 'ShelveDirs_RefVirus' should be equal to that of 'GenomeSeqFiles_RefVirus'")

	def mkdirs(self):
		'''2/ : Return all dirs for db storage and retrieval'''
		if self.UseUcfVirusPPHMMs == True:
			self.HMMERDir_UcfVirus		= f"{self.ShelveDir_UcfVirus}/HMMER"
			HMMER_PPHMMDBDir_UcfVirus	= f"{self.HMMERDir_UcfVirus}/HMMER_PPHMMDB"
			self.HMMER_PPHMMDB_UcfVirus	= f"{HMMER_PPHMMDBDir_UcfVirus}/HMMER_PPHMMDB"
			if not os.path.exists(self.HMMER_PPHMMDB_UcfVirus):
				raise ValueError("Can't find HMMER PPHMM database of unclassified viruses")

		ClassificationResultFile = self.VariableShelveDir_UcfVirus+"/ClassificationResults.txt"
	
		return ClassificationResultFile

	def use_ucf_virus_pphmms(self):
		'''4/ : Scan unclassified viruses against their PPHMM DB to create additional PPHMMSignatureTable, and PPHMMLocationTable'''
		'''Generate PPHMMSignatureTable, and PPHMMLocationTable'''
		'''Make HMMER_hmmscanDir - do not return as regenerated in later functions'''
		HMMER_hmmscanDir_UcfVirus = f"{self.HMMERDir_UcfVirus}/hmmscan_"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMER_hmmscanDir_UcfVirus)
		
		'''Generate PPHMMSignatureTable and PPHMMLocationTable'''
		PPHMMSignatureTable_UcfVirusVSUcfDB, \
		PPHMMLocationTable_UcfVirusVSUcfDB = PPHMMSignatureTable_Constructor(SeqIDLists = self.ucf_genomes["SeqIDLists"],
													GenBankFile			= self.GenomeSeqFile_UcfVirus,
													TranslTableList		= self.ucf_genomes["TranslTableList"],
													SeqLength_Cutoff	= self.SeqLength_Cutoff,
													HMMER_PPHMMDB		= self.HMMER_PPHMMDB_UcfVirus,
													HMMER_hmmscanDir	= HMMER_hmmscanDir_UcfVirus,
													HMMER_N_CPUs		= self.HMMER_N_CPUs,
													HMMER_C_EValue_Cutoff	= self.HMMER_C_EValue_Cutoff,
													HMMER_HitScore_Cutoff	= self.HMMER_HitScore_Cutoff)
		
		'''Delete HMMER_hmmscanDir'''
		_ = subprocess.call("rm -rf %s" %HMMER_hmmscanDir_UcfVirus, shell = True)
		
		'''Update unclassified viruses' PPHMMSignatureTables, and PPHMMLocationTables'''
		self.ucf_annots["PPHMMSignatureTable_Dict"] = {RefVirusGroup: np.hstack((PPHMMSignatureTable_UcfVirusVSRefDB, PPHMMSignatureTable_UcfVirusVSUcfDB)) for RefVirusGroup, PPHMMSignatureTable_UcfVirusVSRefDB in self.ucf_annots["PPHMMSignatureTable_Dict"].items()}
		self.ucf_annots["PPHMMLocationTable_Dict"]  = {RefVirusGroup: np.hstack((PPHMMLocationTable_UcfVirusVSRefDB, PPHMMLocationTable_UcfVirusVSUcfDB)) for RefVirusGroup, PPHMMLocationTable_UcfVirusVSRefDB in self.ucf_annots["PPHMMLocationTable_Dict"].items()}

	def update_globals(self):
		'''5/ : Update global vars with loaded ucf genome data'''
		self.TaxoLabelList_UcfVirus	= list(map('_'.join, list(zip(	["Query%.4d"%i for i in range(len(self.ucf_genomes["SeqIDLists"]))],
						["/".join(SeqIDList) if len(SeqIDList)<=3 else "/".join(SeqIDList[0:3])+"/..." for SeqIDList in self.ucf_genomes["SeqIDLists"]]
						))))
		
		self.N_UcfViruses, self.N_RefVirusGroups = len(self.ucf_genomes["SeqIDLists"]), len(self.ShelveDirs_RefVirus)
		self.MaxSimScoreTable, self.TaxoOfMaxSimScoreTable, self.TaxoAssignmentTable, \
			self.PhyloStatTable = np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0))
		if self.Bootstrap == True:
			self.PairwiseSimilarityScore_CutoffDist_Dict = {Bootstrap_i:{} for Bootstrap_i in range(self.N_Bootstrap)}
			self.MaxSimScoreDistTable = self.TaxoOfMaxSimScoreDistTable	= self.TaxoAssignmentDistTable = self.PhyloStatDistTable \
				 = np.zeros((self.N_Bootstrap, self.N_UcfViruses, 0))

	def classify(self):
		'''6/ : RM << DOCSTRING'''
		RefVirusGroup_i = 0
		for ShelveDir_RefVirus, GenomeSeqFile_RefVirus in zip(self.ShelveDirs_RefVirus, self.GenomeSeqFiles_RefVirus):
			# Parameters / from read_desc_tables
			#"SeqIDLists", "FamilyList", "GenusList", "VirusNameList", "TaxoGroupingList", "TranslTableList",
			# Parameters / from ref virus 
			# "PPHMMSignatureTable", "PPHMMLocationTable", "PPHMMSignatureTable_coo", "PPHMMLocationTable_coo", "GOMSignatureTable", "GOMIDList"

			RefVirusGroup_i += 1
			RefVirusGroup = ShelveDir_RefVirus.split("/")[-1]
			print(f"Classify viruses using {RefVirusGroup} as reference ({RefVirusGroup_i}/{self.N_RefVirusGroups})")
		
			'''Load in genomes variables'''
			VariableShelveDir_RefVirus = ShelveDir_RefVirus+"/Shelves"
			Parameters 				   = retrieve_genome_vars(VariableShelveDir_RefVirus, self.IncludeIncompleteGenomes_RefVirus)		
			ref_virus_genomes 		   = retrieve_ref_virus_vars(VariableShelveDir_RefVirus, self.IncludeIncompleteGenomes_RefVirus)
			Parameters 				   = Parameters | ref_virus_genomes
			N_RefViruses			   = len(Parameters["SeqIDLists"])

			if "PPHMMSignatureTable_coo" in Parameters.keys(): 	
				Parameters["PPHMMSignatureTable"] = Parameters["PPHMMSignatureTable_coo"].toarray()
			if "PPHMMLocationTable_coo" in Parameters.keys():	
				Parameters["PPHMMLocationTable"] = Parameters["PPHMMLocationTable_coo"].toarray()

			if self.UseUcfVirusPPHMMs == True:
				'''Scan reference viruses against the PPHMM database of unclassified viruses to generate additional PPHMMSignatureTable, and PPHMMLocationTable'''
				'''Generate PPHMMSignatureTable, and PPHMMLocationTable'''
				'''Make OR update HMMER_hmmscanDir'''
				HMMER_hmmscanDir_UcfVirus = f"{self.HMMERDir_UcfVirus}/hmmscan_"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMER_hmmscanDir_UcfVirus)
				
				'''Generate PPHMMSignatureTable and PPHMMLocationTable'''
				PPHMMSignatureTable_UcfVirusVSUcfDB, \
				PPHMMLocationTable_UcfVirusVSUcfDB = PPHMMSignatureTable_Constructor(SeqIDLists	= Parameters["SeqIDLists"],
																GenBankFile			= GenomeSeqFile_RefVirus,
																TranslTableList		= Parameters["TranslTableList"],
																SeqLength_Cutoff	= self.SeqLength_Cutoff,
																HMMER_PPHMMDB		= self.HMMER_PPHMMDB_UcfVirus,
																HMMER_hmmscanDir	= HMMER_hmmscanDir_UcfVirus,
																HMMER_N_CPUs		= self.HMMER_N_CPUs,
																HMMER_C_EValue_Cutoff	= self.HMMER_C_EValue_Cutoff,
																HMMER_HitScore_Cutoff	= self.HMMER_HitScore_Cutoff)
				
				'''Delete HMMER_hmmscanDir'''
				_ = subprocess.call("rm -rf %s" %HMMER_hmmscanDir_UcfVirus, shell = True)

				Parameters["PPHMMSignatureTable"]	= np.hstack((Parameters["PPHMMSignatureTable"], PPHMMSignatureTable_UcfVirusVSUcfDB))
				Parameters["PPHMMLocationTable"]	= np.hstack((Parameters["PPHMMLocationTable"],  PPHMMLocationTable_UcfVirusVSUcfDB))				
				UpdatedGOMDB_RefVirus				= GOMDB_Constructor(Parameters["TaxoGroupingList"], Parameters["PPHMMLocationTable"], Parameters["GOMIDList"])
								
				Parameters["GOMSignatureTable"]	= GOMSignatureTable_Constructor(Parameters["PPHMMLocationTable"], UpdatedGOMDB_RefVirus, Parameters["GOMIDList"])
				
				'''Update unclassified viruses' GOMSignatureTable'''
				self.ucf_annots["GOMSignatureTable_Dict"][RefVirusGroup] = GOMSignatureTable_Constructor(self.ucf_annots["PPHMMLocationTable_Dict"][RefVirusGroup], UpdatedGOMDB_RefVirus, Parameters["GOMIDList"])
			
			'''Build the dendrogram, including all sequences'''
			'''Generate TaxoLabelList of reference viruses'''
			TaxoLabelList_RefVirus	= TaxoLabel_Constructor(Parameters["SeqIDLists"], Parameters["FamilyList"], Parameters["GenusList"], Parameters["VirusNameList"])
			
			'''Compute pairwise distances'''
			PPHMMSignatureTable_AllVirus= np.vstack((Parameters["PPHMMSignatureTable"], self.ucf_annots["PPHMMSignatureTable_Dict"][RefVirusGroup]))
			GOMSignatureTable_AllVirus	= np.vstack((Parameters["GOMSignatureTable"], self.ucf_annots["GOMSignatureTable_Dict"][RefVirusGroup]))
			PPHMMLocationTable_AllVirus	= np.vstack((Parameters["PPHMMLocationTable"], self.ucf_annots["PPHMMLocationTable_Dict"][RefVirusGroup]))
			TaxoLabelList_AllVirus		= TaxoLabelList_RefVirus + self.TaxoLabelList_UcfVirus
			SimMat = SimilarityMat_Constructor(PPHMMSignatureTable_AllVirus, GOMSignatureTable_AllVirus, PPHMMLocationTable_AllVirus, self.SimilarityMeasurementScheme, self.p)
			DistMat				= 1 - SimMat
			DistMat[DistMat<0]	= 0
			
			'''Generate dendrogram'''
			VirusDendrogram		= DistMat2Tree(DistMat, TaxoLabelList_AllVirus, self.Dendrogram_LinkageMethod)
			VirusDendrogramFile = f"{self.VariableShelveDir_UcfVirus}/Dendrogram.RefVirusGroup={RefVirusGroup}.IncompleteUcfRefGenomes={str(int(self.IncludeIncompleteGenomes_UcfVirus))+str(int(self.IncludeIncompleteGenomes_RefVirus))}.Scheme={self.SimilarityMeasurementScheme}.Method={self.Dendrogram_LinkageMethod}.p={self.p}.nwk"
			
			with open(VirusDendrogramFile, "w") as VirusDendrogram_txt:
				VirusDendrogram_txt.write(VirusDendrogram)
			
			'''Compute similarity cut off for each taxonomic class'''
			self.PairwiseSimilarityScore_Cutoff_Dict[RefVirusGroup] = PairwiseSimilarityScore_Cutoff_Dict_Constructor(SimMat[:N_RefViruses][:,:N_RefViruses], Parameters["TaxoGroupingList"], self.N_PairwiseSimilarityScores)
			
			'''Propose a taxonomic class to each unclassified virus and evaluate the taxonomic assignments'''
			MaxSimScoreList, TaxoOfMaxSimScoreList, \
			TaxoAssignmentList, PhyloStatList = TaxonomicAssignmentProposerAndEvaluator(SimMat[-self.N_UcfViruses:][:,:N_RefViruses],
																		Parameters["TaxoGroupingList"],
																		VirusDendrogram,
																		TaxoLabelList_RefVirus,
																		self.TaxoLabelList_UcfVirus,
																		self.PairwiseSimilarityScore_Cutoff_Dict[RefVirusGroup])
			
			MaxSimScoreTable		= np.column_stack((MaxSimScoreTable, MaxSimScoreList))
			TaxoOfMaxSimScoreTable	= np.column_stack((TaxoOfMaxSimScoreTable, TaxoOfMaxSimScoreList))
			TaxoAssignmentTable		= np.column_stack((TaxoAssignmentTable, TaxoAssignmentList))
			PhyloStatTable			= np.column_stack((PhyloStatTable, PhyloStatList))
			
			if self.VirusGrouping == True:
				'''(OPT) Virus grouping'''
				self.do_virus_groupings(RefVirusGroup, DistMat, N_RefViruses, Parameters)

			if self.Bootstrap == True:
				'''(OPT) Bootstrap stats'''
				self.construct_bootstrap_distributions(Parameters, RefVirusGroup, PPHMMSignatureTable_AllVirus, PPHMMLocationTable_AllVirus, N_RefViruses, TaxoLabelList_AllVirus, TaxoLabelList_RefVirus, VirusDendrogramFile)

			if self.Heatmap_WithDendrogram == True:
				'''(OPT) Draw hm/dendro'''
				self.heatmap_with_dendrogram(self.Bootstrap)

	def do_virus_groupings(self, RefVirusGroup, DistMat, N_RefViruses, Parameters) -> None:
		'''Accessory fn to CLASSIFY (6):  Group viruses with scipy fcluster/linkage, save results to txt'''
		self.VirusGroupingDict[RefVirusGroup] = {}
		VirusGroupingFile = f"{self.VariableShelveDir_UcfVirus}/VirusGrouping.RefVirusGroup={RefVirusGroup}.IncompleteUcfRefGenomes={str(int(self.IncludeIncompleteGenomes_UcfVirus))+str(int(self.IncludeIncompleteGenomes_RefVirus))}.Scheme={self.SimilarityMeasurementScheme}.Method={self.Dendrogram_LinkageMethod}.p={self.p}.txt"
		
		'''Compute distance cutoff'''
		_, OptDistance_Cutoff, CorrelationScore, Theils_u_TaxoGroupingListGivenPred, \
			Theils_u_PredGivenTaxoGroupingList = VirusGrouping_Estimator(DistMat[:N_RefViruses][:,:N_RefViruses], self.Dendrogram_LinkageMethod, Parameters["TaxoGroupingList"])
		
		'''Virus grouping'''	
		VirusGroupingList = fcluster(linkage(DistMat[np.triu_indices_from(DistMat, k = 1)], self.Dendrogram_LinkageMethod), OptDistance_Cutoff, 'distance')
		self.VirusGroupingDict[RefVirusGroup]["VirusGroupingList"] = VirusGroupingList
		self.VirusGroupingDict[RefVirusGroup]["OptDistance_Cutoff"] = OptDistance_Cutoff
		self.VirusGroupingDict[RefVirusGroup]["CorrelationScore"] = CorrelationScore
		self.VirusGroupingDict[RefVirusGroup]["Theils_u_TaxoGroupingListGivenPred"] = Theils_u_TaxoGroupingListGivenPred
		self.VirusGroupingDict[RefVirusGroup]["Theils_u_PredGivenTaxoGroupingList"] = Theils_u_PredGivenTaxoGroupingList
		
		'''Save result to txt file'''
		np.savetxt(fname = VirusGroupingFile,
					X = np.column_stack((list(map(", ".join, Parameters["SeqIDLists"].tolist()+self.ucf_genomes["SeqIDLists"].tolist())),
								Parameters["FamilyList"].tolist()+[""]*self.N_UcfViruses,
								Parameters["GenusList"].tolist()+[""]*self.N_UcfViruses,
								Parameters["VirusNameList"].tolist()+self.ucf_genomes["VirusNameList"].tolist(),
								Parameters["TaxoGroupingList"].tolist()+[""]*self.N_UcfViruses,
								VirusGroupingList,
								)),
					fmt = '%s',
					delimiter = "\t",
					header = "Sequence identifier\tFamily\tGenus\tVirus name\tClass\tGrouping")
		
		with open(VirusGroupingFile, "a") as VirusGrouping_txt:
			VirusGrouping_txt.write("\n"+
						f"Distance cut off: {OptDistance_Cutoff}\n"
						f"Theil's uncertainty correlation for the reference assignments given the predicted grouping U(Ref|Pred): {Theils_u_TaxoGroupingListGivenPred}\n"
						f"Theil's uncertainty correlation for the predicted grouping given the reference assignments U(Pred|Ref): {Theils_u_PredGivenTaxoGroupingList}\n"
						f"Symmetrical Theil's uncertainty correlation between the reference assignments and the predicted grouping U(Ref, Pred): {CorrelationScore}\n"
						f"U(X|Y) == 1 means that knowing Y implies a perfect knowledge of X, but not vice-versa\n"
						f"U(X,Y) == 1 means that knowing Y implies a perfect knowledge of X and vice-versa\n"
						)

	def construct_bootstrap_distributions(self, Parameters, RefVirusGroup, PPHMMSignatureTable_AllVirus, PPHMMLocationTable_AllVirus, N_RefViruses, TaxoLabelList_AllVirus, TaxoLabelList_RefVirus, VirusDendrogramFile):
		'''Accessory fn to CLASSIFY (6). RM << DOCSTRING'''
		'''Construct result distributions by using the bootstrapping technique'''
		'''Define path to dendrogram distribution file'''
		VirusDendrogramDistFile	= f"{self.VariableShelveDir_UcfVirus}/DendrogramDist.RefVirusGroup={RefVirusGroup}.IncompleteUcfRefGenomes={str(int(self.IncludeIncompleteGenomes_UcfVirus))+str(int(self.IncludeIncompleteGenomes_RefVirus))}.Scheme={self.SimilarityMeasurementScheme}.Method={self.Dendrogram_LinkageMethod}.p={self.p}.nwk"
		if os.path.isfile(VirusDendrogramDistFile):
			os.remove(VirusDendrogramDistFile)
		
		'''Define path to bootstrapped dendrogram file'''
		BootstrappedVirusDendrogramFile	= f"{self.VariableShelveDir_UcfVirus}/BootstrappedDendrogram.RefVirusGroup={RefVirusGroup}.IncompleteUcfRefGenomes={str(int(self.IncludeIncompleteGenomes_UcfVirus))+str(int(self.IncludeIncompleteGenomes_RefVirus))}.Scheme={self.SimilarityMeasurementScheme}.Method={self.Dendrogram_LinkageMethod}.p={self.p}.nwk"
		if os.path.isfile(BootstrappedVirusDendrogramFile):
			os.remove(BootstrappedVirusDendrogramFile)
		
		BootsrappedTaxoAssignmentTable, BootsrappedMaxSimScoreTable, BootsrappedTaxoOfMaxSimScoreTable, \
			BootsrappedPhyloStatTable = np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0)), np.zeros((self.N_UcfViruses, 0))
		N_PPHMMs					   = PPHMMSignatureTable_AllVirus.shape[1]
		for Bootstrap_i in range(0, self.N_Bootstrap):
			'''Bootstrap the data'''
			print(f"Round {Bootstrap_i+1}")

			'''Construct bootstrapped PPHMMSignatureTable and PPHMMLocationTable'''
			PPHMM_IndexList = sorted(np.random.choice(list(range(N_PPHMMs)), N_PPHMMs, replace = True))
			BootstrappedPPHMMSignatureTable = PPHMMSignatureTable_AllVirus[:,PPHMM_IndexList]
			BootstrappedPPHMMLocationTable = PPHMMLocationTable_AllVirus[:,PPHMM_IndexList]
			BootstrappedGOMSignatureTable = None
			if "G" in self.SimilarityMeasurementScheme:
				'''Construct bootstrapped GOMSignatureTable'''
				BootstrappedGOMDB = GOMDB_Constructor(Parameters["TaxoGroupingList"], BootstrappedPPHMMLocationTable[:N_RefViruses], Parameters["GOMIDList"])
				BootstrappedGOMSignatureTable = GOMSignatureTable_Constructor(BootstrappedPPHMMLocationTable, BootstrappedGOMDB,Parameters["GOMIDList"])
			
			'''Construct a dendrogram from the bootstrapped data'''
			BootstrappedSimMat  = SimilarityMat_Constructor(BootstrappedPPHMMSignatureTable, BootstrappedGOMSignatureTable, BootstrappedPPHMMLocationTable, self.SimilarityMeasurementScheme, self.p)
			BootstrappedDistMat = 1 - BootstrappedSimMat
			BootstrappedDistMat[BootstrappedDistMat<0] = 0
			BootstrappedVirusDendrogram = DistMat2Tree(BootstrappedDistMat, TaxoLabelList_AllVirus, self.Dendrogram_LinkageMethod)
			
			with open(VirusDendrogramDistFile, "a") as VirusDendrogramDist_txt:
				VirusDendrogramDist_txt.write(BootstrappedVirusDendrogram+"\n")
			
			'''Compute similarity cut off for each taxonomic class based on the bootstrapped data'''
			self.PairwiseSimilarityScore_CutoffDist_Dict[Bootstrap_i][RefVirusGroup] = PairwiseSimilarityScore_Cutoff_Dict_Constructor (BootstrappedSimMat[:N_RefViruses][:,:N_RefViruses],
																						Parameters["TaxoGroupingList"],
																						self.N_PairwiseSimilarityScores)
			
			'''Propose a taxonomic class to each unclassified virus and evaluate the taxonomic assignments based on the bootstrapped data'''
			BootsrappedMaxSimScoreList, BootsrappedTaxoOfMaxSimScoreList, BootsrappedTaxoAssignmentList, \
				BootsrappedPhyloStatList = TaxonomicAssignmentProposerAndEvaluator(BootstrappedSimMat[-self.N_UcfViruses:][:,:N_RefViruses],
												Parameters["TaxoGroupingList"],
												BootstrappedVirusDendrogram,
												TaxoLabelList_RefVirus,
												self.TaxoLabelList_UcfVirus,
												self.PairwiseSimilarityScore_CutoffDist_Dict[Bootstrap_i][RefVirusGroup])
			
			BootsrappedMaxSimScoreTable			= np.column_stack((BootsrappedMaxSimScoreTable, BootsrappedMaxSimScoreList))
			BootsrappedTaxoOfMaxSimScoreTable	= np.column_stack((BootsrappedTaxoOfMaxSimScoreTable, BootsrappedTaxoOfMaxSimScoreList))
			BootsrappedTaxoAssignmentTable		= np.column_stack((BootsrappedTaxoAssignmentTable, BootsrappedTaxoAssignmentList))
			BootsrappedPhyloStatTable			= np.column_stack((BootsrappedPhyloStatTable, BootsrappedPhyloStatList))
		
		MaxSimScoreDistTable		= np.append(MaxSimScoreDistTable, BootsrappedMaxSimScoreTable.T[...,None],	axis = 2)
		TaxoOfMaxSimScoreDistTable	= np.append(TaxoOfMaxSimScoreDistTable,	BootsrappedTaxoOfMaxSimScoreTable.T[...,None],	axis = 2)
		TaxoAssignmentDistTable		= np.append(TaxoAssignmentDistTable, 	BootsrappedTaxoAssignmentTable.T[...,None],	axis = 2)
		PhyloStatDistTable			= np.append(PhyloStatDistTable,	  BootsrappedPhyloStatTable.T[...,None],		axis = 2)
		
		'''Create a bootstrapped dendrogram'''
		if self.Bootstrap_method == "booster":
			_ = subprocess.Popen(f"booster -i {VirusDendrogramFile} -b {VirusDendrogramDistFile} -o {BootstrappedVirusDendrogramFile} -@ {self.Bootstrap_N_CPUs} ",
										stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()

		elif self.Bootstrap_method == "sumtrees":
			_ = subprocess.Popen(f"sumtrees.py --decimals=2 --no-annotations --preserve-underscores --force-rooted --output-tree-format=newick --output-tree-filepath={BootstrappedVirusDendrogramFile} --target={VirusDendrogramFile} {VirusDendrogramDistFile}",
																									stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()

		else:
			print("WARNING ** Bootstrap dendrogram not constructed: 'Bootstrap_method' can either be 'booster' or 'sumtrees'.")

	def heatmap_with_dendrogram(self, Parameters, VirusDendrogramFile, TaxoLabelList_AllVirus, RefVirusGroup, DistMat, N_RefViruses, TaxoLabelList_RefVirus, TaxoAssignmentList):
		'''Construct GRAViTy heat map with dendrogram'''
		'''Load the tree'''
		if self.Bootstrap == True:
			BootstrappedVirusDendrogramFile	= f"{self.VariableShelveDir_UcfVirus}/BootstrappedDendrogram.RefVirusGroup={RefVirusGroup}.IncompleteUcfRefGenomes={str(int(self.IncludeIncompleteGenomes_UcfVirus))+str(int(self.IncludeIncompleteGenomes_RefVirus))}.Scheme={self.SimilarityMeasurementScheme}.Method={self.Dendrogram_LinkageMethod}.p={self.p}.nwk"
			VirusDendrogram		= Phylo.read(BootstrappedVirusDendrogramFile, "newick")
		else:
			VirusDendrogram		= Phylo.read(VirusDendrogramFile, "newick")
		
		'''Determine virus order'''
		_ 					= VirusDendrogram.ladderize(reverse = True)
		OrderedTaxoLabelList= [Clade.name for Clade in VirusDendrogram.get_terminals()]
		VirusOrder			= [TaxoLabelList_AllVirus.index(TaxoLabel) for TaxoLabel in OrderedTaxoLabelList]
		
		'''Re-order the distance matrix'''
		OrderedDistMat	= DistMat[VirusOrder][:,VirusOrder]
		
		'''Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff'''
		N_InternalNodes = len(VirusDendrogram.get_nonterminals())
		for InternalNode_i in range(N_InternalNodes):
			if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence >= self.Heatmap_DendrogramSupport_Cutoff:
				VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
			else:
				VirusDendrogram.get_nonterminals()[InternalNode_i].confidence=""
		
		'''Colour terminal branches: reference virus's branch is blue, unclassified virus's branch is red'''
		N_Viruses = N_RefViruses + self.N_UcfViruses
		for Virus_i in range(N_Viruses):
			Taxolabel = VirusDendrogram.get_terminals()[Virus_i].name
			if Taxolabel in TaxoLabelList_RefVirus:
				VirusDendrogram.get_terminals()[Virus_i].color="blue"
			elif Taxolabel in self.TaxoLabelList_UcfVirus:
				VirusDendrogram.get_terminals()[Virus_i].color="red"
		
		'''Labels, label positions, and ticks'''
		TaxoGroupingList_AllVirus	= Parameters["TaxoGroupingList"].tolist() + TaxoAssignmentList
		Taxo2ClassDict		= {TaxoLabel: TaxoGrouping for TaxoLabel, TaxoGrouping in zip(TaxoLabelList_AllVirus, TaxoGroupingList_AllVirus)}
		ClassDendrogram		= copy(VirusDendrogram)
		for Clade in ClassDendrogram.find_clades(terminal=True):
			Clade.name = Taxo2ClassDict[Clade.name]
		
		ClassLabelList, LineList = [], [-1]
		TerminalNodeList = [TerminalNode for TerminalNode in ClassDendrogram.get_terminals()]
		while len(TerminalNodeList)!=0:
			FarLeftNode = TerminalNodeList[0]
			for Clade in ([ClassDendrogram]+ClassDendrogram.get_path(FarLeftNode)):
				DescendantNodeList = Clade.get_terminals()
				DescendantClassLabelList = list(set([c.name for c in DescendantNodeList]))
				if len(DescendantClassLabelList)==1:
					ClassLabelList.append(DescendantClassLabelList[0])
					LineList.append(LineList[-1]+len(DescendantNodeList))
					TerminalNodeList = TerminalNodeList[len(DescendantNodeList):]
					break
		
		ClassLabelList	= np.array(ClassLabelList)
		LineList		= np.array(LineList) + 0.5
		TickLocList		= np.array(list(map(np.mean, list(zip(LineList[0:-1],LineList[1:])))))
		
		'''Heat map colour indicators'''
		IndicatorMat_RefVirus	= np.tile([True]*N_RefViruses + [False]*self.N_UcfViruses, N_Viruses).reshape(N_Viruses,-1)
		IndicatorMat_RefVirus	= IndicatorMat_RefVirus * IndicatorMat_RefVirus.T
		IndicatorMat_RefVirus	= IndicatorMat_RefVirus[VirusOrder][:,VirusOrder]
		IndicatorMat_UcfVirus	= np.tile([False]*N_RefViruses + [True]*self.N_UcfViruses, N_Viruses).reshape(N_Viruses,-1)
		IndicatorMat_UcfVirus	= IndicatorMat_UcfVirus * IndicatorMat_UcfVirus.T
		IndicatorMat_UcfVirus	= IndicatorMat_UcfVirus[VirusOrder][:,VirusOrder]
		IndicatorMat_CrossGroup	= ~IndicatorMat_RefVirus * ~IndicatorMat_UcfVirus
		
		'''Masked OrderedDistMat'''
		OrderedDistMat_RefVirus = np.ma.masked_where(~IndicatorMat_RefVirus, OrderedDistMat)
		OrderedDistMat_UcfVirus = np.ma.masked_where(~IndicatorMat_UcfVirus, OrderedDistMat)
		OrderedDistMat_CrossGroup = np.ma.masked_where(~IndicatorMat_CrossGroup, OrderedDistMat)
		
		'''Colour map construction'''
		BlueColour_Dict = {
			'red'  :  ((0., 0., 0.), (1., 1., 1.)),
			'green':  ((0., 0., 0.), (1., 1., 1.)),
			'blue' :  ((0., 1., 1.), (1., 1., 1.))
			}
		
		RedColour_Dict = {
			'red'  :  ((0., 1., 1.), (1., 1., 1.)),
			'green':  ((0., 0., 0.), (1., 1., 1.)),
			'blue' :  ((0., 0., 0.), (1., 1., 1.))
			}
		
		PurpleColour_Dict = {
			'red'  :  ((0., 0.5, 0.5), (1., 1., 1.)),
			'green':  ((0., 0., 0.), (1., 1., 1.)),
			'blue' :  ((0., 1., 1.), (1., 1., 1.))
			}
		
		MyBlues		= LSC('MyBlues', BlueColour_Dict, 1024)
		MyReds 		= LSC('MyReds', RedColour_Dict, 1024)
		MyPurples 	= LSC('MyPurples', PurpleColour_Dict, 1024)
		
		'''Plot configuration'''
		Heatmap_width		= float(12)
		Heatmap_height		= Heatmap_width
		TaxoLable_space		= 1.00

		CBar_Heatmap_gap	= 0.05
		CBar_width			= Heatmap_width
		CBar_height			= 0.50
		CBarLable_space		= 0.25
		
		Dendrogram_width	= Heatmap_width/3
		Dendrogram_height	= Heatmap_height
		Dendrogram_Heatmap_gap	= 0.1
		
		ScaleBar_Dendrogram_gap	= CBar_Heatmap_gap
		ScaleBar_width		= Dendrogram_width
		ScaleBar_height		= CBar_height
		ScaleBarLable_space	= CBarLable_space
		
		Outer_margin		= 0.5
		FontSize			= 6
		
		Fig_width			= Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap + Heatmap_width + TaxoLable_space + Outer_margin
		Fig_height			= Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap + Heatmap_height + TaxoLable_space + Outer_margin
		
		ax_Dendrogram_L		= Outer_margin/Fig_width
		ax_Dendrogram_B		= (Outer_margin + ScaleBarLable_space + ScaleBar_height + ScaleBar_Dendrogram_gap)/Fig_height
		ax_Dendrogram_W		= Dendrogram_width/Fig_width
		ax_Dendrogram_H		= Dendrogram_height/Fig_height
		
		ax_ScaleBar_L		= Outer_margin/Fig_width
		ax_ScaleBar_B		= (Outer_margin + ScaleBarLable_space)/Fig_height
		ax_ScaleBar_W		= ScaleBar_width/Fig_width
		ax_ScaleBar_H		= ScaleBar_height/Fig_height
		
		ax_Heatmap_L		= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
		ax_Heatmap_B		= (Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap)/Fig_height
		ax_Heatmap_W		= Heatmap_width/Fig_width
		ax_Heatmap_H		= Heatmap_height/Fig_height
		
		ax_CBar_L			= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
		ax_CBar_B			= (Outer_margin + CBarLable_space)/Fig_height
		ax_CBar_W			= CBar_width/Fig_width
		ax_CBar_H			= CBar_height/Fig_height
		
		'''Plot the heat map'''
		fig					= plt.figure(figsize = (Fig_width, Fig_height), dpi = 300)
		
		ax_Dendrogram		= fig.add_axes([ax_Dendrogram_L, ax_Dendrogram_B, ax_Dendrogram_W, ax_Dendrogram_H], frame_on = False, facecolor = "white")
		Phylo				.draw(VirusDendrogram, label_func = lambda x: "", do_show = False,  axes = ax_Dendrogram)
		VirusDendrogramDepth= max([v for k,v in VirusDendrogram.depths().items()])
		ax_Dendrogram		.set_xlim([(VirusDendrogramDepth - 1), VirusDendrogramDepth])
		ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
		ax_Dendrogram		.set_axis_off()
		
		ax_ScaleBar		= fig.add_axes([ax_ScaleBar_L, ax_ScaleBar_B, ax_ScaleBar_W, ax_ScaleBar_H], frame_on = False, facecolor = "white")
		ax_ScaleBar		.plot([0,1],[0,0],'k-')
		ScaleBarTicks	= [0, 0.25, 0.5, 0.75, 1]
		for Tick in ScaleBarTicks:
			ax_ScaleBar.plot([Tick, Tick],[-0.05, 0.05],'k-')
		
		ax_ScaleBar		.set_xlim([1, 0])
		ax_ScaleBar		.set_xticks(ScaleBarTicks)
		ax_ScaleBar		.set_xticklabels(list(map(str, ScaleBarTicks)), rotation = 0, size = FontSize)
		ax_ScaleBar		.set_xlabel('Distance', rotation = 0, size = FontSize+2)
		ax_ScaleBar		.xaxis.set_label_position('bottom')
		ax_ScaleBar		.tick_params(	top = 'off',
							bottom = 'off',
							left = 'off',
							right = 'off',
							labeltop = 'off',
							labelbottom = 'on',
							labelleft = 'off',
							labelright = 'off',
							direction = 'out')
		
		ax_Heatmap		= fig.add_axes([ax_Heatmap_L, ax_Heatmap_B, ax_Heatmap_W, ax_Heatmap_H], frame_on = True, facecolor = "white")
		ax_Heatmap		.imshow(OrderedDistMat_RefVirus, cmap = MyBlues, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
		ax_Heatmap		.imshow(OrderedDistMat_UcfVirus, cmap = MyReds, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
		ax_Heatmap		.imshow(OrderedDistMat_CrossGroup, cmap = MyPurples, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')

		for l in LineList:
			ax_Heatmap.axvline(l, color = 'k', lw = 0.2)
			ax_Heatmap.axhline(l, color = 'k', lw = 0.2)
		
		ax_Heatmap		.set_xticks(TickLocList)
		ax_Heatmap		.set_xticklabels(ClassLabelList, rotation = 90, size = FontSize)
		ax_Heatmap		.set_yticks(TickLocList)
		ax_Heatmap		.set_yticklabels(ClassLabelList, rotation = 0, size = FontSize)
		ax_Heatmap		.tick_params(	top = 'on',
							bottom = 'off',
							left = 'off',
							right = 'on',
							labeltop = 'on',
							labelbottom = 'off',
							labelleft = 'off',
							labelright = 'on',
							direction = 'out')
		
		ax_CBar_RefVirus	= fig.add_axes([ax_CBar_L, ax_CBar_B + 2*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
		ax_CBar_RefVirus	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyBlues, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
		ax_CBar_RefVirus	.set_yticks([0.0])
		ax_CBar_RefVirus	.set_yticklabels(["Ref viruses"], rotation = 0, size = FontSize + 2)
		ax_CBar_RefVirus	.tick_params(	top = 'off',
							bottom = 'off',
							left = 'off',
							right = 'on',
							labeltop = 'off',
							labelbottom = 'off',
							labelleft = 'off',
							labelright = 'on',
							direction = 'out')

		ax_CBar_UcfVirus	= fig.add_axes([ax_CBar_L, ax_CBar_B + 1*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
		ax_CBar_UcfVirus	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyReds, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
		ax_CBar_UcfVirus	.set_yticks([0.0])
		ax_CBar_UcfVirus	.set_yticklabels(["Ucf viruses"], rotation = 0, size = FontSize + 2)
		ax_CBar_UcfVirus	.tick_params(	top = 'off',
							bottom = 'off',
							left = 'off',
							right = 'on',
							labeltop = 'off',
							labelbottom = 'off',
							labelleft = 'off',
							labelright = 'on',
							direction = 'out')

		ax_CBar_CrossGroup	= fig.add_axes([ax_CBar_L, ax_CBar_B + 0*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
		ax_CBar_CrossGroup	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyPurples, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
		ax_CBar_CrossGroup	.set_yticks([0.0])
		ax_CBar_CrossGroup	.set_yticklabels(["Ref VS Ucf viruses"], rotation = 0, size = FontSize + 2)
		ax_CBar_CrossGroup	.set_xticks(np.array([0, 0.25, 0.50, 0.75, 1])*(1025)-0.5)
		ax_CBar_CrossGroup	.set_xticklabels(['0', '0.25', '0.50', '0.75', '1'], rotation = 0, size = FontSize)
		ax_CBar_CrossGroup	.set_xlabel("Distance", rotation = 0, size = FontSize + 2)
		ax_CBar_CrossGroup	.tick_params(	top = 'off',
							bottom = 'on',
							left = 'off',
							right = 'on',
							labeltop = 'off',
							labelbottom = 'on',
							labelleft = 'off',
							labelright = 'on',
							direction = 'out')
		
		'''Save the plot to file'''
		HeatmapWithDendrogramFile = f"{self.VariableShelveDir_UcfVirus}/HeatmapWithDendrogram.RefVirusGroup={RefVirusGroup}.IncompleteUcfRefGenomes={str(int(self.IncludeIncompleteGenomes_UcfVirus))+str(int(self.IncludeIncompleteGenomes_RefVirus))}.Scheme={self.SimilarityMeasurementScheme}.Method={self.Dendrogram_LinkageMethod}.p={self.p}.pdf"
		plt.savefig(HeatmapWithDendrogramFile, format = "pdf")

	def main(self):
		'''Classify viruses and evaluate results'''
		section_header("#Classify virus and evaluate the results")

		if self.UseUcfVirusPPHMMs == True:
			'''1/ : (OPT) Guide user to specify input PPHMM file'''
			self.input_test()
		else:
			self.GenomeSeqFiles_RefVirus = [None]*len(self.ShelveDirs_RefVirus)
		
		'''2/ : Make fpaths'''
		ClassificationResultFile = self.mkdirs()
		
		'''3/ : Retrieve variables'''
		self.ucf_genomes = retrieve_genome_vars(self.VariableShelveDir_UcfVirus, self.IncludeIncompleteGenomes_UcfVirus)
		self.ucf_annots  = retrieve_ucf_annots(self.VariableShelveDir_UcfVirus, self.IncludeIncompleteGenomes_UcfVirus, self.IncludeIncompleteGenomes_RefVirus)
		
		if "PPHMMSignatureTable_Dict_coo" in self.ucf_annots.keys():
			self.ucf_annots["PPHMMSignatureTable_Dict"] = {RefVirusGroup:PPHMMSignatureTable_coo.toarray() for RefVirusGroup, PPHMMSignatureTable_coo in self.ucf_annots["PPHMMSignatureTable_Dict_coo"].items()}
		if "PPHMMLocationTable_Dict_coo" in self.ucf_annots.keys():	
			self.ucf_annots["PPHMMLocationTable_Dict"] = {RefVirusGroup:PPHMMLocationTable_coo.toarray() for RefVirusGroup, PPHMMLocationTable_coo in self.ucf_annots["PPHMMLocationTable_Dict_coo"].items()}
		
		if self.UseUcfVirusPPHMMs == True:
			'''(OPT) 4/ : Scan unclassified viruses against PPHMMDB, create additional sig and loc table'''
			self.use_ucf_virus_pphmms()

		'''5/ : Update global vars with loaded data'''
		self.update_globals()

		'''6/ : Classify viruses'''
		self.classify()

		################################################################################
		print("- Pool results from all classfiers, and finalise the taxonomic assignment (and virus grouping)")
		################################################################################
		if VirusGrouping == True:
			FinalisedTaxoAssignmentList	= []
			FinalisedVirusGroupingList	= []
			FinalisedVirusGroupingDict	= {}
			FinalisedVirusGroupingIndexDict	= {}
			for UcfVirus_i in range(N_UcfViruses):
				MaxSimScore_i		= np.argmax(MaxSimScoreTable[UcfVirus_i])
				FinalisedTaxoAssignment	= TaxoAssignmentTable[UcfVirus_i][MaxSimScore_i]
				RefVirusGroup		= ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]
				if RefVirusGroup not in FinalisedVirusGroupingDict: FinalisedVirusGroupingDict[RefVirusGroup] = {}
				if RefVirusGroup not in FinalisedVirusGroupingIndexDict: FinalisedVirusGroupingIndexDict[RefVirusGroup] = 1
				
				if FinalisedTaxoAssignment != "Unclassified":
					FinalisedTaxoAssignmentList.append("%s (%s)"%(FinalisedTaxoAssignment, RefVirusGroup))
					FinalisedVirusGroupingList.append("%s (%s)"%(FinalisedTaxoAssignment, RefVirusGroup))
				else:
					if MaxSimScoreTable[UcfVirus_i][MaxSimScore_i] <= DatabaseAssignmentSimilarityScore_Cutoff:
						FinalisedTaxoAssignmentList.append("Unclassified (NA)")
						FinalisedVirusGroupingList.append("Unclassified (NA)")
					else:
						FinalisedTaxoAssignmentList.append("Unclassified (%s)"%RefVirusGroup)
						VirusGrouping = VirusGroupingDict[RefVirusGroup]["VirusGroupingList"][-N_UcfViruses:][UcfVirus_i]
						if VirusGrouping not in FinalisedVirusGroupingDict[RefVirusGroup]:
							FinalisedVirusGroupingDict[RefVirusGroup][VirusGrouping] = "UTU%.4d"%FinalisedVirusGroupingIndexDict[RefVirusGroup]
							FinalisedVirusGroupingIndexDict[RefVirusGroup] += 1
						
						FinalisedVirusGroupingList.append("%s (%s)"%(FinalisedVirusGroupingDict[RefVirusGroup][VirusGrouping], RefVirusGroup))
		else:
			FinalisedTaxoAssignmentList	= []
			FinalisedVirusGroupingList	= [""]*N_UcfViruses
			for UcfVirus_i in range(N_UcfViruses):
				MaxSimScore_i		= np.argmax(MaxSimScoreTable[UcfVirus_i])
				FinalisedTaxoAssignment	= TaxoAssignmentTable[UcfVirus_i][MaxSimScore_i]
				RefVirusGroup		= ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]
				if FinalisedTaxoAssignment != "Unclassified":
					FinalisedTaxoAssignmentList.append("%s (%s)"%(FinalisedTaxoAssignment, RefVirusGroup))
				else:
					if MaxSimScoreTable[UcfVirus_i][MaxSimScore_i] <= DatabaseAssignmentSimilarityScore_Cutoff:
						FinalisedTaxoAssignmentList.append("Unclassified (NA)")
					else:
						FinalisedTaxoAssignmentList.append("Unclassified (%s)"%RefVirusGroup)
		
		RemainedUcfVirus_IndexList = np.where(list(map(all, MaxSimScoreTable <= DatabaseAssignmentSimilarityScore_Cutoff)))[0]
		SeqIDLists_RemainedUcfVirus = self.ucf_genomes["SeqIDLists"][RemainedUcfVirus_IndexList]
		VirusNameList_RemainedUcfVirus = self.ucf_genomes["VirusNameList"][RemainedUcfVirus_IndexList]
		
		TaxoOfMaxSimScoreRangeTable	= None
		MaxSimScoreRangeTable		= None
		PhyloStatRangeTable		= None
		TaxoAssignmentRangeTable	= None
		FinalisedTaxoAssignmentRangeList= None

		# RM << DIsabled for dev
		if Bootstrap == True:
			TaxoOfMaxSimScoreRangeTable	= np.zeros((N_UcfViruses, 0))
			MaxSimScoreRangeTable		= np.zeros((N_UcfViruses, 0))
			PhyloStatRangeTable		= np.zeros((N_UcfViruses, 0))
			TaxoAssignmentRangeTable	= np.zeros((N_UcfViruses, 0))
			for RefVirusGroup_i in range(N_RefVirusGroups):
				TaxoOfMaxSimScoreRangeList = []
				MaxSimScoreRangeList = []
				PhyloStatRangeList = []
				TaxoAssignmentRangeList = []
				for UcfVirus_i in range(N_UcfViruses):
					TaxoOfMaxSimScoreDist_Counter = sorted(list(Counter(TaxoOfMaxSimScoreDistTable[:, UcfVirus_i, RefVirusGroup_i]).items()), key = operator.itemgetter(1), reverse = True)
					
					TaxoOfMaxSimScoreRangeList.append(", ".join(["%s: %.2f"%(Taxo, Count/float(N_Bootstrap)) for Taxo, Count in TaxoOfMaxSimScoreDist_Counter]))
					MaxSimScoreRangeList.append(", ".join(["%s: %s"%(Taxo, "-".join(map(str,np.around(hpd(x = MaxSimScoreDistTable[np.where(TaxoOfMaxSimScoreDistTable[:, UcfVirus_i, RefVirusGroup_i]==Taxo)[0], UcfVirus_i, RefVirusGroup_i], alpha = 0.05),3)))) for Taxo in zip(*TaxoOfMaxSimScoreDist_Counter)[0]]))
					PhyloStatRangeList.append(", ".join(["%s: %s"%(Taxo, ",".join(["p(%s): %s"%(x[0], x[1]/float(N_Bootstrap)) for x in sorted(list(Counter(PhyloStatDistTable[np.where(TaxoOfMaxSimScoreDistTable[:, UcfVirus_i, RefVirusGroup_i]==Taxo)[0], UcfVirus_i, RefVirusGroup_i]).items()), key = operator.itemgetter(1), reverse = True)])) for Taxo in zip(*TaxoOfMaxSimScoreDist_Counter)[0]]))
					TaxoAssignmentRangeList.append(", ".join(["%s: %.2f"%(Taxo, Count/float(N_Bootstrap)) for Taxo, Count in sorted(list(Counter(TaxoAssignmentDistTable[:, UcfVirus_i, RefVirusGroup_i]).items()), key = operator.itemgetter(1), reverse = True)]))
				
				TaxoOfMaxSimScoreRangeTable = np.column_stack((TaxoOfMaxSimScoreRangeTable, TaxoOfMaxSimScoreRangeList))
				MaxSimScoreRangeTable = np.column_stack((MaxSimScoreRangeTable, MaxSimScoreRangeList))
				PhyloStatRangeTable = np.column_stack((PhyloStatRangeTable, PhyloStatRangeList))
				TaxoAssignmentRangeTable = np.column_stack((TaxoAssignmentRangeTable, TaxoAssignmentRangeList))
			
			FinalisedTaxoAssignmentRangeList = []
			for UcfVirus_i in range(N_UcfViruses):
				BootstrappedFinalisedTaxoAssignmentList	= []
				for Bootstrap_i in range(0, N_Bootstrap):
					MaxSimScore_i		= np.argmax(MaxSimScoreDistTable[Bootstrap_i][UcfVirus_i])
					FinalisedTaxoAssignment	= TaxoAssignmentDistTable[Bootstrap_i][UcfVirus_i][MaxSimScore_i]
					if FinalisedTaxoAssignment != "Unclassified":
						BootstrappedFinalisedTaxoAssignmentList.append("%s (%s)"%(FinalisedTaxoAssignment, ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]))
					else:
						if MaxSimScoreDistTable[Bootstrap_i][UcfVirus_i][MaxSimScore_i] <= DatabaseAssignmentSimilarityScore_Cutoff:
							BootstrappedFinalisedTaxoAssignmentList.append("Unclassified (NA)")
						else:
							BootstrappedFinalisedTaxoAssignmentList.append("Unclassified (%s)"%(ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]))
				
				FinalisedTaxoAssignmentRangeList.append(", ".join(["%s: %.2f"%(Taxo, Count/float(N_Bootstrap)) for Taxo, Count in sorted(list(Counter(BootstrappedFinalisedTaxoAssignmentList).items()), key = operator.itemgetter(1), reverse = True)]))

		print("- Write the results to ClassificationResults.txt")

		if Bootstrap == True:
			try:
				np.savetxt(	fname = ClassificationResultFile,
						X = np.column_stack((	list(map(', '.join, self.ucf_genomes["SeqIDLists"])),
									self.ucf_genomes["VirusNameList"],
									[", ".join(["%s (%s)"%(Taxo, Range) for Taxo, Range in zip(TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange[0], TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange[1])]) for TaxoOfMaxSimScore_TaxoOfMaxSimScoreRange in zip(TaxoOfMaxSimScoreTable, TaxoOfMaxSimScoreRangeTable)],
									[", ".join(["%s (%s)"%(Score, Range) for Score, Range in zip(MaxSimScore_MaxSimScoreRange[0], MaxSimScore_MaxSimScoreRange[1])]) for MaxSimScore_MaxSimScoreRange in zip(np.around(MaxSimScoreTable,3).astype("str"), MaxSimScoreRangeTable)],
									[", ".join(["%s (%s)"%(Stat, Range) for Stat, Range in zip(PhyloStat_PhyloStatRange[0], PhyloStat_PhyloStatRange[1])]) for PhyloStat_PhyloStatRange in zip(PhyloStatTable, PhyloStatRangeTable)],
									[", ".join(["%s (%s)"%(Taxo, Range) for Taxo, Range in zip(TaxoAssignment_TaxoAssignmentRange[0], TaxoAssignment_TaxoAssignmentRange[1])]) for TaxoAssignment_TaxoAssignmentRange in zip(TaxoAssignmentTable, TaxoAssignmentRangeTable)],
									["%s (%s)"%(FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange[0], FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange[1]) for FinalisedTaxoAssignment_FinalisedTaxoAssignmentRange in zip(FinalisedTaxoAssignmentList, FinalisedTaxoAssignmentRangeList)],
									FinalisedVirusGroupingList,
									)),
						fmt = '%s',
						delimiter = "\t",
						header = "Sequence identifier\tSequence description\tCandidate class (class of the best match reference virus)\tSimilarity score*\tSupport from dendrogram**\tEvaluated taxonomic assignment\tThe best taxonomic assignment\tProvisional virus taxonomy")
			except Exception as ex1:
				print(f"ex1 {ex1}")
				breakpoint()
			try:
				with open(ClassificationResultFile, "a") as ClassificationResult_txt:
					ClassificationResult_txt.write(	"\n"+
									"*Similarity score cutoff\n"+
									"\n".join(["%s, %s: %.4f (%s)"%(RefVirusGroup, TaxoGrouping, d["CutOff"], "-".join(np.around(hpd(np.array([PairwiseSimilarityScore_CutoffDist_Dict[BootStrap_i][RefVirusGroup][TaxoGrouping]["CutOff"] for BootStrap_i in range(N_Bootstrap)])),3).astype("str"))) for RefVirusGroup, D in PairwiseSimilarityScore_Cutoff_Dict.items() for TaxoGrouping, d in D.items()])+
									"\n\n"+
									"**Support from dendrogram\n"+
									"NA:the sequence is not similar enough to any of the reference sequences to be assigned to any classes\n"
									"1: the sequence is embedded within a clade of the candidate class\n"
									"2: the sequence has a sister relationship with the candidate class and they are similar enough (the 1st criterion)\n"
									"3: the sequence is 'sandwished' between 2 branches of the candidate class\n"
									"4: the sequence has a paraphyletic relationship with the candidate class (just inside)\n"
									"5: the sequence has a paraphyletic relationship with the candidate class (just outside)\n"
									"6: the candidate class is not supported by the dendrogram\n"
									)
			except Exception as ex2:
				print(f"ex2 {ex2}")
				breakpoint()

		try:
			np.savetxt(	fname = ClassificationResultFile,
					X = np.column_stack((	list(map(', '.join, self.ucf_genomes["SeqIDLists"])),
								self.ucf_genomes["VirusNameList"],
								list(map(', '.join, TaxoOfMaxSimScoreTable)),
								list(map(', '.join, np.around(MaxSimScoreTable,3).astype("str"))),
								list(map(', '.join, PhyloStatTable)),
								list(map(', '.join, TaxoAssignmentTable)),
								FinalisedTaxoAssignmentList,
								FinalisedVirusGroupingList,
								)),
					fmt = '%s',
					delimiter = "\t",
					header = "Sequence identifier\tSequence description\tCandidate class (class of the best match reference virus)\tSimilarity score*\tSupport from dendrogram**\tEvaluated taxonomic assignment\tThe best taxonomic assignment\tProvisional virus taxonomy")
			
			with open(ClassificationResultFile, "a") as ClassificationResult_txt:
				ClassificationResult_txt.write(	"\n"+
								"*Similarity score cutoff\n"+
								"\n".join(["%s, %s: %.4f"%(RefVirusGroup, TaxoGrouping, d["CutOff"]) for RefVirusGroup, D in PairwiseSimilarityScore_Cutoff_Dict.items() for TaxoGrouping, d in D.items()])+
								"\n\n"+
								"**Support from dendrogram\n"+
								"NA:the sequence is not similar enough to any of the reference sequences to be assigned to any classes\n"
								"1: the sequence is embedded within a clade of the candidate class\n"
								"2: the sequence has a sister relationship with the candidate class and they are similar enough (the 1st criterion)\n"
								"3: the sequence is 'sandwished' between 2 branches of the candidate class\n"
								"4: the sequence has a paraphyletic relationship with the candidate class (just inside)\n"
								"5: the sequence has a paraphyletic relationship with the candidate class (just outside)\n"
								"6: the candidate class is not supported by the dendrogram\n"
								)
		except Exception as ex3:
			print(f"ex3 {ex3}")
		
		try:
			if IncludeIncompleteGenomes_UcfVirus == True:
				if IncludeIncompleteGenomes_RefVirus == True:
					################################################################################
					print("- Save variables to VirusClassificationAndEvaluation.AllUcfGenomes.AllRefGenomes.shelve")
					################################################################################
					VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.AllUcfGenomes.AllRefGenomes.shelve"
				elif IncludeIncompleteGenomes_RefVirus == False:
					################################################################################
					print("- Save variables to VirusClassificationAndEvaluation.AllUcfGenomes.CompleteRefGenomes.shelve")
					################################################################################
					VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.AllUcfGenomes.CompleteRefGenomes.shelve"
			elif IncludeIncompleteGenomes_UcfVirus == False:
				if IncludeIncompleteGenomes_RefVirus == True:
					################################################################################
					print("- Save variables to VirusClassificationAndEvaluation.CompleteUcfGenomes.AllRefGenomes.shelve")
					################################################################################
					VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.CompleteUcfGenomes.AllRefGenomes.shelve"
				elif IncludeIncompleteGenomes_RefVirus == False:
					################################################################################
					print("- Save variables to VirusClassificationAndEvaluation.CompleteUcfGenomes.CompleteRefGenomes.shelve")
					################################################################################
					VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.CompleteUcfGenomes.CompleteRefGenomes.shelve"
		except Exception as ex4:
			print(f"ex3 {ex4}")

		try:

			#VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.shelve"
			Parameters = shelve.open(VariableShelveFile_UcfVirus,"n")
			for key in [	
					"TaxoAssignmentTable",
					"MaxSimScoreTable",
					"TaxoOfMaxSimScoreTable",
					"PhyloStatTable",
					"FinalisedTaxoAssignmentList",
					"FinalisedVirusGroupingList",
					"PairwiseSimilarityScore_Cutoff_Dict",
					
					"SeqIDLists_RemainedUcfVirus",
					"VirusNameList_RemainedUcfVirus",
					
					"MaxSimScoreDistTable",
					"TaxoOfMaxSimScoreDistTable",
					"TaxoAssignmentDistTable",
					"PhyloStatDistTable",
					"PairwiseSimilarityScore_CutoffDist_Dict",
					
					"TaxoOfMaxSimScoreRangeTable",
					"MaxSimScoreRangeTable",
					"PhyloStatRangeTable",
					"TaxoAssignmentRangeTable",
					"FinalisedTaxoAssignmentRangeList",
					]:
				try:
					Parameters[key] = locals()[key]
				except TypeError:
					pass

			Parameters.close()

		except Exception as ex5:
			print(f"ex4 {ex5}")
			breakpoint()
		


