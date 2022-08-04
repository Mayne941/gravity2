from sklearn.feature_selection import mutual_info_classif
import numpy as np
import shelve, os, sys
import pickle
#Local functions
from app.utils.ordered_set import OrderedSet

def MutualInformationCalculator (
	ShelveDir,
	IncludeIncompleteGenomes= False,
	
	VirusGroupingFile	= None,
	
	N_Sampling		= 100, 
	SamplingStrategy	= "balance_with_repeat", #None, "balance_without_repeat", "balance_with_repeat"
	SampleSizePerGroup	= 10,
	
	):
	print("################################################################################")
	print("#Compute mutual information score                                              #")
	print("################################################################################")
	'''
	Compute mutual information scores
	---------------------------------------------
	'''
	################################################################################
	print("- Define dir/file paths")
	################################################################################
	print("\tto program output shelve")
	#-------------------------------------------------------------------------------
	VariableShelveDir		= ShelveDir+"/Shelves"

	print("\t\tto Mutual information score directory")
	#-------------------------------------------------------------------------------
	MutualInformationScoreDir	= VariableShelveDir+"/MutualInformationScore"
	if not os.path.exists(MutualInformationScoreDir):
		os.makedirs(MutualInformationScoreDir)
	
	################################################################################
	print("- Retrieve variables")
	################################################################################
	if IncludeIncompleteGenomes == True:
		print("\tfrom ReadGenomeDescTable.AllGenomes.shelve")
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		print("\tfrom ReadGenomeDescTable.CompleteGenomes.shelve")
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.CompleteGenomes.shelve"
	
	#VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.shelve"
	Parameters = shelve.open(VariableShelveFile)
	for key in [
			"SeqIDLists",
			"VirusNameList",
			
			"BaltimoreList",
			"OrderList",
			"FamilyList",
			"SubFamList",
			"GenusList",
			"TaxoGroupingList",
			]:
		globals()[key] = Parameters[key]
		print("\t\t" + key)
	
	Parameters.close()
	
	if IncludeIncompleteGenomes == True:
		print("\tfrom RefVirusAnnotator.AllGenomes.shelve")
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		print("\tfrom RefVirusAnnotator.CompleteGenomes.shelve")
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.CompleteGenomes.shelve"
	
	#VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.shelve"
	#Parameters = shelve.open(VariableShelveFile) # RM <<
	Parameters = pickle.load(open(f"{VariableShelveDir}/RefVirusAnnotator.CompleteGenomes.p", "rb"))
	for key in  [	
			"PPHMMSignatureTable",
			"PPHMMSignatureTable_coo",
			]:
		try:
			globals()[key] = Parameters[key]
			print("\t\t"+key)
		except KeyError:
			pass
	
	#Parameters.close() #RM <<

	if "PPHMMSignatureTable_coo" in list(globals().keys()):	globals()["PPHMMSignatureTable"] = PPHMMSignatureTable_coo.toarray()
	
	print("\tfrom PPHMMDBConstruction.shelve")
	#-------------------------------------------------------------------------------
	VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
	Parameters = shelve.open(VariableShelveFile)
	for key in [	
			"ClusterDescList",
			]:
		globals()[key] = Parameters[key]
		print("\t\t" + key)
	
	Parameters.close()
	
	PPHMMDesc = ["PPHMM|"+ClusterDesc for ClusterDesc in ClusterDescList.astype("str")]
	PPHMMDesc = np.array(PPHMMDesc)
	
	if VirusGroupingFile == None:
		VirusGroupingDict = {"Overvall":TaxoGroupingList}
	else:
		################################################################################
		print("- Read virus grouping file")
		################################################################################
		VirusGroupingTable = []
		with open(VirusGroupingFile, "r") as VirusGrouping_txt:
			VirusGroupingSchemeList = next(VirusGrouping_txt)	#the header
			VirusGroupingSchemeList = VirusGroupingSchemeList.split("\r\n")[0].split("\n")[0].split("\t")
			for Line in VirusGrouping_txt:
				Line = Line.split("\r\n")[0].split("\n")[0].split("\t")
				VirusGroupingTable.append(Line)
		
		VirusGroupingTable = np.array(VirusGroupingTable).T
		VirusGroupingDict = {VirusGroupingScheme:np.array(VirusGroupingList) for VirusGroupingScheme, VirusGroupingList in zip(VirusGroupingSchemeList, VirusGroupingTable)}
		VirusGroupingDict = {"Overvall":TaxoGroupingList}
	
	################################################################################
	print("- Compute mutual information between PPHMM scores and virus classification")
	################################################################################
	ResultDict		= {}
	N_PPHMMs		= len(PPHMMDesc)
	N_VirusGroupingSchemes	= len(VirusGroupingDict)
	VirusGroupingScheme_i	= 1.0
	for VirusGroupingScheme, VirusGroupingList in VirusGroupingDict.items():
		#Compute mutual information between PPHMM scores and virus classification
		#-------------------------------------------------------------------------------
		ResultDict[VirusGroupingScheme]	= {}
		VirusGroupingIDList		= [VirusGroupingID for VirusGroupingID in OrderedSet(VirusGroupingList) if VirusGroupingID!="-"]

		IngroupVirus_IndexList		= np.where([(VirusGroupingID != "-" and not VirusGroupingID.startswith("O_")) for VirusGroupingID in VirusGroupingList.astype("str")])[0]
		IngroupVirus_PPHMM_IndexList	= np.where(np.sum(PPHMMSignatureTable[IngroupVirus_IndexList] != 0, axis = 0) != 0)[0]
		IngroupVirus_PPHMMDesc		= PPHMMDesc[IngroupVirus_PPHMM_IndexList]
		
		PPHMMSignatureTable_Subset	= PPHMMSignatureTable[:,IngroupVirus_PPHMM_IndexList]
		
		SampledMutualInformationTable = []
		for Sampling_i in range(N_Sampling):
			if   SamplingStrategy == None:
				SampledVirus_IndexList = np.where([VirusGroupingID != "-" for VirusGroupingID in VirusGroupingList])[0]
				
			elif SamplingStrategy == "balance_without_repeat":
				SampledVirus_IndexList = sum([list(np.random.permutation(np.where(VirusGroupingList == VirusGroupingID)[0])[:SampleSizePerGroup]) for VirusGroupingID in VirusGroupingIDList],[])
				
			elif SamplingStrategy == "balance_with_repeat":
				SampledVirus_IndexList = sum([list(np.random.choice(a = np.where(VirusGroupingList == VirusGroupingID)[0], size = SampleSizePerGroup, replace = True)) for VirusGroupingID in VirusGroupingIDList],[])
			
			SampledMutualInformationTable.append(mutual_info_classif(PPHMMSignatureTable_Subset[SampledVirus_IndexList], VirusGroupingList[SampledVirus_IndexList]))
			
			#Progress bar
			sys.stdout.write("\033[K" + "Virus grouping scheme = %s (%d/%d): [%-20s] %d/%d samplings" % (VirusGroupingScheme, VirusGroupingScheme_i, N_VirusGroupingSchemes, '='*int(float(Sampling_i+1)/N_Sampling*20), Sampling_i+1, N_Sampling) + "\r")
			sys.stdout.flush()
		
		sys.stdout.write("\033[K")
		sys.stdout.flush()
		
		SampledMutualInformationTable = np.array(SampledMutualInformationTable)
		AverageMutualInformationScoreList = np.mean(SampledMutualInformationTable, axis = 0)
		
		VirusGroupingScheme_i = VirusGroupingScheme_i + 1.0
		
		#Sort the PPHMMs according to mutual information scores
		#-------------------------------------------------------------------------------
		PPHMMOrder = list([zip(*sorted(	[(AverageMutualInformationScore_i, AverageMutualInformationScore) for AverageMutualInformationScore_i, AverageMutualInformationScore in enumerate(AverageMutualInformationScoreList)],
						key = lambda x: x[1],
						reverse = True,
						)
					)][0]
				)

		#ResultDict[VirusGroupingScheme]["AverageMutualInformationScoreList"] = AverageMutualInformationScoreList[PPHMMOrder]
		ResultDict[VirusGroupingScheme]["AverageMutualInformationScoreList"] = [x for _, x in sorted(zip(PPHMMOrder[0], AverageMutualInformationScoreList))]
		#ResultDict[VirusGroupingScheme]["PPHMMDesc"] = IngroupVirus_PPHMMDesc[PPHMMOrder]
		ResultDict[VirusGroupingScheme]["PPHMMDesc"] = [x for _, x in sorted(zip(PPHMMOrder[0], IngroupVirus_PPHMMDesc))]
		#ResultDict[VirusGroupingScheme]["PPHMMSignatureTable"] = PPHMMSignatureTable_Subset[:,PPHMMOrder]
		ResultDict[VirusGroupingScheme]["PPHMMSignatureTable"] = [x for _, x in sorted(zip(PPHMMOrder[0], PPHMMSignatureTable_Subset))]

		#Write the results to file
		#-------------------------------------------------------------------------------
		Dat	= np.column_stack((	list(map(", ".join, SeqIDLists)),
						VirusNameList,
						BaltimoreList,
						OrderList,
						FamilyList,
						SubFamList,
						GenusList,
						TaxoGroupingList,
						VirusGroupingList,
						#PPHMMSignatureTable_Subset[:,PPHMMOrder], # RM <
						sorted(zip(PPHMMOrder[0], PPHMMSignatureTable_Subset))
					))
		''' # RM < DISABLED FOR DEV
#		Dat	= np.vstack((		["Sequence ID", "Virus name", "Baltimore classification group", "Order", "Family", "Subfamily", "Genus", "TaxoGrouping", "Virus grouping"]+IngroupVirus_PPHMMDesc[PPHMMOrder].tolist(),
		Dat	= np.vstack((["Sequence ID", "Virus name", "Baltimore classification group", "Order", "Family", "Subfamily", "Genus", "TaxoGrouping", "Virus grouping"]+sorted(zip(PPHMMOrder[0], IngroupVirus_PPHMMDesc)),
						Dat,
						#["Average mutual information score"]+[""]*8+AverageMutualInformationScoreList[PPHMMOrder].astype(str).tolist(), # RM <
						["Average mutual information score"]+[""]*8+PPHMMOrder[1].astype(str).tolist(),
					))
		
		MutualInformationScoreFile = MutualInformationScoreDir+"/MIScore.Scheme=%s.txt"%VirusGroupingScheme
		np.savetxt(	fname = MutualInformationScoreFile,
				X = Dat,
				fmt = '%s',
				delimiter = "\t")
		
		with open(MutualInformationScoreFile, "a") as MutualInformationScore_txt:
			MutualInformationScore_txt.write("N_Sampling: %d\n"%N_Sampling+
							"Sampling strategy: %s\n"%SamplingStrategy+
							"Sample size per virus group: %s\n"%SampleSizePerGroup)
		'''

	if IncludeIncompleteGenomes == True:
		################################################################################
		print("- Save variables to MutualInformationCalculator.AllGenomes.shelve")
		################################################################################
		VariableShelveFile = VariableShelveDir+"/MutualInformationCalculator.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		################################################################################
		print("- Save variables to MutualInformationCalculator.CompleteGenomes.shelve")
		################################################################################
		VariableShelveFile = VariableShelveDir+"/MutualInformationCalculator.CompleteGenomes.shelve"
	
	#VariableShelveFile = VariableShelveDir+"/MutualInformationCalculator.shelve"
	Parameters = shelve.open(VariableShelveFile,"n")
	for key in ["ResultDict"]:
		try:
			Parameters[key] = locals()[key]
			print("\t" + key)
		except TypeError:
			pass
	
	Parameters.close()


