from app.utils.console_messages import section_header

import numpy as np
from collections import Counter
import os, shelve
import re
'''
Read genome description table (Virus Metadata Resource -- VMR)
---------------------------------------------
The genome description table should be formatted as follows
[Col idx]	= [header]
0		= Realm
1		= Subrealm
2		= Kingdom
3		= Subkingdom
4		= Phylum
5		= Subphylum
6		= Class
7		= Subclass
8		= Order
9		= Suborder
10		= Family
11		= Subfamily
12		= Genus
13		= Subgenus
14		= Species
15		= Exemplar or additional isolate
16		= Virus name (s)
17		= Virus name abbreviation (s)
18		= Virus isolate designation
19		= Virus GENBANK accession
20		= Virus REFSEQ accession
21		= Virus sequence complete
22		= Sort
23		= Study Group approved
24		= Genome Composition
25		= Genetic code table
26		= Baltimore Group
27		= Database
28		= Taxonomic grouping
'''

class ReadGenomeDescTable:
	'''Parse, transform and load input VMR to GRAViTy compatible format, store as shelve'''
	def __init__(self, 
			GenomeDescTableFile, 
			ShelveDir, 
			Database = None, 
			Database_Header = None, 
			TaxoGrouping_Header	= "Family", 
			TaxoGroupingFile = None) -> None:
		self.GenomeDescTableFile = GenomeDescTableFile
		self.VariableShelveDir = ShelveDir + "/Shelves"
		if not os.path.exists(self.VariableShelveDir):
			os.makedirs(self.VariableShelveDir)
		self.Database = Database
		self.Database_Header = Database_Header
		self.TaxoGrouping_Header = TaxoGrouping_Header
		self.TaxoGroupingFile = TaxoGroupingFile
		self.BaltimoreList, self.OrderList, self.FamilyList, \
			self.SubFamList, self.GenusList, self.VirusNameList, \
			self.SeqIDLists, self.SeqStatusList, self.TaxoGroupingList, \
			self.TranslTableList, self.DatabaseList, self.VirusIndexList = [], [], [], [],\
																		   [], [], [], [],\
																		   [], [], [], []
																		   
	def open_table(self) -> None:
		'''Open VMR file and read in relevant data to volatile'''
		print("- Read the GenomeDesc table")
		# RM < Replace with pandas? Existing code is quick and works well...
		with open(self.GenomeDescTableFile, "r") as GenomeDescTable_txt:
			header = next(GenomeDescTable_txt)
			header = header.split("\r\n")[0].split("\n")[0].split("\t")
			
			# RM < new VMR format is not consistent with this
			Baltimore_i		= header.index("Baltimore Group")
			Order_i			= header.index("Order")
			Family_i		= header.index("Family")
			Subfamily_i		= header.index("Subfamily")
			Genus_i			= header.index("Genus")
			VirusName_i		= header.index("Virus name (s)")
			SeqID_i			= header.index("Virus GENBANK accession")
			SeqStatus_i		= header.index("Virus sequence complete")
			TranslTable_i	= header.index("Genetic code table")
			
			if self.Database_Header != None: 
				Database_i	= header.index(self.Database_Header)
			
			if self.TaxoGrouping_Header != None:
				TaxoGrouping_i	= header.index(self.TaxoGrouping_Header)
			else:
				TaxoGrouping_i	= header.index("Family")
			
			for Virus_i, Line in enumerate(GenomeDescTable_txt):
				Line = Line.split("\r\n")[0].split("\n")[0].split("\t")
				SeqIDList = re.findall(r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{4}[0-9]{6,8}|[A-Z]{2}_[0-9]{6}",Line[SeqID_i])

				if SeqIDList != [] or Line[SeqID_i]!="":
					'''Ignore record without sequences'''
					self.VirusIndexList.append(Virus_i)
					self.BaltimoreList.append(Line[Baltimore_i])
					self.OrderList.append(Line[Order_i])
					self.FamilyList.append(Line[Family_i])
					self.SubFamList.append(Line[Subfamily_i])
					self.GenusList.append(Line[Genus_i])
					self.VirusNameList.append(re.sub(r"^\/|\/$","",re.sub(r"[\/ ]{2,}","/",re.sub(r"[^\w^ ^\.^\-]+","/",re.sub(r"[ ]{2,}"," ",Line[VirusName_i]))))) #clean the virus name
					self.TaxoGroupingList.append(Line[TaxoGrouping_i])
					self.SeqStatusList.append(Line[SeqStatus_i])
					
					if SeqIDList != []:
						self.SeqIDLists.append(SeqIDList)
					else:
						self.SeqIDLists.append([Line[SeqID_i]])
					
					try:
						self.TranslTableList.append(int(Line[TranslTable_i]))
					except ValueError:
						print(("Genetic code is not specified. GRAViTy will use the standard code for %s"%self.VirusNameList[-1]))
						self.TranslTableList.append(1)
					if self.Database_Header != None: 
						self.DatabaseList.append(Line[Database_i])

	def update_desc_table(self, flag) -> dict:
		'''Create dictionary in GRAViTy structure for saving to persistent storage'''
		master_data = {}
		if flag == "all_genomes":
			IncludedGenomes_IndexList = np.where(self.DatabaseList==self.Database)[0]
			master_data["DatabaseList"]    = self.DatabaseList[IncludedGenomes_IndexList]
		else:
			IncludedGenomes_IndexList = [SeqStatus_i for SeqStatus_i, SeqStatus in enumerate(self.SeqStatusList) if "Complete" in SeqStatus]
			if self.Database != None: 
				master_data["DatabaseList"] = self.DatabaseList[IncludedGenomes_IndexList]

		self.BaltimoreList = master_data["BaltimoreList"]		= self.BaltimoreList[IncludedGenomes_IndexList]
		self.OrderList 	   = master_data["OrderList"]			= self.OrderList[IncludedGenomes_IndexList]
		self.FamilyList    = master_data["FamilyList"]			= self.FamilyList[IncludedGenomes_IndexList]
		self.SubFamList    = master_data["SubFamList"]			= self.SubFamList[IncludedGenomes_IndexList]
		self.GenusList 	   = master_data["GenusList"]			= self.GenusList[IncludedGenomes_IndexList]
		self.VirusNameList = master_data["VirusNameList"]		= self.VirusNameList[IncludedGenomes_IndexList]		
		self.SeqIDLists    = master_data["SeqIDLists"]			= self.SeqIDLists[IncludedGenomes_IndexList]
		self.SeqStatusList = master_data["SeqStatusList"]		= self.SeqStatusList[IncludedGenomes_IndexList]
		self.TranslTableList = master_data["TranslTableList"]		= self.TranslTableList[IncludedGenomes_IndexList]
		self.TaxoGroupingList = master_data["TaxoGroupingList"]		= self.TaxoGroupingList[IncludedGenomes_IndexList]
		self.DatabaseList  = master_data["DatabaseList"]			= self.DatabaseList
		return master_data

	def save_desc_table(self, complete_desc_table, fname):
		'''Save dictionary in GRAViTy structure to persistent storage'''
		# RM < Shelve > pickle
		parameters = shelve.open(self.VariableShelveDir + fname, "n")
		for key in ["BaltimoreList",
				"OrderList",
				"FamilyList",
				"SubFamList",
				"GenusList",
				"VirusNameList",
				"SeqIDLists",
				"SeqStatusList",
				"TaxoGroupingList",
				"TranslTableList",
				"DatabaseList",
				]:
			try:
				parameters[key] = complete_desc_table[key]
				print("\t" + key)
			except TypeError:
				pass
		parameters.close()

	def entrypoint(self) -> None:
		'''RM < UPDATE DOCSTRING'''
		section_header("Read the GenomeDesc table")

		'''Open & parse input VMR (txt) file'''
		self.open_table()
					
		if self.TaxoGroupingFile != None:
			'''PL1 has option to specify taxonomic grouping level.'''
			# RM < NEEDS TESTING
			self.TaxoGroupingList = []
			with open(self.TaxoGroupingFile, "r") as TaxoGrouping_txt:
				for TaxoGrouping in TaxoGrouping_txt:
					TaxoGrouping = TaxoGrouping.split("\r\n")[0].split("\n")[0]
					self.TaxoGroupingList.append(TaxoGrouping)
			self.TaxoGroupingList = np.array(self.TaxoGroupingList)
			self.TaxoGroupingList = self.TaxoGroupingList[self.VirusIndexList]
		
		'''Generate rows for output'''
		all_desc_table = {}
		self.BaltimoreList  = all_desc_table["BaltimoreList"]	= np.array(self.BaltimoreList)
		self.OrderList		= all_desc_table["OrderList"]		= np.array(self.OrderList)
		self.FamilyList		= all_desc_table["FamilyList"] 		= np.array(self.FamilyList)
		self.SubFamList		= all_desc_table["SubFamList"] 		= np.array(self.SubFamList)
		self.GenusList		= all_desc_table["GenusList"] 		= np.array(self.GenusList)
		self.VirusNameList	= all_desc_table["VirusNameList"]	= np.array(self.VirusNameList)
		
		self.SeqIDLists		.extend([[1],[1,2]])	#making sure that SeqIDLists will be a h list (and not a v list)
		self.SeqIDLists		= np.array(self.SeqIDLists)
		self.SeqIDLists		= self.SeqIDLists[:-2]
		
		# RM < Do a list comp dedupe and not raise exception. Break out to new fn
		SeqIDFlatList 	= [SeqID for SeqIDList in self.SeqIDLists for SeqID in SeqIDList]
		if len(SeqIDFlatList) != len(set(SeqIDFlatList)):
			'''Check if there exist duplicated accession numbers, fail IF TRUE'''
			print("The following accession numbers appear more than once: ")
			print("\n".join([SeqID for SeqID, count in Counter(SeqIDFlatList).items() if count > 1]))
			raise SystemExit("GRAViTy terminated.")

		all_desc_table["SeqIDLists"] = self.SeqIDLists
		
		self.SeqStatusList	  = all_desc_table["SeqStatusList"]		 = np.array(self.SeqStatusList)
		self.TaxoGroupingList = all_desc_table["TaxoGroupingList"] 	 = np.array(self.TaxoGroupingList)
		self.TranslTableList  = all_desc_table["TranslTableList"]	 = np.array(self.TranslTableList)
		self.DatabaseList	  = all_desc_table["DatabaseList"]		 = np.array(self.DatabaseList)

		'''Create ReadGenomeDescTable "all genomes' db'''
		print("- Save variables to ReadGenomeDescTable.AllGenomes.shelve")
		if self.Database != None:
			'''If a database is provided, filter to these specific entries'''
			all_desc_table = self.update_desc_table("all_genomes")

		self.save_desc_table(all_desc_table, "/ReadGenomeDescTable.AllGenomes.shelve")
	
		'''Create ReadGenomeDescTable "complete genomes' db'''
		print("- Save variables to ReadGenomeDescTable.CompleteGenomes.shelve")
		complete_desc_table = self.update_desc_table("complete_genomes")
		self.save_desc_table(complete_desc_table, "/ReadGenomeDescTable.CompleteGenomes.shelve")

