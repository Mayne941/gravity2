from copy import copy
import os, shelve, string, random, subprocess, pickle
from scipy.sparse import coo_matrix

#Local functions
from app.utils.pphmm_signature_table_constructor import PPHMMSignatureTable_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.download_genbank_file import DownloadGenBankFile
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_ref_virus_vars

class UcfVirusAnnotator:
	def __init__(self,
				GenomeSeqFile_UcfVirus,
				ShelveDir_UcfVirus,
				ShelveDirs_RefVirus,
				IncludeIncompleteGenomes_UcfVirus	= True,
				IncludeIncompleteGenomes_RefVirus	= False,
				SeqLength_Cutoff					= 0,
				HMMER_N_CPUs						= 20,
				HMMER_C_EValue_Cutoff				= 1E-3,
				HMMER_HitScore_Cutoff				= 0,
			) -> None:
		'''fpaths'''
		self.GenomeSeqFile_UcfVirus = GenomeSeqFile_UcfVirus
		self.ShelveDir_UcfVirus 	= ShelveDir_UcfVirus
		self.ShelveDirs_RefVirus	= ShelveDirs_RefVirus.split(", ")
		self.VariableShelveDir_UcfVirus			= f"{ShelveDir_UcfVirus}/Shelves"
		self.IncludeIncompleteGenomes_UcfVirus = IncludeIncompleteGenomes_UcfVirus
		self.IncludeIncompleteGenomes_RefVirus = IncludeIncompleteGenomes_RefVirus
		self.SeqLength_Cutoff = SeqLength_Cutoff
		self.HMMER_N_CPUs = HMMER_N_CPUs
		self.HMMER_C_EValue_Cutoff = HMMER_C_EValue_Cutoff
		self.HMMER_HitScore_Cutoff = HMMER_HitScore_Cutoff
		self.VariableShelveFile_UcfVirus = self.VariableShelveDir_RefVirus = {}


	def get_genbank(self) -> None:
		'''2/2: If not present, download GenBank files'''
		print("- Download GenBank file\n GenomeSeqFile doesn't exist. GRAViTy is downloading the GenBank file(s).")
		DownloadGenBankFile(self.GenomeSeqFile_UcfVirus, self.VariableShelveFile_UcfVirus["SeqIDLists"])

	def mkdirs(self, ShelveDir_RefVirus):
		'''Make filepaths for retrieval and storage'''
		HMMERDir_RefVirus			= f"{ShelveDir_RefVirus}/HMMER"
		HMMER_PPHMMDBDir_RefVirus	= f"{HMMERDir_RefVirus}/HMMER_PPHMMDb"
		HMMER_PPHMMDB_RefVirus		= f"{HMMER_PPHMMDBDir_RefVirus}/HMMER_PPHMMDb"
		return HMMERDir_RefVirus, HMMER_PPHMMDB_RefVirus

	def annotate(self) -> None:
		'''3/3: Annotate unclassified viruses via PPHMM and GOM, generate loc/sig tables for unclassified, save to pickle'''
		PPHMMSignatureTable_Dict = PPHMMLocationTable_Dict = GOMSignatureTable_Dict = {}
		N_RefVirusGroups, RefVirusGroup_i = len(self.ShelveDirs_RefVirus), 0
		for ShelveDir_RefVirus in self.ShelveDirs_RefVirus:
			RefVirusGroup_i = RefVirusGroup_i+1
			RefVirusGroup = ShelveDir_RefVirus.split("/")[-1]
			print(f"Annotate unclassified viruses using the PPHMM and GOM databases of the reference viruses: {RefVirusGroup} ({RefVirusGroup_i}/{N_RefVirusGroups})")

			'''Get fpaths'''
			HMMERDir_RefVirus, HMMER_PPHMMDB_RefVirus = self.mkdirs(ShelveDir_RefVirus)

			'''Retrieve variables related to the reference viruses'''
			Parameters = retrieve_ref_virus_vars(f"{ShelveDir_RefVirus}/Shelves", self.IncludeIncompleteGenomes_RefVirus)
			if "GOMDB_coo" in Parameters.keys(): 
				Parameters["GOMDB"] = {GOMID:GOM_coo.toarray() for GOMID, GOM_coo in Parameters["GOMDB_coo"].items()}

			'''Generate PPHMMSignatureTable, PPHMMLocationTable, and GOMSignatureTable for unclassified viruses using the reference PPHMM and GOM database'''
			'''Make HMMER_hmmscanDir'''
			HMMER_hmmscanDir_RefVirus = f"{HMMERDir_RefVirus}/hmmscan_"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMER_hmmscanDir_RefVirus)
			
			'''Generate PPHMMSignatureTable and PPHMMLocationTable'''
			PPHMMSignatureTable, PPHMMLocationTable = PPHMMSignatureTable_Constructor(	
															SeqIDLists				= self.VariableShelveFile_UcfVirus["SeqIDLists"],
															GenBankFile				= self.GenomeSeqFile_UcfVirus,
															TranslTableList			= self.VariableShelveFile_UcfVirus["TranslTableList"],
															SeqLength_Cutoff		= self.SeqLength_Cutoff,
															HMMER_PPHMMDB			= HMMER_PPHMMDB_RefVirus,
															HMMER_hmmscanDir		= HMMER_hmmscanDir_RefVirus,
															HMMER_N_CPUs			= self.HMMER_N_CPUs,
															HMMER_C_EValue_Cutoff	= self.HMMER_C_EValue_Cutoff,
															HMMER_HitScore_Cutoff	= self.HMMER_HitScore_Cutoff
															)
			
			'''Delete HMMER_hmmscanDir'''
			_ = subprocess.call("rm -rf %s" %HMMER_hmmscanDir_RefVirus, shell = True)	
			
			GOMSignatureTable = GOMSignatureTable_Constructor(PPHMMLocationTable = PPHMMLocationTable,
																GOMDB			 = Parameters["GOMDB"],
																GOMIDList		 = Parameters["GOMIDList"])
			
			PPHMMSignatureTable_Dict[RefVirusGroup] = PPHMMSignatureTable
			PPHMMLocationTable_Dict[RefVirusGroup]  = PPHMMLocationTable
			GOMSignatureTable_Dict[RefVirusGroup]   = GOMSignatureTable
		
		'''Generate out fpaths'''
		if self.IncludeIncompleteGenomes_UcfVirus == True:
			if self.IncludeIncompleteGenomes_RefVirus == True:
				VariableShelveFile_UcfVirus = f"{self.VariableShelveDir_UcfVirus}/UcfVirusAnnotator.AllUcfGenomes.AllRefGenomes.p"

			elif self.IncludeIncompleteGenomes_RefVirus == False:
				VariableShelveFile_UcfVirus = f"{self.VariableShelveDir_UcfVirus}/UcfVirusAnnotator.AllUcfGenomes.CompleteRefGenomes.p"

		elif self.IncludeIncompleteGenomes_UcfVirus == False:
			if self.IncludeIncompleteGenomes_RefVirus == True:
				VariableShelveFile_UcfVirus = f"{self.VariableShelveDir_UcfVirus}/UcfVirusAnnotator.CompleteUcfGenomes.AllRefGenomes.p"

			elif self.IncludeIncompleteGenomes_RefVirus == False:
				VariableShelveFile_UcfVirus = f"{self.VariableShelveDir_UcfVirus}/UcfVirusAnnotator.CompleteUcfGenomes.CompleteRefGenomes.p"
		
		'''Construct out dict, dump to pickle'''
		all_ucf_genomes = {} 
		all_ucf_genomes["PPHMMSignatureTable_Dict_coo"] = {RefVirusGroup:coo_matrix(PPHMMSignatureTable) for RefVirusGroup, PPHMMSignatureTable in PPHMMSignatureTable_Dict.items()}
		all_ucf_genomes["PPHMMLocationTable_Dict_coo"]  = {RefVirusGroup:coo_matrix(PPHMMLocationTable) for RefVirusGroup, PPHMMLocationTable in PPHMMLocationTable_Dict.items()}
		all_ucf_genomes["GOMSignatureTable_Dict"] 		= GOMSignatureTable_Dict

		pickle.dump(all_ucf_genomes, open(VariableShelveFile_UcfVirus, "wb"))

	def main(self) -> None:
		'''Generate PPHMM signature table, PPHMM location table, and GOM signature table for unclassified viruses, using PPHMM databases of reference viruses'''
		section_header("#Generate PPHMM signature table, PPHMM location table, and GOM signature table\n#for unclassified viruses, using PPHMM databases of reference viruses")
				
		'''1/3: Retrieve variables'''
		self.VariableShelveFile_UcfVirus = retrieve_genome_vars(self.VariableShelveDir_UcfVirus, self.IncludeIncompleteGenomes_UcfVirus)

		if not os.path.isfile(self.GenomeSeqFile_UcfVirus):
			'''(OPT) 2/3: Get GenBank if not present'''
			self.get_genbank()

		'''3/3: Annotate unclassified viruses using PPHMM & GOM'''
		self.annotate()


