from ete3 import Tree
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import numpy as np
import shelve, subprocess, os, operator, sys, random, string, re, pickle

from app.utils.line_count import LineCount
from app.utils.dist_mat_to_tree import DistMat2Tree
from app.utils.raw_input_with_timeout import raw_input_with_timeout
from app.utils.download_genbank_file import DownloadGenBankFile
from app.utils.console_messages import section_header
from app.utils.orf_identifier import get_orf_trasl_table, no_orf_match
from app.utils.stdout_utils import clean_stdout, progress_bar, error_handler

class PPHMMDBConstruction:
	def __init__(self,
				GenomeSeqFile,
				ShelveDir,
				ProteinLength_Cutoff		= 100,
				IncludeIncompleteGenomes	= True,
				BLASTp_evalue_Cutoff		= 1E-3,
				BLASTp_PercentageIden_Cutoff	= 50,
				BLASTp_QueryCoverage_Cutoff	= 75,
				BLASTp_SubjectCoverage_Cutoff	= 75,
				BLASTp_num_alignments		= 1000000,
				BLASTp_N_CPUs			= 20,
				MUSCLE_GapOpenCost		= -3.0,
				MUSCLE_GapExtendCost		= -0.0,
				ProtClustering_MCLInflation	= 2,
				N_AlignmentMerging		= 0,
				HHsuite_evalue_Cutoff		= 1E-6,
				HHsuite_pvalue_Cutoff		= 0.05,
				HHsuite_N_CPUs			= 10,
				HHsuite_QueryCoverage_Cutoff	= 85,
				HHsuite_SubjectCoverage_Cutoff	= 85,
				PPHMMClustering_MCLInflation	= 5,
				HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging = True,
			) -> None:
		self.GenomeSeqFile = GenomeSeqFile
		self.ShelveDir = ShelveDir
		self.ProteinLength_Cutoff = ProteinLength_Cutoff
		self.IncludeIncompleteGenomes = IncludeIncompleteGenomes
		self.BLASTp_evalue_Cutoff = BLASTp_evalue_Cutoff
		self.BLASTp_PercentageIden_Cutoff = BLASTp_PercentageIden_Cutoff
		self.BLASTp_QueryCoverage_Cutoff = BLASTp_QueryCoverage_Cutoff
		self.BLASTp_SubjectCoverage_Cutoff = BLASTp_SubjectCoverage_Cutoff
		self.BLASTp_num_alignments = BLASTp_num_alignments
		self.BLASTp_N_CPUs = BLASTp_N_CPUs
		self.MUSCLE_GapOpenCost = MUSCLE_GapOpenCost
		self.MUSCLE_GapExtendCost = MUSCLE_GapExtendCost
		self.ProtClustering_MCLInflation = ProtClustering_MCLInflation
		self.N_AlignmentMerging = N_AlignmentMerging
		self.HHsuite_evalue_Cutoff = HHsuite_evalue_Cutoff
		self.HHsuite_pvalue_Cutoff = HHsuite_pvalue_Cutoff
		self.HHsuite_N_CPUs = HHsuite_N_CPUs
		self.HHsuite_QueryCoverage_Cutoff = HHsuite_QueryCoverage_Cutoff
		self.HHsuite_SubjectCoverage_Cutoff = HHsuite_SubjectCoverage_Cutoff
		self.PPHMMClustering_MCLInflation = PPHMMClustering_MCLInflation
		self.HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging = HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging
		self.orf_tranl_table = get_orf_trasl_table()
		self.orf_no_match = no_orf_match()

	def Make_HMMER_PPHMM_DB(self, 
				HMMER_PPHMMDir, 
				HMMER_PPHMMDB, 
				ClustersDir, 
				Cluster_MetaDataDict
			):
		'''RM < DOCSTRING'''
		ClusterSizeList, ClusterSizeByTaxoGroupingList, ClusterSizeByProtList, ClusterTaxoList, ClusterProtSeqIDList, ClusterDescList = [], [], [], [], [], []
		N_PPHMMs = len(Cluster_MetaDataDict)

		for PPHMM_i in range(N_PPHMMs):
			Cluster = Cluster_MetaDataDict[PPHMM_i]["Cluster"]
			DescList = Cluster_MetaDataDict[PPHMM_i]["DescList"]
			TaxoLists = Cluster_MetaDataDict[PPHMM_i]["TaxoLists"]
			
			'''Cluster annotations'''
			ClusterSizeList.append(len(Cluster))
			ClusterSizeByTaxoGroupingList.append(", ".join(["%s: %s"%(TaxoGrouping, N_ProtSeqsInTheTaxoGroup) for TaxoGrouping, N_ProtSeqsInTheTaxoGroup in sorted(list(Counter(list(zip(*TaxoLists))[-1]).items()),
																		key = operator.itemgetter(1),
																		reverse = True
																		)]
								)
							)
			ClusterSizeByProtList.append(", ".join(["%s: %s"%(Desc, N_ProtSeqsWithDesc) for Desc, N_ProtSeqsWithDesc in sorted(	list(Counter(DescList).items()),
																		key = operator.itemgetter(1),
																		reverse = True
																		)]
								)
							)
			
			(ClusterTaxo_UniqueBaltimoreGroup,
			ClusterTaxo_UniqueOrder,
			ClusterTaxo_UniqueFamily,
			ClusterTaxo_UniqueSubFam,
			ClusterTaxo_UniqueGenus,
			ClusterTaxo_UniqueVirusName,
			ClusterTaxo_UniqueTaxoGroup) = [list(set(TaxoList)) for TaxoList in zip(*TaxoLists)]
			ClusterTaxo = []

			for ClusterTaxo_UniqueTaxoLabel in [ClusterTaxo_UniqueBaltimoreGroup, ClusterTaxo_UniqueOrder, ClusterTaxo_UniqueFamily, ClusterTaxo_UniqueSubFam, ClusterTaxo_UniqueGenus]:
				if len(ClusterTaxo_UniqueTaxoLabel) == 1:
					if ClusterTaxo_UniqueTaxoLabel[0] not in ["", "Unassigned", "unassigned"]:
						ClusterTaxo = ClusterTaxo + ClusterTaxo_UniqueTaxoLabel
				else:
					ClusterTaxo = ClusterTaxo + [b"/".join(ClusterTaxo_UniqueTaxoLabel)]
					break
			
			ClusterTaxoList.append(b"; ".join(ClusterTaxo))
			ClusterProtSeqIDList.append(", ".join(Cluster))
			
			ClusterDescCount = sorted(list(Counter(DescList).items()), key = operator.itemgetter(1), reverse = True)
			for ClusterDesc in ClusterDescCount:
				ClusterDesc = ClusterDesc[0]
				if (("Hypothetical" not in ClusterDesc) and ("hypothetical" not in ClusterDesc)):
					break
			
			ClusterDescList.append("%s|%s" %(ClusterDesc, ClusterTaxoList[-1]))
			
			'''Make a PPHMM using HMMER hmmbuild'''
			AlnClusterFile = ClustersDir+"/Cluster_%s.fasta" %PPHMM_i
			HMMER_PPHMMFile = HMMER_PPHMMDir+"/PPHMM_%s.hmm" %PPHMM_i
			_ = subprocess.Popen("hmmbuild %s %s" %(HMMER_PPHMMFile, AlnClusterFile), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			
			'''Modify the DESC line in the HMM file to ClusterDesc|ClusterTaxo'''
			with open(HMMER_PPHMMFile, "r+") as HMMER_PPHMM_txt:
				contents = HMMER_PPHMM_txt.readlines()
				contents.insert(2, f"DESC  {ClusterDescList[-1]}\n")
				contents = "".join(contents)
				HMMER_PPHMM_txt.seek(0)				#Put cursor at the beginning of the file
				HMMER_PPHMM_txt.write(contents)		#Write the contents
				HMMER_PPHMM_txt.truncate()			#Delete everything after the cursor
			
			'''Progress bar'''
			sys.stdout.write("\033[K"+ "Make HMMER PPHMMs: [%-20s] %d/%d PPHHMs" % ('='*int(float(PPHMM_i + 1)/N_PPHMMs*20), PPHMM_i+1, N_PPHMMs) + "\r")
			sys.stdout.flush()
		
		sys.stdout.write("\033[K")
		sys.stdout.flush()
		
		'''Make a HMMER HMM DB with hmmpress'''
		_ = subprocess.Popen("find %s -name '*.hmm' -exec cat {} \; > %s" %(HMMER_PPHMMDir, HMMER_PPHMMDB), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		_ = subprocess.Popen("hmmpress %s" %HMMER_PPHMMDB, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		'''Make a PPHMMDBSummary file'''
		ClusterIDList = ["Cluster_%s"%Cluster_i for Cluster_i in range(N_PPHMMs)]
		np.savetxt(fname = HMMER_PPHMMDB+"_Summary.txt",
				X = np.column_stack((	ClusterIDList,
							ClusterDescList,
							ClusterSizeList,
							ClusterProtSeqIDList,
							ClusterSizeByTaxoGroupingList,
							ClusterSizeByProtList)),
				fmt = '%s',
				delimiter = "\t",
				header = "# Cluster ID\tCluster desc\tSequence number\tProtein ID\tNumber of sequences by class\tNumber of sequences by protein")

		return (np.array(ClusterIDList),
			np.array(ClusterDescList),
			np.array(ClusterSizeList),
			np.array(ClusterProtSeqIDList),
			np.array(ClusterSizeByTaxoGroupingList),
			np.array(ClusterSizeByProtList),
			)

	def mkdirs(self):
		'''Return all directories for db storage'''
		print("Generating directories")

		'''Blast dirs'''
		BLASTMainDir		= self.ShelveDir+"/BLAST"
		if os.path.exists(BLASTMainDir):
			'''Clear previous results'''
			_ = subprocess.call("rm -rf %s" %BLASTMainDir, shell = True)
		os.makedirs(BLASTMainDir)

		BLASTQueryFile		= BLASTMainDir+"/Query.fasta"
		BLASTSubjectFile	= BLASTMainDir+"/Subjects.fasta"
		BLASTOutputFile		= BLASTMainDir+"/BLASTOutput.txt"
		BLASTBitScoreFile	= BLASTMainDir+"/BitScoreMat.txt"
		BLASTProtClusterFile= BLASTMainDir+"/ProtClusters.txt"
		ClustersDir			= BLASTMainDir+"/Clusters"; os.makedirs(ClustersDir)
		
		'''HMMER dirs'''
		HMMERDir			= self.ShelveDir+"/HMMER"
		'''Delete existing HMMER libraries'''
		if os.path.exists(HMMERDir):
			'''Clear previous results'''
			_ = subprocess.call("rm -rf %s" %HMMERDir, shell = True)
		os.makedirs(HMMERDir)
		HMMER_PPHMMDir		= HMMERDir+"/HMMER_PPHMMs"; os.makedirs(HMMER_PPHMMDir)
		HMMER_PPHMMDBDir	= HMMERDir+"/HMMER_PPHMMDB"; os.makedirs(HMMER_PPHMMDBDir)
		HMMER_PPHMMDB		= HMMER_PPHMMDBDir+"/HMMER_PPHMMDB"

		VariableShelveDir 	= self.ShelveDir+"/Shelves"

		'''HHsuite dirs, regardless of if it's enabled'''
		HHsuiteDir	= self.ShelveDir+"/HHsuite"
		if os.path.exists(HHsuiteDir):
			'''Clear previous results'''
			_ = subprocess.call("rm -rf %s" %HHsuiteDir, shell = True)
		os.makedirs(HHsuiteDir)
		HHsuite_PPHMMDir = HHsuiteDir + "/HHsuite_PPHMMs";os.makedirs(HHsuite_PPHMMDir)
		HHsuite_PPHMMDBDir= HHsuiteDir +"/HHsuite_PPHMMDB";os.makedirs(HHsuite_PPHMMDBDir)
		HHsuite_PPHMMDB	= HHsuite_PPHMMDBDir+"/HHsuite_PPHMMDB"
			
		return  BLASTQueryFile, BLASTSubjectFile, BLASTOutputFile, BLASTBitScoreFile, \
				BLASTProtClusterFile, ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDBDir, \
				HMMER_PPHMMDB, VariableShelveDir, HHsuiteDir, HHsuite_PPHMMDir, HHsuite_PPHMMDB

	def retrieve_variables(self, fname):
		'''Read pickle file containing VMR data'''
		return pickle.load(open(fname, "rb"))

	def get_genbank(self, genomes):
		'''Check if GenBank files exist; dl if needed. Index & transform.'''
		if not os.path.isfile(self.GenomeSeqFile):
			print("- Download GenBank file\n GenomeSeqFile doesn't exist. GRAViTy is downloading the GenBank file(s).")
			DownloadGenBankFile(self.GenomeSeqFile, genomes["SeqIDLists"])
		
		print("- Reading GenBank file")
		GenBankDict = SeqIO.index(self.GenomeSeqFile, "genbank")
		return {k.split(".")[0]:v for k,v in GenBankDict.items()}

	def sequence_extraction(self, genomes, GenBankDict):
		'''Extract protein sequences from VMR using GenBank data; manually annotate ORFs if needed'''
		print(f"- Extract/predict protein sequences from virus genomes, excluding proteins with lengthes <{self.ProteinLength_Cutoff} aa")
		ProtList, ProtIDList, N_Viruses, Virus_i = [], [], len(genomes["SeqIDLists"]), 1

		for SeqIDList, TranslTable, BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping in zip(genomes["SeqIDLists"], genomes["TranslTableList"], genomes["BaltimoreList"], genomes["OrderList"], genomes["FamilyList"], genomes["SubFamList"], genomes["GenusList"], genomes["VirusNameList"], genomes["TaxoGroupingList"]):
			assert len(SeqIDList) == 1 # RM < removed a for loop here as SeqIDList only ever seemed to be 1 long - was this just for de-nesting?
			SeqID = SeqIDList[0] 
			GenBankRecord	= GenBankDict[SeqID]
			GenBankID	= GenBankRecord.name
			GenBankFeatures	= GenBankRecord.features

			'''Extract protein sequences'''
			ContainProtAnnotation = 0
			for Feature in GenBankFeatures:
				if(Feature.type == 'CDS' and "protein_id" in Feature.qualifiers and "translation" in Feature.qualifiers):
					ContainProtAnnotation = 1
					try:
						ProtName = Feature.qualifiers["product"][0]
					except KeyError:
						try:
							ProtName = Feature.qualifiers["gene"][0]
						except KeyError:
							try:
								ProtName = Feature.qualifiers["note"][0]
							except KeyError:
								ProtName = "Hypothetical protein"
					ProtID = Feature.qualifiers["protein_id"][0]
					ProtSeq = Feature.qualifiers["translation"][0]
					if len(ProtSeq) >= self.ProteinLength_Cutoff:
						ProtRecord = SeqRecord(	Seq(ProtSeq),
									id = GenBankID+"|"+ProtID,
									name = GenBankID+"|"+ProtID,
									description = ProtName,
									annotations = {'taxonomy':[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping]})
						ProtList.append(ProtRecord)
						ProtIDList.append(GenBankID+"|"+ProtID)

			'''If the genome isn't annotated with any ORFs'''
			if ContainProtAnnotation == 0:	
				'''Identifying ORFs'''
				try:
					Starts = self.orf_tranl_table[TranslTable]
				except KeyError:
					'''No ORF found, use standard code'''
					Starts = self.orf_no_match
				
				'''Generate start and stop codon lists'''
				CodonList = [Base1+Base2+Base3 for Base1 in "TCAG" for Base2 in "TCAG" for Base3 in "TCAG"]
				StartCodonList, StopCodonList = [], []
				for i,j in enumerate(Starts):
					if j == "M":
						StartCodonList.append(CodonList[i])
					if j == "*":
						StopCodonList.append(CodonList[i])		

				GenBankSeq, SeqLength, ORF_i = GenBankRecord.seq, len(GenBankSeq), 0
				for _, nuc in [(+1, GenBankSeq), (-1, GenBankSeq.reverse_complement())]:
					'''Split into multople of 3, get in-frame nucleotide seq, split sequence into codons'''
					for frame in range(3):
						length = 3 * ((SeqLength-frame) // 3)					
						nuc_inframe = nuc[frame:(frame+length)]					
						nuc_codonList = [str(nuc_inframe[i:i+3]) for i in range(0, length, 3)]	
						
						'''Find stop codons'''
						StopCodon_indices = [i for i, codon in enumerate(nuc_codonList) if codon in StopCodonList] 
						Coding_Start_IndexList = np.array([-1]+StopCodon_indices)+1
						Coding_End_IndexList = np.array(StopCodon_indices+[len(nuc_codonList)])
						
						ProtSeqList = []
						for i, j in zip(Coding_Start_IndexList, Coding_End_IndexList):
							for k, codon in enumerate(nuc_codonList[i:j]):
								if codon in StartCodonList:
									ProtSeqList.append(Seq("".join(nuc_codonList[i:j][k:])).translate(table = TranslTable))
									break
						
						for ProtSeq in ProtSeqList:
							'''Exclude protein sequences with <'ProteinLength_Cutoff' aa'''
							if len(ProtSeq) >= self.ProteinLength_Cutoff:
								ProtRecord = SeqRecord(	ProtSeq,
											id = GenBankID+"|ORF%s"%ORF_i,
											name = GenBankID+"|ORF%s"%ORF_i,
											description = "Hypothetical protein",
											annotations = {'taxonomy':[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping]})
								ProtList.append(ProtRecord)
								ProtIDList.append(GenBankID+"|ORF%s"%ORF_i)
								ORF_i = ORF_i + 1

			progress_bar(f"\033[K Extract protein sequences: [{'='*int(Virus_i/N_Viruses*20)}] {Virus_i}/{N_Viruses} viruses \r")
			Virus_i += 1

		clean_stdout()

		return ProtList, np.array(ProtIDList)

	def make_blastp_db(self, BLASTSubjectFile, ProtList):
		'''Make BLAST DB'''
		with open(BLASTSubjectFile, "w") as BLASTSubject_txt:
			SeqIO.write(ProtList, BLASTSubject_txt, "fasta")
		
		_ = subprocess.Popen(f"makeblastdb -in {BLASTSubjectFile} -dbtype prot", stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		error_handler(out, err, "makeblastdb")
		
	def blastp_analysis(self, ProtList, BLASTQueryFile, BLASTSubjectFile, BLASTOutputFile, BLASTBitScoreFile):
		'''Perform ALL-VERSUS-ALL BLASTp analysis'''
		print("Perform ALL-VERSUS-ALL BLASTp analysis")
		BitScoreMat, SeenPair, SeenPair_i, N_ProtSeqs = [], {}, 0, len(ProtList)
		BLASTp_outfmt	= '"6 qseqid sseqid pident qcovs qlen slen evalue bitscore"'
		for ProtSeq_i in range(N_ProtSeqs):
			'''BLAST query fasta file'''
			BLASTQuery = ProtList[ProtSeq_i]
			with open(BLASTQueryFile, "w") as BLASTQuery_txt:
				p = SeqIO.write(BLASTQuery, BLASTQuery_txt, "fasta")
			
			'''Perform BLASTp'''
			_ = subprocess.Popen(f'blastp -query {BLASTQueryFile} -db {BLASTSubjectFile} -out {BLASTOutputFile} -evalue {self.BLASTp_evalue_Cutoff} -outfmt {BLASTp_outfmt} -num_alignments {self.BLASTp_num_alignments} -num_threads {self.BLASTp_N_CPUs}', 	
																		stdout = subprocess.PIPE, stderr = subprocess.PIPE,
																		shell = True)
			out, err = _.communicate()
			error_handler(out, err, f"blastp (protein ID = {ProtList[ProtSeq_i].id})")
			
			'''BitScoreMat conditioned on PIden, QCovs, and SCovs'''
			if os.stat(BLASTOutputFile).st_size != 0: 
				'''if BLAST returns something...'''
				with open(BLASTOutputFile, "r") as BLASTOutput_txt:
					for BLASTHit in BLASTOutput_txt.readlines():
						if BLASTHit == "\n": 
							break
						Line			= BLASTHit.split("\t")
						qseqid			= Line[0]
						sseqid			= Line[1]
						pident			= float(Line[2])
						qcovs			= float(Line[3])
						qlen			= float(Line[4])
						slen			= float(Line[5])
						evalue			= float(Line[6])
						bitscore		= float(Line[7][:-1])
						[SeqID_I, SeqID_II]	= sorted([qseqid, sseqid])
						Pair			= ", ".join([SeqID_I, SeqID_II])
						if ((qseqid != sseqid) and (pident >= self.BLASTp_PercentageIden_Cutoff) and (qcovs >= self.BLASTp_QueryCoverage_Cutoff) and ((qcovs*qlen/slen) >= self.BLASTp_SubjectCoverage_Cutoff)):
							if Pair in SeenPair: 
								'''If the pair has already been seen...'''
								if bitscore > BitScoreMat[SeenPair[Pair]][2]: 
									'''...and if the new bitscore is higher'''
									BitScoreMat[SeenPair[Pair]][2] = bitscore
							else:
								SeenPair[Pair] = SeenPair_i
								BitScoreMat.append([SeqID_I, SeqID_II, bitscore])
								SeenPair_i = SeenPair_i+1
			
			progress_bar(f"\033[K BLASTp: [{'='*int(float(ProtSeq_i+1)/N_ProtSeqs*20)}] {ProtSeq_i+1}/{N_ProtSeqs} proteins \r")
		
		clean_stdout()
		
		BitScoreMat = np.array(BitScoreMat)
		'''Save protein-protein similarity scores (BLASTp bit scores)'''
		np.savetxt(fname = BLASTBitScoreFile,
					X	 = BitScoreMat,
					fmt	 = '%s',
					delimiter= "\t",
					header	 = "SeqID_I\tSeqID_II\tBit score"
					)

	def mcl_clustering(self, ProtIDList, BLASTBitScoreFile, BLASTProtClusterFile):
		'''Use Muscle to do clustering on BLASTp bit scores'''
		print("- Doing protein sequence clustering based on BLASTp bit scores, using the MCL algorithm")
		_ = subprocess.Popen(f"mcl {BLASTBitScoreFile} --abc -o {BLASTProtClusterFile} -I {self.ProtClustering_MCLInflation}", stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		err, out = _.communicate()
		error_handler(out, err, f"Something is wrong with mcl")
		
		'''For each cluster found, pull out seen proteins by ID'''
		SeenProtIDList = []
		with open(BLASTProtClusterFile, 'r') as BLASTProtCluster_txt:
			for Cluster in BLASTProtCluster_txt.readlines():
				SeenProtIDList.extend(Cluster.split("\r\n")[0].split("\n")[0].split("\t"))

		'''Append seen protein IDs to clusters'''
		with open(BLASTProtClusterFile, 'a') as BLASTProtCluster_txt:
			BLASTProtCluster_txt.write("\n".join(list(set(ProtIDList)-set(SeenProtIDList))))

	def make_alignments(self, ProtList, ProtIDList, BLASTProtClusterFile, ClustersDir):
		'''Do protein alignments with muscle align, make cluster alignment annotations'''
		print("- Make protein alignments")
		N_Clusters		= LineCount(BLASTProtClusterFile)+1 #Count the number of clusters
		Cluster_i		= 0
		Cluster_MetaDataDict	= {}
		with open(BLASTProtClusterFile, 'r') as BLASTProtCluster_txt:
			for Cluster in BLASTProtCluster_txt.readlines():
				HitList		= []
				TaxoLists	= []
				DescList	= []
				Cluster		= Cluster.split("\n")[0].split("\t")
				for ProtID in Cluster:
					HitList.append(ProtList[np.where(ProtIDList == ProtID)[0][0]])
					TaxoLists.append(HitList[-1].annotations['taxonomy'])
					DescList.append(HitList[-1].description.replace(", "," ").replace(","," ").replace(": ","_").replace(":","_").replace("; "," ").replace(";"," ").replace(" (","/").replace("(","/").replace(")",""))
				
				'''Cluster file'''
				UnAlnClusterFile = ClustersDir+"/Cluster_%s.fasta" %Cluster_i
				with open(UnAlnClusterFile, "w") as UnAlnClusterTXT:
					p = SeqIO.write(HitList, UnAlnClusterTXT, "fasta")
				
				'''Align cluster using muscle'''
				AlnClusterFile = ClustersDir+"/Cluster_%s.fasta" %Cluster_i
				# RM < Need to research how to replace these dead muscle kwargs
				# RM < Harmonize muslce install location with Dockerfile
				# _ = subprocess.Popen("~/repo/muscle/muscle -align %s -output %s"%(UnAlnClusterFile, #-gapopen %s -gapextend %s" %(	UnAlnClusterFile, # RM << LOTS WRONG HERE
				# 										AlnClusterFile),#,
				# 										#MUSCLE_GapOpenCost,
				# 										#MUSCLE_GapExtendCost),
				# 										stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				_ = subprocess.Popen(f"~/repo/muscle/muscle -align {UnAlnClusterFile} -output {AlnClusterFile}",
											stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				err, out = _.communicate()
				error_handler(out, err, f"Something is wrong with muscle (Cluster_{Cluster_i}):")

				'''Cluster annotations'''
				Cluster_MetaDataDict[Cluster_i] = {	"Cluster":Cluster,
									"DescList":DescList,
									"TaxoLists":TaxoLists,
									"AlignmentLength":AlignIO.read(AlnClusterFile, "fasta").get_alignment_length()
									}
				Cluster_i += 1
				progress_bar(f"\033[K Make protein alignments: [{'='*int(float(Cluster_i)/N_Clusters*20)}] {Cluster_i}/{N_Clusters} alignments \r")
		
		clean_stdout()
		return Cluster_MetaDataDict

	def main(self):
		'''RM < DOCSTRING'''
		section_header("#Build a database of virus protein profile hidden Markov models (PPHMMs) #")

		'''Build db dirs'''
		BLASTQueryFile, BLASTSubjectFile, BLASTOutputFile, BLASTBitScoreFile, \
				BLASTProtClusterFile, ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDBDir, \
				HMMER_PPHMMDB, VariableShelveDir, HHsuiteDir, HHsuite_PPHMMDir, HHsuite_PPHMMDB = self.mkdirs()		
				
		'''Retrieve Variables'''
		if self.IncludeIncompleteGenomes == True:
			VariableShelveFile = VariableShelveDir + "/ReadGenomeDescTable.AllGenomes.p"
		elif self.IncludeIncompleteGenomes == False:
			VariableShelveFile = VariableShelveDir + "/ReadGenomeDescTable.CompleteGenomes.p"

		genomes =  self.retrieve_variables(VariableShelveFile)

		'''Get GenBank files if not exists'''
		GenBankDict = self.get_genbank(genomes)
			
		'''Sequence extraction'''
		ProtList, ProtIDList = self.sequence_extraction(genomes, GenBankDict)

		'''Make BLAST DB'''
		self.make_blastp_db(BLASTSubjectFile, ProtList)

		'''Do BLASTp analysis, save output'''
		self.blastp_analysis(ProtList, BLASTQueryFile, BLASTSubjectFile, BLASTOutputFile, BLASTBitScoreFile)
		
		'''Cluster using Muscle'''
		self.mcl_clustering(ProtIDList, BLASTBitScoreFile, BLASTProtClusterFile)
		
		'''Make Alignments'''
		Cluster_MetaDataDict = self.make_alignments(ProtList, ProtIDList, BLASTProtClusterFile, ClustersDir)
		breakpoint()

		###########################################################
		'''Merge & Make protein alignmens'''
		# RM < Make function
		if self.N_AlignmentMerging != 0:
			if self.N_AlignmentMerging > 0:
				print(f"- Merge protein alignments, {self.N_AlignmentMerging} rounds of merging")
			elif self.N_AlignmentMerging < 0:
				print("- Merge protein alignments until exhausted")
			print("\tMake HHsuite PPHMMs from protein alignments")
			for Cluster_i in range(len(Cluster_MetaDataDict)):
				AlnClusterFile		= ClustersDir+"/Cluster_%s.fasta" %Cluster_i
				HHsuite_PPHMMFile	= HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %Cluster_i
				_ = subprocess.Popen("hhmake -i %s -o %s -seq %s -name Cluster_%s -id 100 -M 50 -v 0" %(	AlnClusterFile,
																HHsuite_PPHMMFile,
																len(Cluster_MetaDataDict[Cluster_i]["Cluster"])+1,
																Cluster_i),
																stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				error_handler(out, err, f"Something is wrong with turning Cluster_{Cluster_i} into a PPHMM by hhmake.")
				progress_bar(f"\033[K Make HHsuite PPHMMs: [{'='*int(float(Cluster_i+1)/len(Cluster_MetaDataDict)*20)}] {Cluster_i+1}/{len(Cluster_MetaDataDict)} PPHMMs \r")

			clean_stdout()
			
			print("\tMake a HHsuite PPHMM DB")
			_ = subprocess.Popen("ffindex_build -s %s_hhm.ffdata %s_hhm.ffindex %s" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB, HHsuite_PPHMMDir), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			
			print("\tMerge protein alignments")
			AlignmentMerging_i_round = 0
			while True:
				if AlignmentMerging_i_round >= self.N_AlignmentMerging and self.N_AlignmentMerging >= 0:
					print("Alignment merging complete")
					break
				
				if self.HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging == True:
					print("\t\tHMMER_PPHMMDB_ForEachRoundOfPPHMMMerging == True. Make a HMMER PPHMM DB. (Round %s)" %AlignmentMerging_i_round)
					_ = self.Make_HMMER_PPHMM_DB(	HMMER_PPHMMDir = HMMER_PPHMMDir,
									HMMER_PPHMMDB = HMMER_PPHMMDBDir+"/HMMER_PPHMMDB_%s" %AlignmentMerging_i_round,
									ClustersDir = ClustersDir,
									Cluster_MetaDataDict = Cluster_MetaDataDict)
					
					_ = subprocess.Popen("find %s -type f -name '*.hmm' -delete" %HMMER_PPHMMDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
					out, err = _.communicate()
				
				print("\t\tRound %s"%(AlignmentMerging_i_round + 1))
				print("\t\t\tDetermine PPHMM-PPHMM similarity scores (ALL-VERSUS-ALL hhsearch)")

				hhsearchDir		= HHsuiteDir+"/hhsearch_"+"".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(hhsearchDir)
				hhsearchOutFile		= hhsearchDir+"/hhsearch.stdout.hhr"
				N_PPHMMs		= LineCount("%s_hhm.ffindex"%HHsuite_PPHMMDB)
				
				SeenPair		= {}
				SeenPair_i		= 0
				PPHMMSimScoreCondensedMat= []
				for PPHMM_i in range(0, N_PPHMMs):
					HHsuite_PPHMMFile = HHsuite_PPHMMDir + "/PPHMM_%s.hhm" %PPHMM_i
					_ = subprocess.Popen("hhsearch -i %s -d %s -o %s -e %s -E %s -z 1 -b 1 -id 100 -global -v 0 -cpu %s" %(	HHsuite_PPHMMFile,
																		HHsuite_PPHMMDB+"_hhm.ffdata",
																		hhsearchOutFile,
																		self.HHsuite_evalue_Cutoff,
																		self.HHsuite_evalue_Cutoff,
																		self.HHsuite_N_CPUs,
																		),
																		stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
					out, err = _.communicate()

					# RM < redo error handler
					# if err != "":
					# 	print("Something is wrong with hhsearching PPHMM %s againt the PPHMM database" %PPHMM_i)
					# 	print("#"*50+"out"+"#"*50)
					# 	print(out)
					# 	print("#"*50+"err"+"#"*50)
					# 	print(err)
					# 	print("_"*100)
					# 	while True:
					# 		Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
					# 		if Input == "N" or Input == "n":
					# 			raise SystemExit ("GRAViTy terminated.")
					# 		elif Input == "Y" or Input == "y":
					# 			print("Continue GRAViTy.")
					# 			break
					# 		else:
					# 			print("Input can only be 'Y' or 'N'.")
					
					with open(hhsearchOutFile, 'r') as hhsearchOut_txt:
						Content		= hhsearchOut_txt.readlines()
						QueryLength	= int(Content[1].split()[1])
						for Line in Content[9:]:
							if Line == "\n":
								break
							else:
								Line  		= Line.replace("("," ").replace(")"," ").split()
								PPHMM_j		= int(Line[1].split("_")[1])
								evalue		= float(Line[3])
								pvalue		= float(Line[4])
								PPHMMSimScore	= float(Line[5])
								Col		= float(Line[7])
								SubjectLength	= int(Line[10])
								qcovs		= Col/QueryLength*100
								scovs		= Col/SubjectLength*100
								if (evalue <= self.HHsuite_evalue_Cutoff and pvalue <= self.HHsuite_pvalue_Cutoff and qcovs >= self.HHsuite_QueryCoverage_Cutoff and scovs >= HHsuite_SubjectCoverage_Cutoff):
									Pair	= ", ".join(sorted(map(str,[PPHMM_i, PPHMM_j])))
									if Pair in SeenPair: 
										'''If the pair has already been seen...'''
										if PPHMMSimScore > PPHMMSimScoreCondensedMat[SeenPair[Pair]][2]: 
											'''...and if the new PPHMMSimScore is higher...'''
											PPHMMSimScoreCondensedMat[SeenPair[Pair]][2] = PPHMMSimScore
									else:
										SeenPair[Pair] = SeenPair_i
										PPHMMSimScoreCondensedMat.append([PPHMM_i, PPHMM_j, PPHMMSimScore])
										SeenPair_i = SeenPair_i+1
					
					_ = subprocess.Popen("rm %s" %hhsearchOutFile, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
					out, err = _.communicate()
					
					'''Progress bar'''
					sys.stdout.write("\033[K" + "hhsearch: [%-20s] %d/%d PPHMMs" % ('='*int(float(PPHMM_i+1)/N_PPHMMs*20), PPHMM_i+1, N_PPHMMs) + "\r")
					sys.stdout.flush()
				
				sys.stdout.write("\033[K")
				sys.stdout.flush()
				
				PPHMMSimScoreCondensedMat	= np.array(PPHMMSimScoreCondensedMat)
				PPHMMSimScoreMat		= np.zeros((N_PPHMMs, N_PPHMMs))
				PPHMMSimScoreMat[list(map(int, PPHMMSimScoreCondensedMat[:,0])), list(map(int, PPHMMSimScoreCondensedMat[:,1]))] = list(map(float, PPHMMSimScoreCondensedMat[:,2]))
				PPHMMSimScoreMat[list(map(int, PPHMMSimScoreCondensedMat[:,1])), list(map(int, PPHMMSimScoreCondensedMat[:,0]))] = list(map(float, PPHMMSimScoreCondensedMat[:,2]))
				PPHMMSimScoreCondensedMat	= np.array([PPHMMSimScorePair for PPHMMSimScorePair in PPHMMSimScoreCondensedMat if PPHMMSimScorePair[0] < PPHMMSimScorePair[1]])
				
				PPHMMSimScoreCondensedMatFile	= hhsearchDir+"/PPHMMSimScoreCondensedMat.txt"
				np.savetxt(	fname	= PPHMMSimScoreCondensedMatFile,
						X	= PPHMMSimScoreCondensedMat,
						fmt	= '%s',
						delimiter= "\t",
						header	= "PPHMM_i\tPPHMM_j\tPPHMMSimScore")
				
				print("\t\t\tCluster PPHMMs based on hhsearch scores, using the MCL algorithm")
				PPHMMClustersFile	= hhsearchDir+"/PPHMMClusters.txt"
				_ = subprocess.Popen("mcl %s --abc -o %s -I %s" %(PPHMMSimScoreCondensedMatFile, PPHMMClustersFile, self.PPHMMClustering_MCLInflation), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				_, out = _.communicate()
				
				SeenProtIDList = []
				with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
					for Cluster in PPHMMClusters_txt.readlines():
						SeenProtIDList.extend(Cluster.split("\n")[0].split("\t"))
				
				with open(PPHMMClustersFile, 'a') as PPHMMClusters_txt:
					PPHMMClusters_txt.write("\n".join(list(set(map(str,list(map(float,list(range(0, N_PPHMMs))))))-set(SeenProtIDList))))
				
				print("\t\t\tCheck if there are alignments to be merged")
				with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
					N_PPHMMs_AfterMerging = len(PPHMMClusters_txt.readlines())
				
				if N_PPHMMs_AfterMerging == N_PPHMMs:
					print("\t\t\t\tNo alignments to be merged. Stop alignment merging process")
					_ = subprocess.Popen("rm -rf %s" %hhsearchDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
					out, err = _.communicate()
					break
				else:
					print("\t\t\t\tMerge %d alignments to make %d alignments" %(N_PPHMMs, N_PPHMMs_AfterMerging))
				
				print("\t\t\tMerge protein alignments and remake HHsuite PPHMMs")
				SelfSimScoreList				= PPHMMSimScoreMat.diagonal()
				PPHMMDissimScoreMat				= 1 - np.transpose(PPHMMSimScoreMat**2/SelfSimScoreList)/SelfSimScoreList
				PPHMMDissimScoreMat[PPHMMDissimScoreMat<0]	= 0
				AfterMergingPPHMM_IndexList			= []
				AfterMergingPPHMM_i				= 1.0
				with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
					for PPHMMCluster in PPHMMClusters_txt.readlines():
						PPHMMCluster = list(map(int, list(map(float, PPHMMCluster.split("\n")[0].split("\t")))))
						AfterMergingPPHMM_IndexList.append(min(PPHMMCluster))
						if len(PPHMMCluster) >= 2:
							PPHMMDissimScoreMat_Subset = PPHMMDissimScoreMat[PPHMMCluster][:,PPHMMCluster]
							PPHMMTreeNewick = DistMat2Tree (DistMat	= PPHMMDissimScoreMat_Subset,
											LeafList= PPHMMCluster,
											Dendrogram_LinkageMethod	= "average")
							PPHMMTreeNewick = Tree(PPHMMTreeNewick)
							_ 		= PPHMMTreeNewick.ladderize()
							PPHMMTreeNewick	= PPHMMTreeNewick.write(format = 9)
							
							while True:
								m = re.search(r"\((\d+),(\d+)\)", PPHMMTreeNewick)
								if not m:
									_ = subprocess.Popen("muscle -in %s -out %s -refine" %(	ClusterFile_i,
																ClusterFile_i),
																stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
									err, out = _.communicate()
									break
								PPHMM_i, PPHMM_j = sorted([int(m.group(1)),int(m.group(2))])
								PPHMMTreeNewick = re.sub(r"\((\d+),(\d+)\)", str(PPHMM_i), PPHMMTreeNewick, count=1)
								
								ClusterFile_i = ClustersDir+"/Cluster_%s.fasta" %PPHMM_i
								ClusterFile_j = ClustersDir+"/Cluster_%s.fasta" %PPHMM_j
								HHsuite_PPHMMFile_j = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_j
								_ = subprocess.Popen("muscle -profile -in1 %s -in2 %s -out %s -gapopen %s -gapextend %s" %(	ClusterFile_i,
																				ClusterFile_j,
																				ClusterFile_i,
																				self.MUSCLE_GapOpenCost,
																				self.MUSCLE_GapExtendCost),
																				stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
								err, out = _.communicate()
								_ = subprocess.Popen("rm %s %s" %(ClusterFile_j, HHsuite_PPHMMFile_j), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
								out, err = _.communicate()
								
								Cluster_MetaDataDict[PPHMM_i]["Cluster"]	= Cluster_MetaDataDict[PPHMM_i]["Cluster"] + Cluster_MetaDataDict[PPHMM_j]["Cluster"]
								Cluster_MetaDataDict[PPHMM_i]["DescList"]	= Cluster_MetaDataDict[PPHMM_i]["DescList"] + Cluster_MetaDataDict[PPHMM_j]["DescList"]
								Cluster_MetaDataDict[PPHMM_i]["TaxoLists"]	= Cluster_MetaDataDict[PPHMM_i]["TaxoLists"] + Cluster_MetaDataDict[PPHMM_j]["TaxoLists"]
								del Cluster_MetaDataDict[PPHMM_j]
							
							HHsuite_PPHMMFile_i = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_i
							Cluster_MetaDataDict[PPHMM_i]["AlignmentLength"]	= AlignIO.read(ClusterFile_i, "fasta").get_alignment_length()
							_ = subprocess.Popen("hhmake -i %s -o %s -v 0 -seq %s -name Cluster_%s -id 100 -M 50" %(	ClusterFile_i,
																			HHsuite_PPHMMFile_i,
																			len(Cluster_MetaDataDict[PPHMM_i]["Cluster"])+1,
																			PPHMM_i),
																			stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
							out, err = _.communicate()

							# RM < Redo error handler
							# if err != "":
							# 	print("Something is wrong with constructing a PPHMM from cluster %s" %PPHMM_i)
							# 	print("#"*50+"out"+"#"*50)
							# 	print(out)
							# 	print("#"*50+"err"+"#"*50)
							# 	print(err)
							# 	print("_"*100)
							# 	while True:
							# 		Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
							# 		if Input == "N" or Input == "n":
							# 			raise SystemExit("GRAViTy terminated.")
							# 		elif Input == "Y" or Input == "y":
							# 			print("Continue GRAViTy.")
							# 			break
							# 		else:
							# 			print("Input can only be 'Y' or 'N'.")
							
						elif len(PPHMMCluster) == 1:
							pass
						
						'''Progress bar'''
						sys.stdout.write("\033[K" + "Merge alignments and make new PPHMMs: [%-20s] %d/%d PPHHMs" % ('='*int(AfterMergingPPHMM_i/N_PPHMMs_AfterMerging*20), AfterMergingPPHMM_i, N_PPHMMs_AfterMerging) + "\r")
						sys.stdout.flush()
						AfterMergingPPHMM_i = AfterMergingPPHMM_i + 1
				
				sys.stdout.write("\033[K")
				sys.stdout.flush()
				
				print("\t\t\tRename protein alignments and their associated PPHMMs")
				AfterMergingPPHMM_IndexList	= sorted(AfterMergingPPHMM_IndexList)
				AfterMergingPPHMM_i		= 0
				for PPHMM_i in AfterMergingPPHMM_IndexList:
					Cluster_MetaDataDict[AfterMergingPPHMM_i] = Cluster_MetaDataDict.pop(PPHMM_i)
					
					ClusterFile_i = ClustersDir+"/Cluster_%s.fasta" %PPHMM_i
					ClusterFile_j = ClustersDir+"/Cluster_%s.fasta" %AfterMergingPPHMM_i
					_ = subprocess.Popen("mv %s %s" %(ClusterFile_i, ClusterFile_j), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
					out, err = _.communicate()
					
					HHsuite_PPHMMFile_i = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_i
					HHsuite_PPHMMFile_j = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %AfterMergingPPHMM_i
					_ = subprocess.Popen("mv %s %s" %(HHsuite_PPHMMFile_i, HHsuite_PPHMMFile_j), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
					out, err = _.communicate()
					
					with open(HHsuite_PPHMMFile_j, "r+") as HHsuite_PPHMM_txt:
						contents = HHsuite_PPHMM_txt.readlines()
						contents[1] = "NAME  Cluster_%s\n" %AfterMergingPPHMM_i
						contents = "".join(contents)
						HHsuite_PPHMM_txt.seek(0)		#Put cursor at the beginning of the file
						HHsuite_PPHMM_txt.write(contents)	#Write the contents
						HHsuite_PPHMM_txt.truncate()		#Delete everything after the cursor
					
					AfterMergingPPHMM_i = AfterMergingPPHMM_i + 1
				
				print("\t\t\tRebuild the HHsuite PPHMM database\n")
				_ = subprocess.Popen("rm %s_hhm.ffdata %s_hhm.ffindex" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				
				if list(set([f.split(".")[-1] for f in os.listdir(HHsuite_PPHMMDir)]))!=["hhm"]:
					print("There are some other files/folders other than HHsuite PPHMMs in the folder %s. Remove them first." %HHsuite_PPHMMDir)
				
				_ = subprocess.Popen("ffindex_build -s %s_hhm.ffdata %s_hhm.ffindex %s" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB, HHsuite_PPHMMDir), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				
				_ = subprocess.Popen("rm -rf %s" %hhsearchDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				
				AlignmentMerging_i_round = AlignmentMerging_i_round + 1
			
			print("\tAlignment merging is done. Delete HHsuite shelve directory")
			_ = subprocess.Popen("rm -rf %s" %HHsuiteDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
		
		'''Make HMMER DB'''
		# RM < make function
		print("- Make HMMER PPHMMDB and its summary file")
		(ClusterIDList,
		ClusterDescList,
		ClusterSizeList,
		ClusterProtSeqIDList,
		ClusterSizeByTaxoGroupingList,
		ClusterSizeByProtList) = self.Make_HMMER_PPHMM_DB(	HMMER_PPHMMDir = HMMER_PPHMMDir,
								HMMER_PPHMMDB = HMMER_PPHMMDB,
								ClustersDir = ClustersDir,
								Cluster_MetaDataDict = Cluster_MetaDataDict)
		VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
		# RM < Pickle
		Parameters = shelve.open(VariableShelveFile,"n")
		for key in ["ClusterIDList",
				"ClusterDescList",
				"ClusterSizeList",
				"ClusterProtSeqIDList",
				"ClusterSizeByTaxoGroupingList",
				"ClusterSizeByProtList",
				]:
			try:
				Parameters[key] = locals()[key]
				print("\t"+key)
			except TypeError:
				pass
		
		Parameters.close()

