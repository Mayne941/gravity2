from app.src.read_genome_desc_table import ReadGenomeDescTable
from app.src.pphmmdb_construction import PPHMMDBConstruction
from app.src.ucf_virus_annotator import UcfVirusAnnotator
from app.src.virus_classification import VirusClassificationAndEvaluation

from app.utils.str_to_bool import str2bool
from app.utils.benchmark import benchmark_start, benchmark_end
from app.utils.generate_logs import Log_Generator_Pl2

import optparse, os, time

class Pipeline_II:
	def __init__(self, payload):
		self.options = payload
		'''Create directories'''
		if not os.path.exists(self.options["ShelveDir_UcfVirus"]):
			os.makedirs(self.options["ShelveDir_UcfVirus"])
		'''Logs'''
		self.log_gen = Log_Generator_Pl2(self.options, self.options["ShelveDir_UcfVirus"])

	def main(self):
		actual_start = time.time()
		# RM < Update for API
		logs = self.log_gen.entrypoint()

		'''Catch bad database flags'''
		if (self.options['Database'] != None and self.options['Database_Header'] == None):
			raise optparse.OptionValueError(f"You have specified DATABASE as {self.options['Database']}, 'Database_Header' cannot be 'None'")
		if (self.options['Database'] == None and self.options['Database_Header'] != None):
			Proceed = input (f"You have specified 'Database_Header' as {self.options['Database_Header']}, but 'Database' is 'None'. GRAViTy will analyse all genomes. Do you want to proceed? [Y/n]: ")
			if Proceed != "Y":
				raise SystemExit("GRAViTy terminated.")
		
		'''I: Fire Read Genom Desc Table'''
		[print(log_text) for log_text in logs[0]]
		start = benchmark_start("ReadGenomeDescTable")
		ReadGenomeDescTable(
			GenomeDescTableFile	= self.options['GenomeDescTableFile_UcfVirus'],
			ShelveDir			= self.options['ShelveDir_UcfVirus'],
			Database			= self.options['Database'],
			Database_Header		= self.options['Database_Header'],
			)
		benchmark_end("ReadGenomeDescTable", start)
		
		'''II: Fire PPHMDB Constructor'''
		if str2bool(self.options['UseUcfVirusPPHMMs']) == True:
			[print(log_text) for log_text in logs[1]]				
			start = benchmark_start("PPHMMDBConstruction")
			PPHMMDBConstruction (
				GenomeSeqFile = self.options['GenomeSeqFile_UcfVirus'],
				ShelveDir = self.options['ShelveDir_UcfVirus'],
				ProteinLength_Cutoff = self.options['ProteinLength_Cutoff'],
				IncludeIncompleteGenomes = str2bool(self.options['IncludeProteinsFromIncompleteGenomes']),
				BLASTp_evalue_Cutoff = self.options['BLASTp_evalue_Cutoff'],
				BLASTp_PercentageIden_Cutoff = self.options['BLASTp_PercentageIden_Cutoff'],
				BLASTp_QueryCoverage_Cutoff = self.options['BLASTp_QueryCoverage_Cutoff'],
				BLASTp_SubjectCoverage_Cutoff = self.options['BLASTp_SubjectCoverage_Cutoff'],
				BLASTp_num_alignments = self.options['BLASTp_num_alignments'],
				BLASTp_N_CPUs = self.options['BLASTp_N_CPUs'],
				MUSCLE_GapOpenCost = self.options['MUSCLE_GapOpenCost'],
				MUSCLE_GapExtendCost = self.options['MUSCLE_GapExtendCost'],
				ProtClustering_MCLInflation = self.options['ProtClustering_MCLInflation'],
				N_AlignmentMerging = self.options['N_AlignmentMerging'],
				HHsuite_evalue_Cutoff = self.options['HHsuite_evalue_Cutoff'],
				HHsuite_pvalue_Cutoff = self.options['HHsuite_pvalue_Cutoff'],
				HHsuite_N_CPUs = self.options['HHsuite_N_CPUs'],
				HHsuite_QueryCoverage_Cutoff = self.options['HHsuite_QueryCoverage_Cutoff'],
				HHsuite_SubjectCoverage_Cutoff = self.options['HHsuite_SubjectCoverage_Cutoff'],
				PPHMMClustering_MCLInflation = self.options['PPHMMClustering_MCLInflation_ForAlnMerging'],
				HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging = str2bool(self.options['HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging']),
				)
			benchmark_end("PPHMMDBConstruction", start)

		'''III: Fire UCF Virus Annotator'''
		[print(log_text) for log_text in logs[2]]
		start = benchmark_start("UcfVirusAnnotator")
		UcfVirusAnnotator (
			GenomeSeqFile_UcfVirus = self.options['GenomeSeqFile_UcfVirus'],
			ShelveDir_UcfVirus = self.options['ShelveDir_UcfVirus'],
			ShelveDirs_RefVirus = self.options['ShelveDirs_RefVirus'],
			IncludeIncompleteGenomes_UcfVirus = str2bool(self.options['AnnotateIncompleteGenomes_UcfVirus']),
			IncludeIncompleteGenomes_RefVirus = str2bool(self.options['UsingDatabaseIncludingIncompleteRefViruses']),
			SeqLength_Cutoff = 0,
			HMMER_N_CPUs = self.options['HMMER_N_CPUs'],
			HMMER_C_EValue_Cutoff = self.options['HMMER_C_EValue_Cutoff'],
			HMMER_HitScore_Cutoff = self.options['HMMER_HitScore_Cutoff'],
			)
		benchmark_end("UcfVirusAnnotator", start)

		'''IV: Fire Classifiers'''
		[print(log_text) for log_text in logs[3]]
		start = benchmark_start("VirusClassificationAndEvaluation")
		VirusClassificationAndEvaluation (
			ShelveDir_UcfVirus = self.options['ShelveDir_UcfVirus'],
			ShelveDirs_RefVirus = self.options['ShelveDirs_RefVirus'],
			IncludeIncompleteGenomes_UcfVirus = str2bool(self.options['AnnotateIncompleteGenomes_UcfVirus']),
			IncludeIncompleteGenomes_RefVirus = str2bool(self.options['UsingDatabaseIncludingIncompleteRefViruses']),
			UseUcfVirusPPHMMs = str2bool(self.options['UseUcfVirusPPHMMs']),
			GenomeSeqFile_UcfVirus = self.options['GenomeSeqFile_UcfVirus'],
			GenomeSeqFiles_RefVirus = self.options['GenomeSeqFiles_RefVirus'],
			SeqLength_Cutoff = 0,
			HMMER_N_CPUs = int(self.options['HMMER_N_CPUs']),
			HMMER_C_EValue_Cutoff = float(self.options['HMMER_C_EValue_Cutoff']),
			HMMER_HitScore_Cutoff = float(self.options['HMMER_HitScore_Cutoff']),
			SimilarityMeasurementScheme = self.options['SimilarityMeasurementScheme'],
			p = float(self.options['p']),
			Dendrogram_LinkageMethod = self.options['Dendrogram_LinkageMethod'],
			DatabaseAssignmentSimilarityScore_Cutoff = float(self.options['DatabaseAssignmentSimilarityScore_Cutoff']),
			N_PairwiseSimilarityScores = int(self.options['N_PairwiseSimilarityScores']),
			Heatmap_WithDendrogram = str2bool(self.options['Heatmap_WithDendrogram']),
			Heatmap_DendrogramSupport_Cutoff = float(self.options['Heatmap_DendrogramSupport_Cutoff']),
			Bootstrap = str2bool(self.options['Bootstrap']),
			N_Bootstrap = int(self.options['N_Bootstrap']),
			Bootstrap_method = self.options['Bootstrap_method'],
			Bootstrap_N_CPUs = int(self.options['Bootstrap_N_CPUs']),
			VirusGrouping = str2bool(self.options['VirusGrouping']),
			)
		benchmark_end("VirusClassificationAndEvaluation", start)	

		print(f"Time to complete: {time.time() - actual_start}")

if __name__ == '__main__':
	pl = Pipeline_II()
	pl.main()



