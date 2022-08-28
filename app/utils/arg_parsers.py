from app.utils.check_input import check_FILEPATH, check_FILEPATHS, check_PERCENT, check_PROB, check_POS, check_POSINTEGER, check_NONNEG, check_NONNEGINTEGER, check_NONPOS, check_N_AlignmentMerging

import optparse
import textwrap
import multiprocessing
'''
DEPRECATED FROM APP VERSION 1.0 (PYTHON3)
'''

class IndentedHelpFormatterWithNL(optparse.IndentedHelpFormatter):
	def format_description(self, description):
		if not description: return ""
		desc_width = self.width - self.current_indent
		indent = " "*self.current_indent
		bits = description.split('\n')
		formatted_bits = [
			textwrap.fill(bit,
				desc_width,
				initial_indent=indent,
				subsequent_indent=indent)
			for bit in bits]
		result = "\n".join(formatted_bits) + "\n"
		return result
	
	def format_option(self, option):
		result = []
		opts   = self.option_strings[option]
		opt_width = self.help_position - self.current_indent - 2
		if len(opts) > opt_width:
			opts = "%*s%s\n" % (self.current_indent, "", opts)
			indent_first = self.help_position
		else: # start help on same line as opts
			opts = "%*s%-*s	" % (self.current_indent, "", opt_width, opts)
			indent_first = 0
		result.append(opts)
		if option.help:
			help_text = self.expand_default(option)
			# Everything is the same up through here
			help_lines = []
			for para in help_text.split("\n"):
				help_lines.extend(textwrap.wrap(para, self.help_width))
			# Everything is the same after here
			result.append("%*s%s\n" % (
				indent_first, "", help_lines[0]))
			result.extend(["%*s%s\n" % (self.help_position, "", line)
				for line in help_lines[1:]])
		elif opts[-1] != "\n":
			result.append("\n")
		return "".join(result)

def generate_pipeline_i_arguments():
	parser = optparse.OptionParser(	usage = "usage: %prog [options]",
					version = "%prog 1.1.0",
					formatter = IndentedHelpFormatterWithNL())
	
	parser.add_option(	'--GenomeDescTableFile',
				dest 	= "GenomeDescTableFile",
				help 	= "Full path to the Virus Metadata Resource (VMR) tab delimited file, wth headers. "
					  "VMR can be downloaded from https://talk.ictvonline.org/taxonomy/vmr/",
				metavar	= "FILEPATH",
				type	= "string",
				action	= "callback",
				callback= check_FILEPATH,
	)
	parser.add_option(	'--ShelveDir',
				dest	= "ShelveDir",
				help	= "Full path to the shelve directory, storing GRAViTy outputs.",
				metavar	= "DIRECTORYPATH",
				type	= "string",
	)
	parser.add_option(	'--Database',
				dest	= "Database",
				default	= None,
				help	= "GRAViTy will only analyse genomes that are labelled with DATABASE in the database column. "
					  "The database column can be specified by the DATABASE HEADER argument. "
					  "If 'None', all entries are analysed. [default: %default]",
				metavar	= "DATABASE",
				type	= "string",
	)
	parser.add_option(	'--Database_Header',
				dest	= "Database_Header",
				default	= None,
				help	= "The header of the database column. Cannot be 'None' if DATABASE is specified. [default: %default]",
				metavar	= "DATABASE COLUMN HEADER",
				type	= "string",
	)
	parser.add_option(	'--TaxoGrouping_Header',
				dest	= "TaxoGrouping_Header",
				default	= "Family",
				help	= "The header of the Taxonomic grouping column. [default: %default]",
				metavar	= "TAXOGROUPING COLUMN HEADER",
				type	= "string",
	)
	parser.add_option(	'--TaxoGroupingFile',
				dest 	= "TaxoGroupingFile",
				default	= None,
				help 	= "It is possible that the user might want to associate different viruses with different taxonomic assignment levels, "
					  "e.g. family assignments for some viruses, and subfamily or genus assignments for some other viruses, etc. "
					  "To accomodate this need, the user can either add a column in the VMR file, and use --TaxoGrouping_Header to specify the column (see --TaxoGrouping_Header). "
					  "Alternatively, the user can provide a file (with no header) that contains a single column of taxonomic groupings for all viruses in the order that appears in the VMR file. "
					  "The user can specify the full path to the taxonomic grouping file using this options. "
					  "If this option is used, it will override the one specified by --TaxoGrouping_Header. [default: %default]",
				metavar	= "FILEPATH",
				type	= "string",
				action	= "callback",
				callback= check_FILEPATH,
	)
	parser.add_option(	'--GenomeSeqFile',
				dest	= "GenomeSeqFile",
				help	= "Full path to the genome sequence GenBank file. "
					  "If the file doesn't exist, GRAViTy will download the sequences from the NCBI database using accession numbers specified in the VMR file, 'Virus GENBANK accession' column",
				metavar	= "FILEPATH",
				type	= "string",
	)
	
	ProtExtractionGroup 	= optparse.OptionGroup(	parser		= parser,
							title		= "Protein extraction options",
							description	= "Protein sequences are extracted from the genome sequence GenBank File. "
									  "If the file contains protein sequences, GRAViTy will use those annotations. "
									  "Otherwise, protein sequences will be inferred using the genetic code specified in the VMR file, 'Genetic code table' column")
	ProtExtractionGroup.add_option(	'--ProteinLength_Cutoff',
					dest	= "ProteinLength_Cutoff",
					default	= 100,
					help	= "Proteins with length < LENGTH aa will be ignored [default: %default]",
					metavar	= "LENGTH",
					type	= "int",
					action	= "callback",
					callback= check_NONNEGINTEGER,
	)
	ProtExtractionGroup.add_option(	'--IncludeProteinsFromIncompleteGenomes',
					dest	= "IncludeProteinsFromIncompleteGenomes",
					default	= True,
					help	= "Include protein sequences from incomplete genomes to the database if True. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	parser.add_option_group(ProtExtractionGroup)
	
	ProtClusteringGroup	= optparse.OptionGroup( parser		= parser,
							title		= "Protein clustering options",
							description	= "Protein sequences are clustered based on ALL-VERSUS-ALL BLASTp hit scores. "
									  "In summary, a protein database is constructed, and all protein sequences (queries) are searched against the database one by one using BLASTp. "
									  "Protein-protein pairwise similarity scores (bit scores) are collected, and subsequently used for protein clustering by the MCL algorithm. "
									  "Sequneces in each cluster are then aligned by using MUSCLE, and turned into a protein profile hidden Markov model (PPHMM).")
	ProtClusteringGroup.add_option(	'--BLASTp_evalue_Cutoff',
					dest	= "BLASTp_evalue_Cutoff",
					default	= 1E-3,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with an E-value > E-VALUE will be ignored. [default: %default]",
					metavar	= "E-VALUE",
					type	= "float",
					action	= "callback",
					callback= check_POS,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_PercentageIden_Cutoff',
					dest	= "BLASTp_PercentageIden_Cutoff",
					default	= 50,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with a percentage identity < PERCENTAGE IDENTITY will be ignored [default: %default]",
					metavar	= "PERCENTAGE IDENTITY",
					type	= "float",
					action	= "callback",
					callback= check_PERCENT,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_QueryCoverage_Cutoff',
					dest	= "BLASTp_QueryCoverage_Cutoff",
					default	= 75,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with a query coverage < COVERAGE will be ignored [default: %default]",
					metavar	= "COVERAGE",
					type	= "float",
					action	= "callback",
					callback= check_PERCENT,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_SubjectCoverage_Cutoff',
					dest	= "BLASTp_SubjectCoverage_Cutoff",
					default	= 75,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with a subject coverage < COVERAGE will be ignored [default: %default]",
					metavar	= "COVERAGE",
					type	= "float",
					action	= "callback",
					callback= check_PERCENT,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_num_alignments',
					dest	= "BLASTp_num_alignments",
					default	= 1000000,
					help	= "Maximum number of sequences to be considered in a BLASTp search. [default: %default]",
					metavar	= "BLAST ALIGNMENT NUMBER",
					type	= "int",
					action	= "callback",
					callback= check_POSINTEGER,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_N_CPUs',
					dest	= "BLASTp_N_CPUs",
					default	= multiprocessing.cpu_count(),
					help	= "The number of threads (CPUs) to use in the BLASTp search. [default: %default - all threads]",
					metavar	= "THREADS",
					type	= "int",
					action	= "callback",
					callback= check_POSINTEGER,
	)
	
	ProtClusteringGroup.add_option(	'--MUSCLE_GapOpenCost',
					dest	= "MUSCLE_GapOpenCost",
					default	= -3.0,
					help	= "MUSCLE gap opening panelty for aligning protein sequences. [default: %default]",
					metavar	= "gap opening panelty",
					type	= "float",
					action	= "callback",
					callback= check_NONPOS,
	)
	ProtClusteringGroup.add_option(	'--MUSCLE_GapExtendCost',
					dest	= "MUSCLE_GapExtendCost",
					default	= -0.00,
					help	= "MUSCLE gap extension panelty for aligning protein sequences. [default: %default]",
					metavar	= "gap extension panelty",
					type	= "float",
					action	= "callback",
					callback= check_NONPOS,
	)
	
	ProtClusteringGroup.add_option(	'--ProtClustering_MCLInflation',
					dest	= "ProtClustering_MCLInflation",
					default	= 2,
					help	= "Cluster granularity. Increasing INFLATION will increase cluster granularity. [default: %default]",
					metavar	= "INFLATION",
					type	= "float",
					action	= "callback",
					callback= check_POS,
	)
	
	parser.add_option_group(ProtClusteringGroup)
	
	MergingProtAlnsGroup	= optparse.OptionGroup( parser		= parser,
							title		= "Protein alignment merging options",
							description	= "Protein alignments can be merged to reduce the number of protein profile hidden markov models (PPHMMs) in the PPHMM database. "
									  "In summary, a database of PPHMMs is constructed, and all PPHMMs (queries) are searched against the database one by one using HHsuite (hhsearch, see 'HHsuite options' for more options). "
									  "PPHMM-PPHMM pairwise similarity scores are collected, and subsequently used for PPHMM clustering by MCL algorithm. "
									  "Sequneces in the same PPHMM cluster are then re-aligned using MUSCLE.")
	MergingProtAlnsGroup.add_option('--N_AlignmentMerging',
					dest	= "N_AlignmentMerging",
					default	= 0,
					help	= "Number of rounds of alignment merging. ROUND == 0 means no merging. ROUND == -1 means merging until exhausted. [default: %default]",
					metavar	= "ROUND",
					type	= "int",
					action	= "callback",
					callback= check_N_AlignmentMerging,
	)
	
	MergingProtAlnsGroup.add_option('--PPHMMClustering_MCLInflation_ForAlnMerging',
					dest	= "PPHMMClustering_MCLInflation_ForAlnMerging",
					default	= 5,
					help	= "Cluster granularity. Increasing INFLATION will increase cluster granularity. [default: %default]",
					metavar	= "INFLATION",
					type	= "float",
					action	= "callback",
					callback= check_POS,
	)
	
	MergingProtAlnsGroup.add_option('--HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging',
					dest	= "HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging",
					default	= True,
					help	= "Make a HMMER PPHMM DB for each round of protein merging if True. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	parser.add_option_group(MergingProtAlnsGroup)
	
	HHsuiteGroup	= optparse.OptionGroup( parser		= parser,
						title		= "HHsuite options",
						description	= "HHsuite (hhsearch) is used to search PPHMMs against a database of PPHMMs. "
								  "Here are some options regarding hhsearch. ")
	HHsuiteGroup.add_option('--HHsuite_evalue_Cutoff',
				dest	= "HHsuite_evalue_Cutoff",
				default	= 1E-6,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with an E-value > E-VALUE will be ignored. [default: %default]",
				metavar = "E-VALUE",
				type	= "float",
				action	= "callback",
				callback= check_POS,
	)
	HHsuiteGroup.add_option('--HHsuite_pvalue_Cutoff',
				dest	= "HHsuite_pvalue_Cutoff",
				default	= 0.05,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with a p-value > P-VALUE will be ignored. [default: %default]",
				metavar	= "P-VALUE",
				type	= "float",
				action	= "callback",
				callback= check_PROB,
	)
	HHsuiteGroup.add_option('--HHsuite_N_CPUs',
				dest	= "HHsuite_N_CPUs",
				default	= multiprocessing.cpu_count(),
				help	= "Number of threads (CPUs) to use in the hhsearch search. [default: %default - all threads]",
				metavar	= "THREADS",
				type	= "int",
				action	= "callback",
				callback= check_POSINTEGER,
	)
	HHsuiteGroup.add_option('--HHsuite_QueryCoverage_Cutoff',
				dest	= "HHsuite_QueryCoverage_Cutoff",
				default	= 85,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with a query coverage < COVERAGE will be ignored [default: %default]",
				metavar	= "COVERAGE",
				type	= "float",
				action	= "callback",
				callback= check_PERCENT,
	)
	HHsuiteGroup.add_option('--HHsuite_SubjectCoverage_Cutoff',
				dest	= "HHsuite_SubjectCoverage_Cutoff",
				default	= 85,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with a subject coverage < COVERAGE will be ignored [default: %default]",
				metavar	= "COVERAGE",
				type	= "float",
				action	= "callback",
				callback= check_PERCENT,
	)
	parser.add_option_group(HHsuiteGroup)
	
	VirusAnnotationGroup 	= optparse.OptionGroup(	parser		= parser,
							title		= "Reference virus annotation options",
							description	= "Virus genomes are 6-framed translated and are scanned against the PPHMM database, using HMMER (hmmscan, see 'HMMER options' for more options). "
									  "Two types of information are collected: PPHMM hit scores (PPHMM signatures) and hit locations (PPHMM LOCATION signatures).")

	VirusAnnotationGroup.add_option('--AnnotateIncompleteGenomes',
					dest	= "AnnotateIncompleteGenomes",
					default	= False,
					help	= "Annotate all viral genomes if True, otherwise only complete genomes. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	parser.add_option_group(VirusAnnotationGroup)
	
	HMMERGroup	= optparse.OptionGroup( parser		= parser,
						title		= "HMMER options",
						description	= "HMMER (hmmscan) is used to search virus (translated) genomes against a database of PPHMMs. "
								  "Here are some options regarding hmmscan. ")
	
	HMMERGroup.add_option(	'--HMMER_N_CPUs',
				dest	= "HMMER_N_CPUs",
				default	= multiprocessing.cpu_count(),
				help	= "Number of threads (CPUs) to use in the hmmscan search. [default: %default - all threads]",
				metavar	= "THREADS",
				type	= "int",
				action	= "callback",
				callback= check_POSINTEGER,
	)
	HMMERGroup.add_option(	'--HMMER_C_EValue_Cutoff',
				dest	= "HMMER_C_EValue_Cutoff",
				default	= 1E-3,
				help	= "Threshold for HMM-protein similarity detection. "
					  "A hit with an E-value > E-VALUE will be ignored. [default: %default]",
				metavar	= "E-VALUE",
				type	= "float",
				action	= "callback",
				callback= check_POS,
	)
	HMMERGroup.add_option(	'--HMMER_HitScore_Cutoff',
				dest	= "HMMER_HitScore_Cutoff",
				default	= 0,
				help	= "Threshold for HMM-protein similarity detection. "
					  "A hit with a score < SCORE will be ignored. [default: %default]",
				metavar	= "SCORE",
				type	= "float",
				action	= "callback",
				callback= check_NONNEG,
	)	
	parser.add_option_group(HMMERGroup)
	
	RemoveSingletonPPHMMGroup 	= optparse.OptionGroup(	parser		= parser,
								title		= "Remove singleton PPHMM options",
								description	= "Remove singleton PPHMMs (PPHMM that show similarity to only one virus) from the database. "
										  "Note that some singleton PPHMMs may be informative, in particular those that show similarity to the only representative of its taxonomic group (N = 1). "
										  "Sometimes, a taxonomic group only has a few members. Singleton PPHMMs associated with those taxonomic groups can also be informative. "
										  "Some viruses (which may belong to a large taxonomic group) may exhibit similarity to singleton PPHMMs excluively. "
										  "Thus removing singleton PPHMMs can be dangerous. "
										  "This action modifies the PPHMM database permanently. Use with caution. "
										  )
	RemoveSingletonPPHMMGroup.add_option(	'--RemoveSingletonPPHMMs',
						dest	= "RemoveSingletonPPHMMs",
						default	= False,
						help	= "Remove singleton PPHMMs from the database if True. [default: %default]",
						metavar	= "BOOLEAN",
						type	= "choice",
						choices	= ["True", "False",],
	)
	RemoveSingletonPPHMMGroup.add_option(	'--N_VirusesOfTheClassToIgnore',
						dest	= "N_VirusesOfTheClassToIgnore",
						default	= 1,
						help	= "When 'RemoveSingletonPPHMMs' == TRUE, singleton PPHMMs are removed from the PPHMM database "
							  "only if they show similarity to viruses that belong to taxonomic groups with more than NUMBER members. [default: %default]",
						metavar	= "NUMBER",
						type	= "int",
						action	= "callback",
						callback= check_POSINTEGER,
	)
	parser.add_option_group(RemoveSingletonPPHMMGroup)
	
	PPHMMSortingGroup 	= optparse.OptionGroup(	parser		= parser,
							title		= "PPHMM sorting options",
							description	= "Sort PPHMMs in the table and the database based on their pairiwse similarity and their presence/absence patterns in the viruses. "
									  "PPHMM pairwise similarity scores were determined by HHsuite (hhsearch, see 'HHsuite options' for more options), and PPHMM clustering was performed by the MCL algorithm. "
									  "This action modifies the PPHMM database permanently. Use with caution. ")
	PPHMMSortingGroup.add_option(	'--PPHMMSorting',
					dest	= "PPHMMSorting",
					default	= False,
					help	= "Sort PPHMMs if True. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)	
	PPHMMSortingGroup.add_option(	'--PPHMMClustering_MCLInflation_ForPPHMMSorting',
					dest	= "PPHMMClustering_MCLInflation_ForPPHMMSorting",
					default	= 2,
					help	= "Cluster granularity. Increasing INFLATION will increase cluster granularity. [default: %default]",
					metavar	= "INFLATION",
					type	= "float",
					action	= "callback",
					callback= check_POS,
	)
	parser.add_option_group(PPHMMSortingGroup)
	
	SimilarityMeasurementGroup = optparse.OptionGroup(	parser		= parser,
								title		= "Virus (dis)similarity measurement options",
								description	= "Unlike in many traditional approaches that represent a virus using its molecular sequences, "
										  "in GRAViTy, a virus is represented by 3 signatures: a PPHMM signature, a PPHMM LOCATION signature, and a GOM signature. "
										  "For each virus pair, generalised Jaccard similarity scores (GJ) between their PPHMM signatures (GJ_P), PPHMM LOCATION signatures (GJ_L), and GOM signatures (GJ_G) are computed. "
										  "An overall similarity between a pair of viruses is computed based on these 3 quantities. "
										  "Various schemes of overall similarity measurement are implemented in GRAViTy, specified by using 'SimilarityMeasurementScheme' argument. "
										  "See help text of 'SimilarityMeasurementScheme' for more details. "
										  "The degree of overall similarity (S) is ranging between 0 and 1. "
										  "An overall distance (D) between a pair of viruses is 1 - S, thus again ranging between 0 and 1."
										  "In addition, GRAViTy allows users to transform the overall distance as well, governed by parameter 'p'. "
										  "See the help text of 'p' for more details.",
										  )
	SimilarityMeasurementGroup.add_option(	'--SimilarityMeasurementScheme',
						dest	= "SimilarityMeasurementScheme",
						default	= "PG",
						help	= "Virus similarity measurement SCHEMEs. [default: %default]\n "
							  "If SCHEME = 'P', an overall similarity between two viruses is their GJ_P.\n "
							  "If SCHEME = 'L', an overall similarity between two viruses is their GJ_L.\n "
							  "If SCHEME = 'G', an overall similarity between two viruses is their GJ_G.\n "
							  "If SCHEME = 'PG', an overall similarity between two viruses is a geometric mean - or a 'composite generalised Jaccard score' (CGJ) - of their GJ_P and GJ_G.\n "
							  "If SCHEME = 'PL', an overall similarity between two viruses is a geometric mean - or a 'composite generalised Jaccard score' (CGJ) - of their GJ_P and GJ_L.",
						metavar	= "SCHEME",
						type	= "choice",
						choices	= ["P", "G", "L", "PG", "PL"],
	)
	SimilarityMeasurementGroup.add_option(	'--p',
						dest	= "p",
						default	= 1,
						help	= "Distance transformation P coefficient, 0 <= P. [default: %default]\n "
							  "D = 1 - S**P.\n "
							  "If P = 1, no distance transformation is applied.\n "
							  "If P > 1, shallow branches will be stretched out so shallow splits can be seen more clearly.\n "
							  "If p < 1, deep branches are stretch out so that deep splits can be seen more clearly.\n "
							  "If p = 0, the entire dendrogram will be reduced to a single branch (root), with all taxa forming a polytomy clade at the tip of the dendrogram.\n "
							  "If p = Inf, the resultant dendrogram will be star-like, with all taxa forming a polytomy clade at the root.",
						metavar	= "P",
						type	= "float",
						action	= "callback",
						callback= check_NONNEG,
	)
	parser.add_option_group(SimilarityMeasurementGroup)
	
	DendrogramConstructionGroup = optparse.OptionGroup(	parser		= parser,
								title		= "Dendrogram construction options",
								description	= "GRAViTy can generate a (bootstrapped) dendrogram based on the pairwise distance matrix using a hierarchical clustering algorithm. "
										  "Various algorithms of hierarchical clustering are implemented in GRAViTy, specified by using 'Dendrogram_LinkageMethod' argument. "
										  "See the help text of 'Dendrogram_LinkageMethod' for more details.")
	DendrogramConstructionGroup.add_option(	'--Dendrogram',
						dest	= "Dendrogram",
						default	= True,
						help	= "Construct dendrogram if True. [default: %default]",
						metavar	= "BOOLEAN",
						type	= "choice",
						choices	= ["True", "False",],
	)
	DendrogramConstructionGroup.add_option(	'--Dendrogram_LinkageMethod',
						dest	= "Dendrogram_LinkageMethod",
						default	= "average",
						help	= "LINKAGE for dendrogram construction. [default: %default]\n "
							  "If LINKAGE = 'single', the nearest point algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'complete', the farthest point algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'average', the UPGMA algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'weighted', the WPGMA algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'centroid', the UPGMC algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'median', the WPGMC algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'ward', the incremental algorithm is used to cluster viruses and compute cluster distances.",
						metavar	= "LINKAGE",
						type	= "choice",
						choices	= ["single", "complete", "average", "weighted", "centroid", "median", "ward"],
	)
	
	DendrogramConstructionGroup.add_option(	'--Bootstrap',
						dest	= "Bootstrap",
						default	= True,
						help	= "Perform bootstrapping if True. [default: %default]",
						metavar	= "BOOLEAN",
						type	= "choice",
						choices	= ["True", "False",],
	)
	DendrogramConstructionGroup.add_option(	'--N_Bootstrap',
						dest	= "N_Bootstrap",
						default	= 10,
						help	= "The number of pseudoreplicate datasets by resampling. [default: %default]",
						metavar	= "NUMBER",
						type	= "int",
						action	= "callback",
						callback= check_POSINTEGER,
	)
	DendrogramConstructionGroup.add_option(	'--Bootstrap_method',
						dest	= "Bootstrap_method",
						default	= "booster",
						help	= "Two METHODs for tree summary construction are implemented in GRAViTy. [default: %default]\n "
							  "If METHOD = 'sumtrees', SumTrees (Sukumaran, J & MT Holder, 2010, Bioinformatics; https://dendropy.org/programs/sumtrees.html) will be used to summarize non-parameteric bootstrap support for splits on the best estimated dendrogram. The calculation is based on the standard Felsenstein bootstrap method.\n "
							  "If METHOD = 'booster', BOOSTER (Lemoine et al., 2018, Nature; https://booster.pasteur.fr/) will be used. With large trees and moderate phylogenetic signal, BOOSTER tends to be more informative than the standard Felsenstein bootstrap method.",
						metavar	= "METHOD",
						type	= "choice",
						choices	= ["booster", "sumtrees"],
	)
	DendrogramConstructionGroup.add_option(	'--Bootstrap_N_CPUs',
						dest	= "Bootstrap_N_CPUs",
						default	= multiprocessing.cpu_count(),
						help	= "Number of threads (CPUs) to use in tree summary. Only used when 'Bootstrap_method' == 'booster' [default: %default - all threads]",
						metavar	= "THREADS",
						type	= "int",
						action	= "callback",
						callback= check_POSINTEGER,
	)
	parser.add_option_group(DendrogramConstructionGroup)
	
	HeatmapConstructionGroup = optparse.OptionGroup(	parser		= parser,
								title		= "Heatmap construction options",
								description	= "GRAViTy can generate a heatmap (with the dendrogram) to represent the pairwise (dis)similarity matrix")
	HeatmapConstructionGroup.add_option(	'--Heatmap',
						dest	= "Heatmap",
						default	= False,
						help	= "Construct (dis)similarity heatmap if True. [default: %default]",
						metavar	= "BOOLEAN",
						type	= "choice",
						choices	= ["True", "False",],
	)
	HeatmapConstructionGroup.add_option(	'--Heatmap_VirusOrderScheme',
						dest	= "Heatmap_VirusOrderScheme",
						default	= None,
						help	= "Full path to the virus order file. The indices of the genome entries start from 0 [default: %default].",
						metavar	= "FILEPATH",
						type	= "string",
						action	= "callback",
						callback= check_FILEPATH,
	)
	
	HeatmapConstructionGroup.add_option(	'--Heatmap_WithDendrogram',
						dest	= "Heatmap_WithDendrogram",
						default	= True,
						help	= "Construct (dis)similarity heatmap with dendrogram if True. [default: %default]",
						metavar	= "BOOLEAN",
						type	= "choice",
						choices	= ["True", "False",],
	)
	HeatmapConstructionGroup.add_option(	'--Heatmap_DendrogramFile',
						dest	= "Heatmap_DendrogramFile",
						default	= None,
						help	= "Full path to the dendrogram file. If 'None', the dendrogram will be estimated by GRAViTy [default: %default]",
						metavar = "FILE",
						type	= "string",
						action	= "callback",
						callback= check_FILEPATH,
	)
	HeatmapConstructionGroup.add_option(	'--Heatmap_DendrogramSupport_Cutoff',
						dest	= "Heatmap_DendrogramSupport_Cutoff",
						default	= 0.75,
						help	= "Threshold for the BOOTSTRAP SUPPORT to be shown on the dendrogram on the heatmap. [default: %default]",
						metavar	= "BOOTSTRAP SUPPORT",
						type	= "float",
						action	= "callback",
						callback= check_PROB,
	)
	parser.add_option_group(HeatmapConstructionGroup)
	
	VirusGroupingGroup = optparse.OptionGroup(	parser		= parser,
							title		= "Virus grouping options",
							description	= "GRAViTy can estimate the distance cutoff that best separates the input reference taxonomic groupings, and report virus groups suggested by the estimated cutoff. "
									  "Theil's uncertainty correlation for the reference taxonomic grouping given the predicted grouping, and vice versa, are reported. "
									  "Symmetrical Theil's uncertainty correlation between the reference and predicted taxonomic grouping are also reported.")
	VirusGroupingGroup.add_option(	'--VirusGrouping',
					dest	= "VirusGrouping",
					default	= True,
					help	= "Perform virus grouping if True. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	parser.add_option_group(VirusGroupingGroup)
	
	VirusGroupingForMICalGroup = optparse.OptionGroup(	parser		= parser,
								title		= "Virus grouping for mutual information calculation options",
								description	= "GRAViTy can calculate mutual information between (various schemes of) taxonomic groupings and values of PPHMM scores "
										  "to determine which PPHMMs are highly (or weakly) correlated with the virus taxonomic scheme(s).")
	VirusGroupingForMICalGroup.add_option(	'--VirusGroupingSchemesFile',
						dest	= "VirusGroupingFile",
						default	= None,
						help	= "Fill path to the virus grouping scheme file. [default: %default] "
							  "The file contains column(s) of arbitrary taxonomic grouping scheme(s) that users want to investigate. The file may look something like: \n\n"
							  "Scheme 1\tScheme 2\t...\n"
							  "A\tX\t...\n"
							  "A\tX\t...\n"
							  "B\tX\t...\n"
							  "B\tX\t...\n"
							  "C\tY\t...\n"
							  "C\tY\t...\n"
							  "...\t...\t...\n\n"
							  "If 'None', the taxonomic grouping as specified in 'Taxonomic grouping' column in the VMR will be used. "
							  "Note that the file must contain headers.",
						metavar	= "FILEPATH",
						type	= "string",
						action	= "callback",
						callback= check_FILEPATH,
	)
	parser.add_option_group(VirusGroupingForMICalGroup)
	
	SamplingGroup = optparse.OptionGroup(	parser		= parser,
						title		= "Virus sampling options",
						description	= "")
	SamplingGroup.add_option(	'--N_Sampling',
					dest = "N_Sampling",
					default = 10,
					help = "The number of mutual information scores sample size. [default: %default]",
					metavar = "NUMBER",
					type = "int",
					action = "callback",
					callback = check_POSINTEGER,
	)
	SamplingGroup.add_option(	'--SamplingStrategy',
					dest = "SamplingStrategy",
					default = "balance_with_repeat",
					help = "Virus sampling scheme. [default: %default]",
					metavar = "SCHEME",
					type = "choice",
					choices = [None, "balance_without_repeat", "balance_with_repeat"],
	)
	SamplingGroup.add_option(	'--SampleSizePerGroup',
					dest = "SampleSizePerGroup",
					default = 10,
					help = "If 'SamplingStrategy' != None, this option specifies the number of viruses to be sampled per taxonomic group [default: %default]",
					metavar = "NUMBER",
					type = "int",
					action = "callback",
					callback = check_POSINTEGER,
	)
	parser.add_option_group(SamplingGroup)
	
	return parser

def generate_pipeline_ii_arguments():
	parser = optparse.OptionParser(	usage = "usage: %prog [options]",
					version = "%prog 1.1.0",
					formatter = IndentedHelpFormatterWithNL())
	
	parser.add_option(	'--GenomeDescTableFile_UcfVirus',
				dest 	= "GenomeDescTableFile_UcfVirus",
				help 	= "Full path to the Virus Metadata Resource-like (VMR-like) tab delimited file of unclassified viruses, wth headers. "
					  "VMR can be downloaded from https://talk.ictvonline.org/taxonomy/vmr/",
				metavar	= "FILEPATH",
				type	= "string",
				action	= "callback",
				callback= check_FILEPATH,
	)
	parser.add_option(	'--ShelveDir_UcfVirus',
				dest	= "ShelveDir_UcfVirus",
				help	= "Full path to the shelve directory of unclassified viruses, storing GRAViTy outputs.",
				metavar	= "DIRECTORYPATH",
				type	= "string",
	)
	parser.add_option(	'--ShelveDirs_RefVirus',
				dest	= "ShelveDirs_RefVirus",
				help	= "Full path(s) to the shelve director(y/ies) of reference virus(es). "
					  "For example: 'path/to/shelve/ref1, path/to/shelve/ref2, ...'",
				metavar	= "DIRECTORYPATH",
				type	= "string",
	)
	
	parser.add_option(	'--Database',
				dest	= "Database",
				default	= None,
				help	= "GRAViTy will only analyse genomes that are labelled with DATABASE in the database column. "
					  "The database column can be specified by the DATABASE HEADER argument. "
					  "If 'None', all entries are analysed. [default: %default]",
				metavar	= "DATABASE",
				type	= "string",
	)
	parser.add_option(	'--Database_Header',
				dest	= "Database_Header",
				default	= None,
				help	= "The header of the database column. Cannot be 'None' if DATABASE is specified. [default: %default]",
				metavar	= "DATABASE COLUMN HEADER",
				type	= "string",
	)

	parser.add_option(	'--GenomeSeqFile_UcfVirus',
				dest	= "GenomeSeqFile_UcfVirus",
				help	= "Full path to the genome sequence GenBank file of unclassified viruses.",
				metavar	= "FILEPATH",
				type	= "string",
	)
	parser.add_option(	'--UseUcfVirusPPHMMs',
				dest	= "UseUcfVirusPPHMMs",
				default	= True,
				help	= "Annotate reference and unclassified viruses using the PPHMM database derived from unclassified viruses if True. [default: %default]",
				metavar	= "BOOLEAN",
				type	= "choice",
				choices	= ["True", "False",],
	)
	parser.add_option(	'--GenomeSeqFiles_RefVirus',
				dest	= "GenomeSeqFiles_RefVirus",
				default	= None,
				help	= "Full path(s) to the genome sequence GenBank file(s) of reference viruses. "
					  "For example: 'path/to/GenBank/ref1, path/to/GenBank/ref2, ...' "
					  "This cannot be 'None' if UseUcfVirusPPHMMs = True. ",
				metavar	= "FILEPATH",
				type	= "string",
				action	= "callback",
				callback= check_FILEPATHS,
	)
	
	ProtExtractionGroup 	= optparse.OptionGroup(	parser		= parser,
							title		= "Protein extraction options",
							description	= "Protein sequences are extracted from the genome sequence GenBank File. "
									  "If the file contains protein sequences, GRAViTy will use those annotations. "
									  "Otherwise, protein sequences will be inferred using the genetic code specified in the VMR file, 'Genetic code table' column")
	ProtExtractionGroup.add_option(	'--ProteinLength_Cutoff',
					dest	= "ProteinLength_Cutoff",
					default	= 100,
					help	= "Proteins with length < LENGTH aa will be ignored [default: %default]",
					metavar	= "LENGTH",
					type	= "int",
					action	= "callback",
					callback= check_NONNEGINTEGER,
	)
	ProtExtractionGroup.add_option(	'--IncludeProteinsFromIncompleteGenomes',
					dest	= "IncludeProteinsFromIncompleteGenomes",
					default	= True,
					help	= "Include protein sequences from incomplete genomes to the database if True. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	parser.add_option_group(ProtExtractionGroup)
	
	ProtClusteringGroup	= optparse.OptionGroup( parser		= parser,
							title		= "Protein clustering options",
							description	= "Protein sequences are clustered based on ALL-VERSUS-ALL BLASTp hit scores. "
									  "In summary, a protein database is constructed, and all protein sequences (queries) are searched against the database one by one using BLASTp. "
									  "Protein-protein pairwise similarity scores (bit scores) are collected, and subsequently used for protein clustering by the MCL algorithm. "
									  "Sequneces in each cluster are then aligned by using MUSCLE, and turned into a protein profile hidden Markov model (PPHMM).")
	ProtClusteringGroup.add_option(	'--BLASTp_evalue_Cutoff',
					dest	= "BLASTp_evalue_Cutoff",
					default	= 1E-3,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with an E-value > E-VALUE will be ignored. [default: %default]",
					metavar	= "E-VALUE",
					type	= "float",
					action	= "callback",
					callback= check_POS,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_PercentageIden_Cutoff',
					dest	= "BLASTp_PercentageIden_Cutoff",
					default	= 50,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with a percentage identity < PERCENTAGE IDENTITY will be ignored. [default: %default]",
					metavar	= "PERCENTAGE IDENTITY",
					type	= "float",
					action	= "callback",
					callback= check_PERCENT,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_QueryCoverage_Cutoff',
					dest	= "BLASTp_QueryCoverage_Cutoff",
					default	= 75,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with a query coverage < COVERAGE will be ignored. [default: %default]",
					metavar	= "COVERAGE",
					type	= "float",
					action	= "callback",
					callback= check_PERCENT,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_SubjectCoverage_Cutoff',
					dest	= "BLASTp_SubjectCoverage_Cutoff",
					default	= 75,
					help	= "Threshold for protein sequence similarity detection. "
						  "A hit with a subject coverage < COVERAGE will be ignored. [default: %default]",
					metavar	= "COVERAGE",
					type	= "float",
					action	= "callback",
					callback= check_PERCENT,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_num_alignments',
					dest	= "BLASTp_num_alignments",
					default	= 1000000,
					help	= "Maximum number of sequences to be considered in a BLASTp search. [default: %default]",
					metavar	= "BLAST ALIGNMENT NUMBER",
					type	= "int",
					action	= "callback",
					callback= check_POSINTEGER,
	)
	ProtClusteringGroup.add_option(	'--BLASTp_N_CPUs',
					dest	= "BLASTp_N_CPUs",
					default	= multiprocessing.cpu_count(),
					help	= "The number of threads (CPUs) to use in the BLASTp search. [default: %default - all threads]",
					metavar	= "THREADS",
					type	= "int",
					action	= "callback",
					callback= check_POSINTEGER,
	)
	
	ProtClusteringGroup.add_option(	'--MUSCLE_GapOpenCost',
					dest	= "MUSCLE_GapOpenCost",
					default	= -3.0,
					help	= "MUSCLE gap opening panelty for aligning protein sequences. [default: %default]",
					metavar	= "gap opening panelty",
					type	= "float",
					action	= "callback",
					callback= check_NONPOS,
	)
	ProtClusteringGroup.add_option(	'--MUSCLE_GapExtendCost',
					dest	= "MUSCLE_GapExtendCost",
					default	= -0.00,
					help	= "MUSCLE gap extension panelty for aligning protein sequences. [default: %default]",
					metavar	= "gap extension panelty",
					type	= "float",
					action	= "callback",
					callback= check_NONPOS,
	)
	
	ProtClusteringGroup.add_option(	'--ProtClustering_MCLInflation',
					dest	= "ProtClustering_MCLInflation",
					default	= 2,
					help	= "Cluster granularity. Increasing INFLATION will increase cluster granularity. [default: %default]",
					metavar	= "INFLATION",
					type	= "float",
					action	= "callback",
					callback= check_POS,
	)
	parser.add_option_group(ProtClusteringGroup)
	
	MergingProtAlnsGroup	= optparse.OptionGroup( parser		= parser,
							title		= "Protein alignment merging options",
							description	= "Protein alignments can be merged to reduce the number of protein profile hidden markov models (PPHMMs) in the PPHMM database. "
									  "In summary, a database of PPHMMs is constructed, and all PPHMMs (queries) are searched against the database one by one using HHsuite (hhsearch, see 'HHsuite options' for more options). "
									  "PPHMM-PPHMM pairwise similarity scores are collected, and subsequently used for PPHMM clustering by MCL algorithm. "
									  "Sequneces in the same PPHMM cluster are then re-aligned using MUSCLE.")
	MergingProtAlnsGroup.add_option('--N_AlignmentMerging',
					dest	= "N_AlignmentMerging",
					default	= 0,
					help	= "Number of rounds of alignment merging. ROUND == 0 means no merging. ROUND == -1 means merging until exhausted. [default: %default]",
					metavar	= "ROUND",
					type	= "int",
					action	= "callback",
					callback= check_N_AlignmentMerging,
	)
	
	MergingProtAlnsGroup.add_option('--PPHMMClustering_MCLInflation_ForAlnMerging',
					dest	= "PPHMMClustering_MCLInflation_ForAlnMerging",
					default	= 5,
					help	= "Cluster granularity. Increasing INFLATION will increase cluster granularity. [default: %default]",
					metavar	= "INFLATION",
					type	= "float",
					action	= "callback",
					callback= check_POS,
	)
	
	MergingProtAlnsGroup.add_option('--HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging',
					dest	= "HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging",
					default	= True,
					help	= "Make a HMMER PPHMM DB for each round of protein merging if True. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	parser.add_option_group(MergingProtAlnsGroup)
	
	HHsuiteGroup	= optparse.OptionGroup( parser		= parser,
						title		= "HHsuite options",
						description	= "HHsuite (hhsearch) is used to search PPHMMs against a database of PPHMMs. "
								  "Here are some options regarding hhsearch. ")
	HHsuiteGroup.add_option('--HHsuite_evalue_Cutoff',
				dest	= "HHsuite_evalue_Cutoff",
				default	= 1E-6,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with an E-value > E-VALUE will be ignored. [default: %default]",
				metavar = "E-VALUE",
				type	= "float",
				action	= "callback",
				callback= check_POS,
	)
	HHsuiteGroup.add_option('--HHsuite_pvalue_Cutoff',
				dest	= "HHsuite_pvalue_Cutoff",
				default	= 0.05,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with a p-value > P-VALUE will be ignored. [default: %default]",
				metavar	= "P-VALUE",
				type	= "float",
				action	= "callback",
				callback= check_PROB,
	)
	HHsuiteGroup.add_option('--HHsuite_N_CPUs',
				dest	= "HHsuite_N_CPUs",
				default	= multiprocessing.cpu_count(),
				help	= "Number of threads (CPUs) to use in the hhsearch search. [default: %default - all threads]",
				metavar	= "THREADS",
				type	= "int",
				action	= "callback",
				callback= check_POSINTEGER,
	)
	HHsuiteGroup.add_option('--HHsuite_QueryCoverage_Cutoff',
				dest	= "HHsuite_QueryCoverage_Cutoff",
				default	= 85,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with a query coverage < COVERAGE will be ignored [default: %default]",
				metavar	= "COVERAGE",
				type	= "float",
				action	= "callback",
				callback= check_PERCENT,
	)
	HHsuiteGroup.add_option('--HHsuite_SubjectCoverage_Cutoff',
				dest	= "HHsuite_SubjectCoverage_Cutoff",
				default	= 85,
				help	= "Threshold for PPHMM similarity detection. "
					  "A hit with a subject coverage < COVERAGE will be ignored [default: %default]",
				metavar	= "COVERAGE",
				type	= "float",
				action	= "callback",
				callback= check_PERCENT,
	)
	parser.add_option_group(HHsuiteGroup)
	
	VirusAnnotationUsingRefDBGroup 	= optparse.OptionGroup(	parser		= parser,
								title		= "Unclassified virus annotation options",
								description	= "Virus genomes are 6-framed translated and are scanned against the reference PPHMM database(s), using HMMER (hmmscan). "
										  "Two types of information are collected: PPHMM hit scores (PPHMM signatures) and hit locations (PPHMM LOCATION signatures). ")
	
	VirusAnnotationUsingRefDBGroup.add_option('--AnnotateIncompleteGenomes_UcfVirus',
					dest	= "AnnotateIncompleteGenomes_UcfVirus",
					default	= True,
					help	= "Annotate all unclassified viruses using reference PPHMM database(s) if True, otherwise only complete genomes. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	VirusAnnotationUsingRefDBGroup.add_option('--UsingDatabaseIncludingIncompleteRefViruses',
					dest	= "UsingDatabaseIncludingIncompleteRefViruses",
					default	= False,
					help	= "Annotate unclassified viruses using the PPHMM and GOM databases derived from all reference viruses if True, "
						  "otherwise using those derived from complete reference genomes only. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)

	parser.add_option_group(VirusAnnotationUsingRefDBGroup)
	
	HMMERGroup	= optparse.OptionGroup( parser		= parser,
						title		= "HMMER options",
						description	= "HMMER (hmmscan) is used to search virus (translated) genomes against a database of PPHMMs. "
								  "Here are some options regarding hmmscan. ")
	
	HMMERGroup.add_option(	'--HMMER_N_CPUs',
				dest	= "HMMER_N_CPUs",
				default	= multiprocessing.cpu_count(),
				help	= "Number of threads (CPUs) to use in the hmmscan search. [default: %default - all threads]",
				metavar	= "THREADS",
				type	= "int",
				action	= "callback",
				callback= check_POSINTEGER,
	)
	HMMERGroup.add_option(	'--HMMER_C_EValue_Cutoff',
				dest	= "HMMER_C_EValue_Cutoff",
				default	= 1E-3,
				help	= "Threshold for HMM-protein similarity detection. "
					  "A hit with an E-value > E-VALUE will be ignored. [default: %default]",
				metavar	= "E-VALUE",
				type	= "float",
				action	= "callback",
				callback= check_POS,
	)
	HMMERGroup.add_option(	'--HMMER_HitScore_Cutoff',
				dest	= "HMMER_HitScore_Cutoff",
				default	= 0,
				help	= "Threshold for HMM-protein similarity detection. "
					  "A hit with a score < SCORE will be ignored. [default: %default]",
				metavar	= "SCORE",
				type	= "float",
				action	= "callback",
				callback= check_NONNEG,
	)	
	parser.add_option_group(HMMERGroup)
	
	SimilarityMeasurementGroup = optparse.OptionGroup(	parser		= parser,
								title		= "Virus (dis)similarity measurement options",
								description	= "Unlike in many traditional approaches that represent a virus using its molecular sequences, "
										  "in GRAViTy, a virus is represented by 3 signatures: a PPHMM signature, a PPHMM LOCATION signature, and a GOM signature. "
										  "For each virus pair, generalised Jaccard similarity scores (GJ) between their PPHMM signatures (GJ_P), PPHMM LOCATION signatures (GJ_L), and GOM signatures (GJ_G) are computed. "
										  "An overall similarity between a pair of viruses is computed based on these 3 quantities. "
										  "Various schemes of overall similarity measurement are implemented in GRAViTy, specified by using 'SimilarityMeasurementScheme' argument. "
										  "See help text of 'SimilarityMeasurementScheme' for more details. "
										  "The degree of overall similarity (S) is ranging between 0 and 1. "
										  "An overall distance (D) between a pair of viruses is 1 - S, thus again ranging between 0 and 1."
										  "In addition, GRAViTy allows users to transform the overall distance as well, governed by parameter 'p'. "
										  "See the help text of 'p' for more details.",
										  )
	SimilarityMeasurementGroup.add_option(	'--SimilarityMeasurementScheme',
						dest	= "SimilarityMeasurementScheme",
						default	= "PG",
						help	= "Virus similarity measurement SCHEMEs. [default: %default]\n "
							  "If SCHEME = 'P', an overall similarity between two viruses is their GJ_P.\n "
							  "If SCHEME = 'L', an overall similarity between two viruses is their GJ_L.\n "
							  "If SCHEME = 'G', an overall similarity between two viruses is their GJ_G.\n "
							  "If SCHEME = 'PG', an overall similarity between two viruses is a geometric mean - or a 'composite generalised Jaccard score' (CGJ) - of their GJ_P and GJ_G.\n "
							  "If SCHEME = 'PL', an overall similarity between two viruses is a geometric mean - or a 'composite generalised Jaccard score' (CGJ) - of their GJ_P and GJ_L.",
						metavar	= "SCHEME",
						type	= "choice",
						choices	= ["P", "G", "L", "PG", "PL"],
	)
	SimilarityMeasurementGroup.add_option(	'--p',
						dest	= "p",
						default	= 1,
						help	= "Distance transformation P coefficient, 0 <= P. [default: %default]\n "
							  "D = 1 - S**P.\n "
							  "If P = 1, no distance transformation is applied.\n "
							  "If P > 1, shallow branches will be stretched out so shallow splits can be seen more clearly.\n "
							  "If p < 1, deep branches are stretch out so that deep splits can be seen more clearly.\n "
							  "If p = 0, the entire dendrogram will be reduced to a single branch (root), with all taxa forming a polytomy clade at the tip of the dendrogram.\n "
							  "If p = Inf, the resultant dendrogram will be star-like, with all taxa forming a polytomy clade at the root.",
						metavar	= "P",
						type	= "float",
						action	= "callback",
						callback= check_NONNEG,
	)
	parser.add_option_group(SimilarityMeasurementGroup)
	
	VirusClassificationGroup = optparse.OptionGroup(	parser		= parser,
								title		= "Virus classification options",
								description	= "In summary, GRAViTy proposes a candidate taxonomic group for an unclassified virus by identifying the most similar virus(es) in the refernece database(s). "
										  "To validate the candidate taxonomic assignment, GRAViTy employs two-step evaluation protocol.\n\n"
										  
										  "In the first step, GRAViTy checks whether or not the unclassified virus is 'similar enough' to the proposed candidate group, of which the CGJ similarity threshold is group specific. "
										  "To estimate the threshold for a particular taxonomic group, GRAViTy builds distributions of its inter-group and intra-group CGJ similarity scores. "
										  "GRAViTy then computes the score that best separate the two distributions using the support vector machine (SVM) algorithm. "
										  "The SVM algorithm used is implemented in SVC function, available from Scikit-learn python library, with 'balanced' class weight option. "
										  "If the observed CGJ similarity is less than the threshold, the candidate taxonomic assignment is rejected, and the sample is relabelled as 'Unclassified', "
										  "otherwise, the second step of the evaluation will be employed to further evaluate the candidate taxonomic assignment.\n\n"
										  
										  "In the second step, a dendrogram containing reference viruses and the unclassified virus is used, and the evaluator will look at its neighbourhood. "
										  "The taxonomic proposal will be accepted if any of the following conditions are met:\n"
										  "\ti)	the sister clade is composed entirely of the members of the proposed candidate taxonomic group\n"
										  "\tii)	the immediate out group is composed entirely of the members of the proposed candidate taxonomic group\n"
										  "\tiii)	one of the two basal branches of its sister clade leading to a clade that is composed entirely of the members of the proposed candidate taxonomic group\n"
										  "To best estimate the placement of viruses, if multiple unclassified viruses are to be analysed at the same time, "
										  "a dendrogram containing all unclassified viruses will be used. "
										  "Furthermore, since there can be multiple reference databases, there are possibilities that a virus might be assigned to multiple taxonomic groups belonging to different databases. "
										  "In such cases, the finalised taxonomic assignment is the one associated with the highest CGJ similarity score.\n\n"
										  
										  "GRAViTy can also perform bootstrapping to evaluate the uncertainty of the proposed taxonomic group.")
	VirusClassificationGroup.add_option(	'--Dendrogram_LinkageMethod',
						dest	= "Dendrogram_LinkageMethod",
						default	= "average",
						help	= "LINKAGE for dendrogram construction. [default: %default]\n "
							  "If LINKAGE = 'single', the nearest point algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'complete', the farthest point algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'average', the UPGMA algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'weighted', the WPGMA algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'centroid', the UPGMC algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'median', the WPGMC algorithm is used to cluster viruses and compute cluster distances.\n "
							  "If LINKAGE = 'ward', the incremental algorithm is used to cluster viruses and compute cluster distances.",
						metavar	= "LINKAGE",
						type	= "choice",
						choices	= ["single", "complete", "average", "weighted", "centroid", "median", "ward"],
	)
	VirusClassificationGroup.add_option(	'--Bootstrap',
						dest	= "Bootstrap",
						default	= True,
						help	= "Perform bootstrapping if True. [default: %default]",
						metavar	= "BOOLEAN",
						type	= "choice",
						choices	= ["True", "False",],
	)
	VirusClassificationGroup.add_option(	'--N_Bootstrap',
						dest	= "N_Bootstrap",
						default	= 10,
						help	= "The number of pseudoreplicate datasets by resampling. [default: %default]",
						metavar	= "NUMBER",
						type	= "int",
						action	= "callback",
						callback= check_POSINTEGER,
	)
	VirusClassificationGroup.add_option(	'--Bootstrap_method',
						dest	= "Bootstrap_method",
						default	= "booster",
						help	= "Two METHODs for dendrogram summary construction are implemented in GRAViTy. [default: %default] "
							  "If METHOD = 'sumtrees', SumTrees (Sukumaran, J & MT Holder, 2010, Bioinformatics; https://dendropy.org/programs/sumtrees.html) will be used to summarize non-parameteric bootstrap support for splits on the best estimated dendrogram. The calculation is based on the standard Felsenstein bootstrap method.\n "
							  "If METHOD = 'booster', BOOSTER (Lemoine et al., 2018, Nature; https://booster.pasteur.fr/) will be used. With large trees and moderate phylogenetic signal, BOOSTER tends to be more informative than the standard Felsenstein bootstrap method.",
						metavar	= "METHOD",
						type	= "choice",
						choices	= ["booster", "sumtrees"],
	)
	VirusClassificationGroup.add_option(	'--Bootstrap_N_CPUs',
						dest	= "Bootstrap_N_CPUs",
						default	= multiprocessing.cpu_count(),
						help	= "Number of threads (CPUs) to use in dendrogram summary. Only used when 'Bootstrap_method' == 'booster' [default: %default - all threads]",
						metavar	= "THREADS",
						type	= "int",
						action	= "callback",
						callback= check_POSINTEGER,
	)
	
	VirusClassificationGroup.add_option(	'--DatabaseAssignmentSimilarityScore_Cutoff',
						dest	= "DatabaseAssignmentSimilarityScore_Cutoff",
						default	= 0.01,
						help	= "Threshold to determine if the unclassified virus at least belongs to a particular database. [default: %default] "
							  "For example, an unclassified virus is assigned to the family 'X' in the reference 'Baltimore group X' database, with the (greatest) similarity score of 0.1. "
							  "This score might be too low to justify that the virus is a member of the family 'X', and fail the similarity threshold test. "
							  "However, since the similarity score of 0.1 > %default, GRAViTy will make a guess that it might still be a virus of the 'Baltimore group X' database, under the default setting.",
						metavar	= "SCORE",
						type	= "int",
						action	= "callback",
						callback= check_POS,
	)
	VirusClassificationGroup.add_option(	'--N_PairwiseSimilarityScores',
						dest	= "N_PairwiseSimilarityScores",
						default	= 10000,
						help	= "Number of data points in the distributions of intra- and inter-group similarity scores used to estimate the similarity threshold. [default: %default]",
						metavar	= "NUMBER",
						type	= "int",
						action	= "callback",
						callback= check_POSINTEGER,
	)	
	parser.add_option_group(VirusClassificationGroup)
	
	HeatmapConstructionGroup = optparse.OptionGroup(	parser		= parser,
								title		= "Heatmap construction options",
								description	= "GRAViTy can generate a heatmap (with the dendrogram) to represent the pairwise (dis)similarity matrix")
	HeatmapConstructionGroup.add_option(	'--Heatmap_WithDendrogram',
						dest	= "Heatmap_WithDendrogram",
						default	= True,
						help	= "Construct (dis)similarity heatmap with dendrogram if True. [default: %default]",
						metavar	= "BOOLEAN",
						type	= "choice",
						choices	= ["True", "False",],
	)
	HeatmapConstructionGroup.add_option(	'--Heatmap_DendrogramSupport_Cutoff',
						dest	= "Heatmap_DendrogramSupport_Cutoff",
						default	= 0.75,
						help	= "Threshold for the BOOTSTRAP SUPPORT to be shown on the dendrogram on the heatmap. [default: %default]",
						metavar	= "BOOTSTRAP SUPPORT",
						type	= "float",
						action	= "callback",
						callback= check_PROB,
	)
	parser.add_option_group(HeatmapConstructionGroup)
	
	VirusGroupingGroup = optparse.OptionGroup(	parser		= parser,
							title		= "Virus grouping options",
							description	= "GRAViTy can estimate the distance cutoff that best separates the input reference taxonomic groupings, and report virus groups suggested by the estimated cutoff. "
									  "Theil's uncertainty correlation for the reference taxonomic grouping given the predicted grouping, and vice versa, are reported. "
									  "Symmetrical Theil's uncertainty correlation between the reference and predicted taxonomic grouping are also reported.")
	VirusGroupingGroup.add_option(	'--VirusGrouping',
					dest	= "VirusGrouping",
					default	= True,
					help	= "Perform virus grouping if True. [default: %default]",
					metavar	= "BOOLEAN",
					type	= "choice",
					choices	= ["True", "False",],
	)
	parser.add_option_group(VirusGroupingGroup)

	return parser