class Log_Generator_Pl1:
    '''Generate list of lists containing input config, for printing and saving to persistent log'''

    def __init__(self, options, fpath) -> None:
        self.options = options
        self.fpath = fpath

    def text_gen_start(self) -> list:
        return ["GRAViTy V2 adaptation by Mayne, R., Aiewsakun, P., Simmonds., P. et al. (2022)",
                "Based on original GRAViTy software https://github.com/PAiewsakun/GRAViTy",
                "Aiewsakun, P., Simmonds, P. The genomic underpinnings of eukaryotic virus taxonomy: creating a sequence-based framework for family-level virus classification. Microbiome 6, 38 (2018). https://doi.org/10.1186/s40168-018-0422-7"]

    def text_gen_pre_read(self) -> list:
        return ["Input for ReadGenomeDescTable:",
                "="*100,
                "Main input",
                "-"*50,
                f"GenomeDescTableFile: {self.options['GenomeDescTableFile']}",
                f"ShelveDir: {self.options['ShelveDir']}",
                f"Database: {self.options['Database']}",
                f"Database_Header: {self.options['Database_Header']}",
                f"TaxoGrouping_Header: {self.options['TaxoGrouping_Header']}",
                f"TaxoGroupingFile: {self.options['TaxoGroupingFile']}",
                "="*100
                ]

    def text_gen_pre_pphmdb(self) -> list:
        return ["Input for PPHMMDBConstruction:",
                "="*100,
                "Main input",
                "-"*50,
                f"GenomeSeqFile: {self.options['GenomeSeqFile']}",
                f"ShelveDir: {self.options['ShelveDir']}",
                "Protein extraction options",
                "-"*50,
                f"ProteinLength_Cutoff: {self.options['ProteinLength_Cutoff']}",
                f"IncludeProteinsFromIncompleteGenomes: {self.options['IncludeProteinsFromIncompleteGenomes']}",
                f"Protein clustering options",
                "-"*50,
                f"BLASTp_evalue_Cutoff: {self.options['BLASTp_evalue_Cutoff']}",
                f"BLASTp_PercentageIden_Cutoff: {self.options['BLASTp_PercentageIden_Cutoff']}",
                f"BLASTp_QueryCoverage_Cutoff: {self.options['BLASTp_QueryCoverage_Cutoff']}",
                f"BLASTp_SubjectCoverage_Cutoff: {self.options['BLASTp_SubjectCoverage_Cutoff']}",
                f"BLASTp_num_alignments: {self.options['BLASTp_num_alignments']}",
                f"N_CPUs: {self.options['N_CPUs']}",
                f"MUSCLE_GapOpenCost: {self.options['MUSCLE_GapOpenCost']}",
                f"MUSCLE_GapExtendCost: {self.options['MUSCLE_GapExtendCost']}",
                f"ProtClustering_MCLInflation: {self.options['ProtClustering_MCLInflation']}",
                f"Protein alignment merging options",
                "-"*50,
                f"N_AlignmentMerging: {self.options['N_AlignmentMerging']}",
                f"HHsuite_evalue_Cutoff: {self.options['HHsuite_evalue_Cutoff']}",
                f"HHsuite_pvalue_Cutoff:{self.options['HHsuite_pvalue_Cutoff']}",
                f"HHsuite_QueryCoverage_Cutoff: {self.options['HHsuite_QueryCoverage_Cutoff']}",
                f"HHsuite_SubjectCoverage_Cutoff: {self.options['HHsuite_SubjectCoverage_Cutoff']}",
                f"PPHMMClustering_MCLInflation_ForAlnMerging: {self.options['PPHMMClustering_MCLInflation_ForAlnMerging']}",
                f"HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging: {self.options['HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging']}",
                "="*100
                ]

    def text_gen_pre_refvirusannotator(self) -> list:
        return ["Input for RefVirusAnnotator:",
                "="*100,
                "Main input",
                "-"*50,
                f"GenomeSeqFile: {self.options['GenomeSeqFile']}",
                f"ShelveDir: {self.options['ShelveDir']}",
                "Reference virus annotation options",
                "-"*50,
                f"AnnotateIncompleteGenomes: {self.options['AnnotateIncompleteGenomes']}",
                f"N_CPUs: {self.options['N_CPUs']}",
                f"HMMER_C_EValue_Cutoff: {self.options['HMMER_C_EValue_Cutoff']}",
                f"HMMER_HitScore_Cutoff: {self.options['HMMER_HitScore_Cutoff']}",
                "Remove singleton PPHMM options",
                "-"*50,
                f"RemoveSingletonPPHMMs: {self.options['RemoveSingletonPPHMMs']}",
                f"N_VirusesOfTheClassToIgnore: {self.options['N_VirusesOfTheClassToIgnore']}",
                "PPHMM sorting options",
                "-"*50,
                f"PPHMMSorting: {self.options['PPHMMSorting']}",
                f"HHsuite_evalue_Cutoff: {self.options['HHsuite_evalue_Cutoff']}",
                f"HHsuite_pvalue_Cutoff: {self.options['HHsuite_pvalue_Cutoff']}",
                f"HHsuite_QueryCoverage_Cutoff: {self.options['HHsuite_QueryCoverage_Cutoff']}",
                f"HHsuite_SubjectCoverage_Cutoff: {self.options['HHsuite_SubjectCoverage_Cutoff']}",
                f"PPHMMClustering_MCLInflation_ForPPHMMSorting: {self.options['PPHMMClustering_MCLInflation_ForPPHMMSorting']}",
                "="*100
                ]

    def text_gen_pre_dendrogram(self) -> list:
        return ["Input for GRAViTyDendrogramAndHeatmapConstruction:",
                "="*100,
                "Main input",
                "-"*50,
                f"ShelveDir: {self.options['ShelveDir']}",
                f"AnnotateIncompleteGenomes: {self.options['AnnotateIncompleteGenomes']}",
                "Virus (dis)similarity measurement options",
                "-"*50,
                f"SimilarityMeasurementScheme: {self.options['SimilarityMeasurementScheme']}",
                f"p: {self.options['p']}",
                "Dendrogram construction options",
                "-"*50,
                f"Dendrogram: {self.options['Dendrogram']}",
                f"Dendrogram_LinkageMethod: {self.options['Dendrogram_LinkageMethod']}",
                f"Bootstrap: {self.options['Bootstrap']}",
                f"N_Bootstrap: {self.options['N_Bootstrap']}",
                f"Bootstrap_method: {self.options['Bootstrap_method']}",
                f"N_CPUs: {self.options['N_CPUs']}",
                "Heatmap construction options",
                "-"*50,
                f"Heatmap: {self.options['Heatmap']}",
                f"Heatmap_VirusOrderScheme: {self.options['Heatmap_VirusOrderScheme']}",
                f"Heatmap_WithDendrogram: {self.options['Heatmap_WithDendrogram']}",
                f"Heatmap_DendrogramFile: {self.options['Heatmap_DendrogramFile']}",
                f"Heatmap_DendrogramSupport_Cutoff: {self.options['Heatmap_DendrogramSupport_Cutoff']}",
                "Virus grouping options",
                "-"*50,
                f"VirusGrouping: {self.options['VirusGrouping']}",
                "="*100
                ]

    def text_gen_pre_mutalinfo(self) -> list:
        return ["Input for MutualInformationCalculator:",
                "="*100,
                "Main input",
                "-"*50,
                f"ShelveDir: {self.options['ShelveDir']}",
                f"AnnotateIncompleteGenomes: {self.options['AnnotateIncompleteGenomes']}",
                "Virus grouping for mutual information calculation options",
                "-"*50,
                f"VirusGroupingFile: {self.options['VirusGroupingFile']}",
                "Virus sampling options",
                "-"*50,
                f"N_Sampling: {self.options['N_Sampling']}",
                f"SamplingStrategy: {self.options['SamplingStrategy']}",
                f"SampleSizePerGroup: {self.options['SampleSizePerGroup']}",
                "="*100
                ]

    def entrypoint(self) -> list:
        print_title()
        logs = [self.text_gen_start()]
        logs.append(self.text_gen_pre_read())
        logs.append(self.text_gen_pre_pphmdb())
        logs.append(self.text_gen_pre_refvirusannotator())
        logs.append(self.text_gen_pre_dendrogram())
        logs.append(self.text_gen_pre_mutalinfo())

        with open(self.fpath + "/input_parameter_log.txt", "w") as f:
            for list_type in logs:
                [f.write(f"{item}\n") for item in list_type]

        return logs


class Log_Generator_Pl2:
    '''Generate list of lists containing input config, for printing and saving to persistent log'''

    def __init__(self, options, fpath) -> None:
        self.options = options
        self.fpath = fpath

    def text_gen_start(self) -> list:
        return ["GRAViTy V2 adaptation by Mayne, R., Aiewsakun, P., Simmonds., P. et al. (2022)",
                "Based on original GRAViTy software https://github.com/PAiewsakun/GRAViTy",
                "Aiewsakun, P., Simmonds, P. The genomic underpinnings of eukaryotic virus taxonomy: creating a sequence-based framework for family-level virus classification. Microbiome 6, 38 (2018). https://doi.org/10.1186/s40168-018-0422-7"]

    def text_gen_pre_read(self) -> list:
        return ["Input for ReadGenomeDescTable:",
                "="*100,
                "Main input",
                "-"*50,
                f"GenomeDescTableFile_UcfVirus: {self.options['GenomeDescTableFile_UcfVirus']}",
                f"ShelveDir_UcfVirus: {self.options['ShelveDir_UcfVirus']}",
                f"Database: {self.options['Database']}",
                f"Database_Header: {self.options['Database_Header']}",
                "="*100
                ]

    def text_gen_pre_pphmdb(self) -> list:
        return ["Input for PPHMMDBConstruction:",
                "="*100,
                "Main input",
                "-"*50,
                f"GenomeSeqFile_UcfVirus: {self.options['GenomeSeqFile_UcfVirus']}",
                f"ShelveDir_UcfVirus: {self.options['ShelveDir_UcfVirus']}",
                "Protein extraction options",
                "-"*50,
                f"ProteinLength_Cutoff: {self.options['ProteinLength_Cutoff']}",
                f"IncludeProteinsFromIncompleteGenomes: {self.options['IncludeProteinsFromIncompleteGenomes']}",
                f"Protein clustering options",
                "-"*50,
                f"BLASTp_evalue_Cutoff: {self.options['BLASTp_evalue_Cutoff']}",
                f"BLASTp_PercentageIden_Cutoff: {self.options['BLASTp_PercentageIden_Cutoff']}",
                f"BLASTp_QueryCoverage_Cutoff: {self.options['BLASTp_QueryCoverage_Cutoff']}",
                f"BLASTp_SubjectCoverage_Cutoff: {self.options['BLASTp_SubjectCoverage_Cutoff']}",
                f"BLASTp_num_alignments: {self.options['BLASTp_num_alignments']}",
                f"N_CPUs: {self.options['N_CPUs']}",
                f"MUSCLE_GapOpenCost: {self.options['MUSCLE_GapOpenCost']}",
                f"MUSCLE_GapExtendCost: {self.options['MUSCLE_GapExtendCost']}",
                f"ProtClustering_MCLInflation: {self.options['ProtClustering_MCLInflation']}",
                "Protein alignment merging options",
                "-"*50,
                f"N_AlignmentMerging: {self.options['N_AlignmentMerging']}",
                f"HHsuite_evalue_Cutoff: {self.options['HHsuite_evalue_Cutoff']}",
                f"HHsuite_pvalue_Cutoff: {self.options['HHsuite_pvalue_Cutoff']}",
                f"HHsuite_QueryCoverage_Cutoff: {self.options['HHsuite_QueryCoverage_Cutoff']}",
                f"HHsuite_SubjectCoverage_Cutoff: {self.options['HHsuite_SubjectCoverage_Cutoff']}",
                f"PPHMMClustering_MCLInflation_ForAlnMerging: {self.options['PPHMMClustering_MCLInflation_ForAlnMerging']}",
                f"HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging: {self.options['HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging']}",
                "="*100
                ]

    def text_gen_pre_ucfvirusannotator(self) -> list:
        return [f"Input for UcfVirusAnnotator:",
                "="*100,
                "Main input",
                "-"*50,
                f"GenomeSeqFile_UcfVirus: {self.options['GenomeSeqFile_UcfVirus']}",
                f"ShelveDir_UcfVirus: {self.options['ShelveDir_UcfVirus']}",
                f"ShelveDirs_RefVirus: {self.options['ShelveDirs_RefVirus']}",
                "Unclassified virus annotation options",
                "-"*50,
                f"AnnotateIncompleteGenomes_UcfVirus: {self.options['AnnotateIncompleteGenomes_UcfVirus']}",
                f"UsingDatabaseIncludingIncompleteRefViruses: {self.options['UsingDatabaseIncludingIncompleteRefViruses']}",
                f"N_CPUs: {self.options['N_CPUs']}",
                f"HMMER_C_EValue_Cutoff: {self.options['HMMER_C_EValue_Cutoff']}",
                f"HMMER_HitScore_Cutoff: {self.options['HMMER_HitScore_Cutoff']}",
                "="*100
                ]

    def text_gen_pre_classifier(self) -> list:
        return ["Input for VirusClassificationAndEvaluation:",
                "="*100,
                "Main input",
                "-"*50,
                f"ShelveDir_UcfVirus: {self.options['ShelveDir_UcfVirus']}",
                f"ShelveDirs_RefVirus: {self.options['ShelveDirs_RefVirus']}",
                f"AnnotateIncompleteGenomes_UcfVirus: {self.options['AnnotateIncompleteGenomes_UcfVirus']}",
                f"UsingDatabaseIncludingIncompleteRefViruses: {self.options['UsingDatabaseIncludingIncompleteRefViruses']}",
                f"Virus annotation (with PPHMM database derived from unclassified viruses) options",
                "-"*50,
                f"UseUcfVirusPPHMMs: {self.options['UseUcfVirusPPHMMs']}",
                f"GenomeSeqFile_UcfVirus: {self.options['GenomeSeqFile_UcfVirus']}",
                f"GenomeSeqFiles_RefVirus: {self.options['GenomeSeqFiles_RefVirus']}",
                f"N_CPUs: {self.options['N_CPUs']}",
                f"HMMER_C_EValue_Cutoff: {self.options['HMMER_C_EValue_Cutoff']}",
                f"HMMER_HitScore_Cutoff: {self.options['HMMER_HitScore_Cutoff']}",
                "Virus (dis)similarity measurement options",
                "-"*50,
                f"SimilarityMeasurementScheme: {self.options['SimilarityMeasurementScheme']}",
                f"p: {self.options['p']}",
                "Virus classification options",
                "-"*50,
                f"Dendrogram_LinkageMethod: {self.options['Dendrogram_LinkageMethod']}",
                f"Bootstrap: {self.options['Bootstrap']}",
                f"N_Bootstrap: {self.options['N_Bootstrap']}",
                f"Bootstrap_method: {self.options['Bootstrap_method']}",
                f"DatabaseAssignmentSimilarityScore_Cutoff: {self.options['DatabaseAssignmentSimilarityScore_Cutoff']}",
                f"N_PairwiseSimilarityScores: {self.options['N_PairwiseSimilarityScores']}",
                "Heatmap construction options",
                "-"*50,
                f"Heatmap_WithDendrogram: {self.options['Heatmap_WithDendrogram']}",
                f"Heatmap_DendrogramSupport_Cutoff: {self.options['Heatmap_DendrogramSupport_Cutoff']}",
                "Virus grouping options",
                "-"*50,
                f"VirusGrouping: {self.options['VirusGrouping']}",
                "="*100,
                "&"*100
                ]

    def entrypoint(self) -> list:
        print_title()
        logs = [self.text_gen_start()]
        logs.append(self.text_gen_pre_read())
        logs.append(self.text_gen_pre_pphmdb())
        logs.append(self.text_gen_pre_ucfvirusannotator())
        logs.append(self.text_gen_pre_classifier())

        with open(self.fpath + "/input_parameter_log.txt", "w") as f:
            for list_type in logs:
                [f.write(f"{item}\n") for item in list_type]

        return logs


def print_title():
    print('''
   _____ _____       __      ___ _______     __      _____
  / ____|  __ \     /\ \    / (_)__   __|    \ \    / /__ \
 | |  __| |__) |   /  \ \  / / _   | |_   _   \ \  / /   ) |
 | | |_ |  _  /   / /\ \ \/ / | |  | | | | |   \ \/ /   / /
 | |__| | | \ \  / ____ \  /  | |  | | |_| |    \  /   / /_
  \_____|_|  \_\/_/    \_\/   |_|  |_|\__, |     \/   |____|
                                       __/ |
                                      |___/
									  		''')
