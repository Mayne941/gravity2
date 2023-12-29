from fastapi import Query
from typing import Union, Literal
from pydantic import BaseModel, Field, FilePath, DirectoryPath


'''PRIMITIVES'''
'''Utility fns'''
class Data_vmr_name(BaseModel):
    vmr_name: str = Query('latest_vmr.csv',
                          description="Filename for new VMR.")

class Data_save_path(BaseModel):
    save_path: DirectoryPath = Query('./data/',
                                     description="Path to save the Virus Metadata Resource (VMR).")

class Data_vmr_filter_params(BaseModel):
    save_name: str = Query('latest_vmr_first_pass_filter.csv',
                        description="Filename of output VMR file")
    filter_level: str = Query('Family',
                              description="Specify which taxo grouping to filter by. Used in combination with filter_name to generate a VMR of all viral genomes within a specific taxo grouping.")
    filter_name: str = Query('Caudoviricetes',
                             description="Used to specify which taxo grouping name to filter by: used in combination with filter_level to generate a VMR of all viral genomes within a specific taxo grouping.")

class Data_vmr_filter_threshold(BaseModel):
    filter_threshold: int = Query(10,
                                  description="If filter = true, how many members should be in each taxo grouping? If < threshold members in family, resolve at genus level, and likewise for species.")

class Data_fasta_fns(BaseModel):
    fasta_fname: str = Query('my_fasta',
                             description="Filename for fasta file to convert.")
    genbank_fname: str = Query('data/my_genbank.gb',
                               description="Filename for output genbank. ")
    vmr_fname: str = Query('data/my_vmr.csv',
                           description="Filename for VMR-like document to produce alongside genbank file.")

class Data_genome_segment_fpath(BaseModel):
    combined_seq_name: str = Query("my_seq.gb",
                                   description="Give a name to the virus whose genome segments you are concatenating")

'''E2E'''
class Data_e2e_unique_fns(BaseModel):
    ExperimentName: str = Query("myexp",
                                description="Exeriment name will prefix all output folders (and files not contained in the main experiment folders, such as genbank files)")
    SkipFirstPass: bool = Query(False,
                                        description="Opt to skip first pass through both pipelines. N.b. first pass parameters still need to be filled in with valid entries.")
    GenomeDescTableFile_FirstPass: FilePath = Query('./data/latest_vmr_firstpass.csv',
                                          description="FIRST PASS PIPELINE 1: Full path to the Virus Metadata Resource (VMR) tab delimited file, wth headers. VMR can be downloaded using the scrape endpoint.")
    GenomeDescTableFile_SecondPass: FilePath = Query('./data/latest_vmr_secondpass.csv',
                                          description="SECOND PASS PIPELINE 1: Full path to the Virus Metadata Resource (VMR) tab delimited file, wth headers. VMR can be downloaded using the scrape endpoint.")
    TaxoGroupingFile_FirstPass: Union[FilePath, None] = Query(None,
                                                    description="FIRST PASS PIPELINE 1: It is possible that the user might want to associate different viruses with different taxonomic assignment levels, e.g. family assignments for some viruses, and subfamily or genus assignments for some other viruses, etc. To accomodate this need, the user can either add a column in the VMR file, and use --TaxoGrouping_Header to specify the column (see --TaxoGrouping_Header). Alternatively, the user can provide a file (with no header) that contains a single column of taxonomic groupings for all viruses in the order that appears in the VMR file. The user can specify the full path to the taxonomic grouping file using this options. If this option is used, it will override the one specified by --TaxoGrouping_Header.")
    TaxoGroupingFile_SecondPass: Union[FilePath, None] = Query(None,
                                                    description="SECOND PASS PIPELINE 1: It is possible that the user might want to associate different viruses with different taxonomic assignment levels, e.g. family assignments for some viruses, and subfamily or genus assignments for some other viruses, etc. To accomodate this need, the user can either add a column in the VMR file, and use --TaxoGrouping_Header to specify the column (see --TaxoGrouping_Header). Alternatively, the user can provide a file (with no header) that contains a single column of taxonomic groupings for all viruses in the order that appears in the VMR file. The user can specify the full path to the taxonomic grouping file using this options. If this option is used, it will override the one specified by --TaxoGrouping_Header.")


'''PL1'''
class Data_pl1_unique_params(BaseModel):
    TaxoGrouping_Header: Literal["Taxonomic grouping", "Family"] = Query('Taxonomic grouping',
                                                                         description="The header of the Taxonomic grouping column.")
    AnnotateIncompleteGenomes: bool = Query(False,
                                            description="Annotate all unclassified viruses using reference PPHMM database(s) if True, otherwise only complete genomes.")
    RemoveSingletonPPHMMs: bool = Query(False,
                                        description="Remove singleton PPHMMs from the database if True.")
    N_VirusesOfTheClassToIgnore: int = Field(1, gt=0,
                                             description="When 'RemoveSingletonPPHMMs' == TRUE, singleton PPHMMs are removed from the PPHMM database only if they show similarity to viruses that belong to taxonomic groups with more than NUMBER members.")
    PPHMMSorting: bool = Query(False,
                               description="Sort PPHMMs if True.")
    PPHMMClustering_MCLInflation_ForPPHMMSorting: int = Field(2, gt=0,
                                                              description="Cluster granularity. Increasing INFLATION will increase cluster granularity.")
    VirusGroupingFile: Union[FilePath, None] = Query(None, # RM < TODO Cull?
                                                     description="Fill path to the virus grouping scheme file. The file contains column(s) of arbitrary taxonomic grouping scheme(s) that users want to investigate. Note that file must contain headers, see docs for further info. If 'None', the taxonomic grouping as specified in 'Taxonomic grouping' column in the VMR will be used.")
    N_Sampling: int = Field(10, gt=0,
                            description="The number of mutual information scores sample size.")
    SamplingStrategy: Literal[None, "balance_without_repeat", "balance_with_repeat"] = Query('balance_with_repeat',
                                                                                             description="Virus sampling scheme.")
    SampleSizePerGroup: int = Field(10, gt=0,
                                    description="If 'SamplingStrategy' != None, this option specifies the number of viruses to be sampled per taxonomic group")


'''PL2'''
class Data_pl2_unique_params(BaseModel):
    GenomeDescTableFile_UcfVirus: FilePath = Query("data/unclassified_viruses_vmr.csv",
                                                   description="Full path to the Virus Metadata Resource-like (VMR-like) tab delimited file of unclassified viruses, wth headers. VVMR can be downloaded using the scrape endpoint")
    UseUcfVirusPPHMMs: bool = Query(True,
                                    description="Annotate reference and unclassified viruses using the PPHMM database derived from unclassified viruses if True.")
    AnnotateIncompleteGenomes_UcfVirus: bool = Query(True,
                                                     description="Annotate all unclassified viruses using reference PPHMM database(s) if True, otherwise only complete genomes.")
    UsingDatabaseIncludingIncompleteRefViruses: bool = Query(False,
                                                             description="Annotate unclassified viruses using the PPHMM and GOM databases derived from all reference viruses if True, otherwise using those derived from complete reference genomes only.")
    DatabaseAssignmentSimilarityScore_Cutoff: float = Field(0.01, gt=0,
                                                            description="Threshold to determine if the unclassified virus at least belongs to a particular database. For example, an unclassified virus is assigned to the family 'X' in the reference 'Baltimore group X' database, with the (greatest) similarity score of 0.1. This score might be too low to justify that the virus is a member of the family 'X', and fail the similarity threshold test. However, since the similarity score of 0.1 > %default, GRAViTy will make a guess that it might still be a virus of the 'Baltimore group X' database, under the default setting.")
    N_PairwiseSimilarityScores: int = Field(10000, gt=0,
                                            description="Number of data points in the distributions of intra- and inter-group similarity scores used to estimate the similarity threshold. ")

'''General pipeline'''
class Data_common_pipeline_params(BaseModel):
    genbank_email: str = Query('name@provider.com',
                               description="A valid email address is required to download genbank files.")
    Database: Union[str, None] = Query(None,  # RM < TODO Cull?
                                       description="GRAViTy will only analyse genomes that are labelled with DATABASE in the database column. The database column can be specified by the DATABASE HEADER argument. If 'None', all entries are analysed.")
    Database_Header: Union[str, None] = Query(None, # RM < TODO Cull?
                                              description="The header of the database column. Cannot be 'None' if DATABASE is specified.")
    ProteinLength_Cutoff: int = Field(100, gt=0,
                                      description="Proteins with length < LENGTH aa will be ignored")
    IncludeProteinsFromIncompleteGenomes: bool = Query(True,
                                                       description="Include protein sequences from incomplete genomes to the database if True.")
    Mash_p_val_cutoff: float = Field(0.05,
                                        description="P value threshold below which results in Mash analysis will be ignored.")
    Mash_sim_score_cutoff: float = Field(0.5,
                                              description="Similarity score threshold below which results in Mash analysis will be ignored.")
    ProtClustering_MCLInflation: int = Field(2, gt=0,
                                             description="Cluster granularity. Increasing INFLATION will increase cluster granularity.")
    N_AlignmentMerging: int = Query(0,    # RM < TODO Cull?
                                               description="Number of rounds of alignment merging. ROUND == 0 means no merging. ROUND == -1 means merging until exhausted, ROUND int = n rounds.")
    PPHMMClustering_MCLInflation_ForAlnMerging: int = Field(5, gt=0,
                                                            description="Cluster granularity. Increasing INFLATION will increase cluster granularity.")
    HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging: bool = Query(True,
                                                           description="Make a HMMER PPHMM DB for each round of protein merging if True.")
    HHsuite_evalue_Cutoff: float = Field(1e-06, gt=0,
                                         description="Threshold for PPHMM similarity detection. A hit with an E-value > E-VALUE will be ignored.")
    HHsuite_pvalue_Cutoff: float = Field(0.05, ge=0, le=1,
                                         description="Threshold for PPHMM similarity detection. A hit with a p-value > P-VALUE will be ignored.")
    HHsuite_QueryCoverage_Cutoff: float = Field(85.0, ge=0, le=100,
                                                description="Threshold for PPHMM similarity detection. A hit with a query coverage < COVERAGE will be ignored.")
    HHsuite_SubjectCoverage_Cutoff: float = Field(85.0, ge=0, le=100,
                                                  description="Threshold for PPHMM similarity detection. A hit with a subject coverage < COVERAGE will be ignored")
    HMMER_C_EValue_Cutoff: float = Field(0.001, gt=0,
                                         description="Threshold for HMM-protein similarity detection. A hit with an E-value > E-VALUE will be ignored.")
    HMMER_HitScore_Cutoff: int = Field(0, ge=0,
                                       description="Threshold for HMM-protein similarity detection. A hit with a score < SCORE will be ignored.")
    p: float = Field(1, ge=0,
                     description="Distance transformation P coefficient, 0 <= P. D = 1 - S**P. If P = 1, no distance transformation is applied. If P > 1, shallow branches will be stretched out so shallow splits can be seen more clearly. If p < 1, deep branches are stretch out so that deep splits can be seen more clearly. If p = 0, the entire dendrogram will be reduced to a single branch (root), with all taxa forming a polytomy clade at the tip of the dendrogram. If p = Inf, the resultant dendrogram will be star-like, with all taxa forming a polytomy clade at the root.")
    Dendrogram_LinkageMethod: Literal["single", "complete", "average", "weighted", "centroid", "median", "ward"] = Query('average',
                                                                                                                         description="LINKAGE for dendrogram construction. If LINKAGE = 'single', the nearest point algorithm is used to cluster viruses and compute cluster distances. If LINKAGE = 'complete', the farthest point algorithm is used to cluster viruses and compute cluster distances. If LINKAGE = 'average', the UPGMA algorithm is used to cluster viruses and compute cluster distances. If LINKAGE = 'weighted', the WPGMA algorithm is used to cluster viruses and compute cluster distances. If LINKAGE = 'centroid', the UPGMC algorithm is used to cluster viruses and compute cluster distances. If LINKAGE = 'median', the WPGMC algorithm is used to cluster viruses and compute cluster distances. If LINKAGE = 'ward', the incremental algorithm is used to cluster viruses and compute cluster distances.")
    Bootstrap: bool = Query(True,
                            description="Perform bootstrapping if True.")
    N_Bootstrap: int = Field(10, gt=0,
                             description="The number of pseudoreplicate datasets by resampling.")
    Bootstrap_method: Literal["booster", "sumtrees"] = Query('sumtrees',
                                                             description="Two METHODs for tree summary construction are implemented in GRAViTy. If METHOD = 'sumtrees', SumTrees (Sukumaran, J & MT Holder, 2010, Bioinformatics; https://dendropy.org/programs/sumtrees.html) will be used to summarize non-parameteric bootstrap support for splits on the best estimated dendrogram. The calculation is based on the standard Felsenstein bootstrap method. If METHOD = 'booster', BOOSTER (Lemoine et al., 2018, Nature; https://booster.pasteur.fr/) will be used. With large trees and moderate phylogenetic signal, BOOSTER tends to be more informative than the standard Felsenstein bootstrap method.")
    VirusGrouping: bool = Query(True,
                                description="Perform virus grouping if True.")
    SimilarityMeasurementScheme: Literal["P", "G", "L", "PG", "PL"] = Query("PG",
                                                                            description="Virus similarity measurement SCHEMEs. If SCHEME = 'P', an overall similarity between two viruses is their GJ_P. If SCHEME = 'L', an overall similarity between two viruses is their GJ_L. If SCHEME = 'G', an overall similarity between two viruses is their GJ_G. If SCHEME = 'PG', an overall similarity between two viruses is a geometric mean - or a 'composite generalised Jaccard score' (CGJ) - of their GJ_P and GJ_G. If SCHEME = 'PL', an overall similarity between two viruses is a geometric mean - or a 'composite generalised Jaccard score' (CGJ) - of their GJ_P and GJ_L.")
    Heatmap_WithDendrogram: bool = Query(True,
                                         description="Construct (dis)similarity heatmap with dendrogram if True.")
    Heatmap_DendrogramSupport_Cutoff: float = Field(0.75, ge=0, le=1,
                                                    description="Threshold for the BOOTSTRAP SUPPORT to be shown on the dendrogram on the heatmap.")


'''ENDPOINT MODEL COLLECTIONS'''
'''Utility fns'''
class Endp_data_scrape_data(Data_save_path, Data_vmr_name):
    pass

class Endp_data_first_pass_taxon_filter(Data_save_path, Data_vmr_name, Data_vmr_filter_params, Data_vmr_filter_threshold):
    pass

class Endp_data_second_pass_filter(Data_save_path, Data_vmr_name, Data_vmr_filter_params):
    pass

class Endp_data_fasta_to_gb(Data_fasta_fns):
    pass

class CombineGenomeSegs(Data_fasta_fns, Data_genome_segment_fpath):
    pass


'''Pipelines'''
class E2e_data(Data_common_pipeline_params, Data_pl1_unique_params, Data_pl2_unique_params, Data_e2e_unique_fns):
    '''PL1 (no unique PL2/General)'''
    pass

class Pipeline_i_data(Data_pl1_unique_params, Data_common_pipeline_params):
    GenomeDescTableFile: FilePath = Query('./data/latest_vmr.csv',
                                          description="Full path to the Virus Metadata Resource (VMR) tab delimited file, wth headers. VMR can be downloaded using the scrape endpoint.")
    ExpDir: str = Query('./output/myexperiment_pipeline_1',
                           description="Full path to the shelve directory, storing GRAViTy outputs. Makes new dir if not exists.")
    TaxoGroupingFile: Union[FilePath, None] = Query(None,
                                                    description="It is possible that the user might want to associate different viruses with different taxonomic assignment levels, e.g. family assignments for some viruses, and subfamily or genus assignments for some other viruses, etc. To accomodate this need, the user can either add a column in the VMR file, and use --TaxoGrouping_Header to specify the column (see --TaxoGrouping_Header). Alternatively, the user can provide a file (with no header) that contains a single column of taxonomic groupings for all viruses in the order that appears in the VMR file. The user can specify the full path to the taxonomic grouping file using this options. If this option is used, it will override the one specified by --TaxoGrouping_Header.")
    GenomeSeqFile: str = Query('./output/ref_sequences.gb',
                               description="Full path to the genome sequence GenBank file. If the file doesn't exist, GRAViTy will download the sequences from the NCBI database using accession numbers specified in the VMR file, 'Virus GENBANK accession' column")


class Pipeline_ii_data(Data_pl2_unique_params, Data_common_pipeline_params):
    GenomeDescTableFile: FilePath = Query('data/latest_vmr.csv',
                                          description="Full path to the Virus Metadata Resource (VMR) tab delimited file, wth headers. VMR can be downloaded using the scrape endpoint")
    ShelveDir_UcfVirus: str = Query('output/myexperiment_pipeline_2',
                                    description="Full path to the shelve directory of unclassified viruses, storing GRAViTy outputs.")
    ShelveDirs_RefVirus: str = Query('output/myexperiment_pipeline_1',
                                     description="Full path(s) to the shelve director(y/ies) of reference virus(es). For example: 'path/to/shelve/ref1, path/to/shelve/ref2, ...'")
    GenomeSeqFile_UcfVirus: str = Query('output/unclassified_seqs.gb',
                                        description="Full path to the genome sequence GenBank file of unclassified viruses.")
    GenomeSeqFiles_RefVirus: str = Query('output/ref_sequences.gb',
                                         description="Full path(s) to the genome sequence GenBank file(s) of reference viruses. For example: 'path/to/GenBank/ref1, path/to/GenBank/ref2, ...' This cannot be 'None' if UseUcfVirusPPHMMs = True. ")


'''Deprecated'''
# class FirstPassBaltimoreFilter(BaseModel):
#     save_path: DirectoryPath = Query('./data/',
#                                      description="Path to save the Virus Metadata Resource (VMR).")
#     vmr_name: str = Query('latest_vmr.csv',
#                           description="Filename for VMR to edit.")
#     save_name: str = Query('latest_vmr_first_pass_filter.csv',
#                            description="Filename of output VMR file")
#     filter_threshold: int = Query(10,
#                                   description="If filter = true, how many members should be in each taxo grouping? If < threshold members in family, resolve at genus level, and likewise for species.")
#     baltimore_filter: Literal["none", "dsDNA", "ssDNA", "RNA"] = Query("RNA",
#                                                                         description="Do an extra filter to further reduce size of first pass database. Currently supports: `none`, `dsDNA`, `ssDNA`, `RNA`")

# class VmrFilter(BaseModel, Data_save_path, Data_vmr_name):
#     save_path: DirectoryPath = Query('./data/',
#                                      description="Path to save the Virus Metadata Resource (VMR).")
#     vmr_name: str = Query('latest_vmr.csv',
#                           description="Filename for VMR to edit.")
#     save_name: str = Query('latest_vmr_filtered.csv',
#                            description="Filename of output VMR file")
#     filter_level: str = Query('Family',
#                               description="Specify which taxo grouping to filter by. Filter will remove all except first instance from this taxon.")
