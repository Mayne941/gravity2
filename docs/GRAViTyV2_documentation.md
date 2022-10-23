# GRAViTyV2 Documentation
Contents
1. Function descriptions
1. Example workflow

## Function Descriptions
### Pipeline I
Pipeline I (PL1) is for analysing and creating "databases" of information pertaining to these analyses, on "reference" viruses, i.e. genomes for which we have previously generated taxonomic data. PL1 comprises the following stages:
1. Read VMR data
1. Construct protein-protein hidden markov models (PPHMMs)
1. Annotate reference viruses
1. Do mutual information calculation functions and generate visualisations (optional)
The GRAViTyV2 API allows for running the entire pipeline, or entering from any specific point within it.

#### Minimum requirements to start PL1:
1. A table of reference virus descriptions ("GenomeDescTableFile"), in .csv format (.txt may also work but is not recommended). This table must either be an ICTV Virus Metadata Resource (VMR), or a document with a similar strucutre. As a minimum, your table will need the following columns: "Species", "Virus GENBANK accession", "Genome coverage". VMR scrape and VMR-like document generation from input FASTA data are included, see section "VMR Utilities".
1. A valid email address ("genbank_email") to download genomic data from GenBank, or otherwise genomic data in a .gb file. N.b. you must specify a file path for the GenBank file ("GenomeSeqFile"): if it doesn't exist, GRAViTyV2 will automatically attempt to download genomic data; if the file does exist, it will attempt to read data in from it.
1. A directory for saving all of your experimental data ("ShelveDir").

#### Options also exist for:
1. Specifying a database within your reference virus table to analyse, ignoring others, e.g. Baltimore Group ("Database", "Database_header").
1. Using a pre-defined grouping ("TaxoGrouping_Header", "TaxoGrouping_File")
1. Fine-tuning of all numerical procedures within.
1. Generating visualisations.

#### Output data:
1. Details of input parameters, to allow for full experiment reproducibility are stored to ./{experiment directory}
1. Intermediate data (for allowing re-runs or starting pipeline mid-way through) are saved to ./{experiment directory}/{BLAST}{HHsuite}}{HMMER}
1. Experiment results are saved to ./{experiment directory}/Shelves. This includes the mutual information score, dendrograms and visualisations. Pickle files are also saved here, which contain raw, serialised output, such that all experimental data can be retrieved if necessary.

### Pipeline II
Pipeline II (Pl2) is for comparing unclassified virus genomes against databases compiled via PL1. It comprises these stages:
1. Read in unclassified virus table
1. Compile protein-protein hidden markov models (PPHMMs, optional)
1. Annotate unclassified viruses
1. Prepare statistics and visualisations (optional)

#### Minimum requirements to start PL1:
1. Valid email address ("genbank_email") for downloading GenBank data OR .gb file ("GenomeSeqFile_UcfVirus"). A function for converting FASTA data to gb and generating a VMR-like table is included, see section "VMR Utilities". N.b. you must specify a file path for the GenBank file ("GenomeSeqFile"): if it doesn't exist, GRAViTyV2 will automatically attempt to download genomic data; if the file does exist, it will attempt to read data in from it.
1. A GenBank file from your PL1 run ("GenomeSeqFile_RefVirus")
1. A directory for PL1 output ("ShelveDir_RefVirus")
1. A directory to save output to ("ShelveDir_UcfVirus")
1. VMR or VMR-like documents for your PL1 run ("GenomeDescTableFile") and to correspond to your unclassified genomes file ("GenomeDescTableFile_UcfVirus").

#### Options also exist for:
1. Specifying a database within your reference virus table to analyse, ignoring others, e.g. Baltimore Group ("Database", "Database_header").
1. Using a pre-defined grouping ("TaxoGrouping_Header", "TaxoGrouping_File")
1. Fine-tuning of all numerical procedures within.
1. Generating visualisations.

#### Output data:
1. Details of input parameters, to allow for full experiment reproducibility are stored to ./{experiment directory}
1. Intermediate data (for allowing re-runs or starting pipeline mid-way through) are saved to ./{experiment directory}/{BLAST}{HHsuite}}{HMMER}
1. Experiment results are saved to ./{experiment directory}/Shelves. This includes the classification results, dendrograms and visualisations. Pickle files are also saved here, which contain raw, serialised output, such that all experimental data can be retrieved if necessary.

## VMR Utilities
### Scrape VMR
Fetch virus metadata resource (VMR) file from ICTV website, save as a csv file. Some basic cleaning is done at this stage, e.g. some new columns synthesised to extract Baltimore group, removal of incomplete genomes.

### Construct first pass VMR
Filter a VMR file by a specific rule, for significantly decreasing the length of a VMR to make a granular, "first pass" run. The purpose of this is to build a representative dataset with one example genome frome each taxon. Save results as csv.

Filter rules:
- Threshold: function will start at the Family level; if n samples within family < threshold, split at the genus level. If n samples within genus < threshold, extract species name. This ensures that both monotypic and polytypic genomes are adequately represented.
- Additional filter: optional, we recommend building first pass databases at RNA, dsDNA and ssDNA levels. Hence each unclassified virus should be run against each of these three first pass databases to get an approximate classification.

### Construct second pass VMR
For use after the first pass: classification results from first pass PL2 will indicate which existing taxa are most similar to unknown samples. This function may then be used to filter by a specific label, to construct a VMR with all individual species from this taxon. Output as CSV.

- Filter level: column name in VMR, e.g. Family
- Filter name: specific grouping's name, e.g. "Caulimoviridae"

Depending on the complexity of your second pass, you may need to run this function several times and combine their output, e.g. if you are comparing multiple, distantly related unclassified viruses and want to compare them to several taxa in one run.

### Convert FASTA to genbank and VMR
Input a FASTA file, output a corresponding .gb (GenBank) file and associated VMR-like .csv file.

## Example Workflow: Get taxonomy for single unknown RNA virus, from FASTA input
All necessary data to run this example workflow are found in ./docs/example. N.b. parameters can't be copy-pasted exactly: ensure you check all parameters, including genbank_email and n_cpus.
1. Get latest VMR with scrape_vmr endpoint
```
{
  "save_path": "./data/",
  "vmr_name": "latest_vmr.csv"
}
```
2. Construct a first pass VMR with construct_first_pass_vmr endpoint
```
{
  "save_path": "./data/",
  "vmr_name": "latest_vmr.csv",
  "save_name": "latest_vmr_rna_filter.csv",
  "filter_threshold": 5,
  "additional_filter": "RNA"
}
```
4. Create a first pass PL1 database for filtered RNA viruses with pipeline_i_full endpoint. This took about 15 minutes to run on our dev machine.
```
{
  "genbank_email": "name@provider.com",
  "GenomeDescTableFile": "./data/latest_vmr_rna_filter.csv",
  "ShelveDir": "./output/pipeline_1_rna",
  "Database": null,
  "Database_Header": null,
  "TaxoGrouping_Header": "Taxonomic grouping",
  "TaxoGroupingFile": null,
  "GenomeSeqFile": "./output/pipeline_1_rna.gb",
  "ProteinLength_Cutoff": 100,
  "IncludeProteinsFromIncompleteGenomes": true,
  "BLASTp_evalue_Cutoff": 0.001,
  "BLASTp_PercentageIden_Cutoff": 50,
  "BLASTp_QueryCoverage_Cutoff": 75,
  "BLASTp_SubjectCoverage_Cutoff": 75,
  "BLASTp_num_alignments": 1000000,
  "N_CPUs": 24,
  "MUSCLE_GapOpenCost": -3,
  "MUSCLE_GapExtendCost": 0,
  "ProtClustering_MCLInflation": 2,
  "N_AlignmentMerging": 0,
  "PPHMMClustering_MCLInflation_ForAlnMerging": 5,
  "HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging": true,
  "HHsuite_evalue_Cutoff": 0.000001,
  "HHsuite_pvalue_Cutoff": 0.05,
  "HHsuite_QueryCoverage_Cutoff": 85,
  "HHsuite_SubjectCoverage_Cutoff": 85,
  "AnnotateIncompleteGenomes": false,
  "HMMER_C_EValue_Cutoff": 0.001,
  "HMMER_HitScore_Cutoff": 0,
  "RemoveSingletonPPHMMs": false,
  "N_VirusesOfTheClassToIgnore": 1,
  "PPHMMSorting": false,
  "PPHMMClustering_MCLInflation_ForPPHMMSorting": 2,
  "p": 1,
  "Dendrogram": true,
  "Dendrogram_LinkageMethod": "average",
  "Bootstrap": false,
  "N_Bootstrap": 10,
  "Bootstrap_method": "booster",
  "Heatmap": false,
  "Heatmap_VirusOrderScheme": null,
  "Heatmap_WithDendrogram": true,
  "Heatmap_DendrogramFile": null,
  "Heatmap_DendrogramSupport_Cutoff": 0.75,
  "VirusGrouping": true,
  "VirusGroupingFile": null,
  "N_Sampling": 10,
  "SamplingStrategy": "balance_with_repeat",
  "SampleSizePerGroup": 10,
  "SimilarityMeasurementScheme": "PG"
}
```
![First pass PL1 Heatmap](/docs/example/HeatmapWithDendrogram.IncompleteGenomes=0.Scheme=PG.Method=average.p=1.0.support_cutoff=0.75.png)
*First pass PL1 heatmap, showing representative examples from each Family in the RNA viruses where Family has at least 5 members. The dendrogram demonstrates that this is an extremely coarse representation. We will now compare our unknown virus to these reference viruses to get an approximate classification and allow for a more granular, second pass to be run.*

5. Construct genbank file and VMR-like table for unknown virus, using convert_fasta_to_genbank and a fasta file "my_genome".
```
{
  "save_path": "./docs/example",
  "fasta_fname": "my_fasta",
  "genbank_fname": "my_genbank.gb",
  "vmr_fname": "my_vmr.csv"
}
```
6. Compare unknown virus vs first pass database with pipeline_ii_full endpoint. Several optional parameters are disabled which allows for an extremely brief run time (8 s on our dev machine).
```
{
  "genbank_email": "name@provider.com",
  "GenomeDescTableFile": "data/latest_vmr_rna_filter.csv",
  "ShelveDir_UcfVirus": "output/my_experiment_firstpass",
  "ShelveDirs_RefVirus": "output/pipeline_1_rna",
  "GenomeDescTableFile_UcfVirus": "docs/example/my_vmr.csv",
  "Database": null,
  "Database_Header": null,
  "GenomeSeqFile_UcfVirus": "docs/example/my_genbank.gb",
  "GenomeSeqFiles_RefVirus": "output/pipeline_1_rna.gb",
  "UseUcfVirusPPHMMs": false,
  "ProteinLength_Cutoff": 100,
  "IncludeProteinsFromIncompleteGenomes": true,
  "BLASTp_evalue_Cutoff": 0.001,
  "BLASTp_PercentageIden_Cutoff": 50,
  "BLASTp_QueryCoverage_Cutoff": 75,
  "BLASTp_SubjectCoverage_Cutoff": 75,
  "BLASTp_num_alignments": 1000000,
  "N_CPUs": 24,
  "MUSCLE_GapOpenCost": -3,
  "MUSCLE_GapExtendCost": 0,
  "ProtClustering_MCLInflation": 2,
  "N_AlignmentMerging": 0,
  "PPHMMClustering_MCLInflation_ForAlnMerging": 5,
  "HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging": true,
  "HHsuite_evalue_Cutoff": 0.000001,
  "HHsuite_pvalue_Cutoff": 0.05,
  "HHsuite_QueryCoverage_Cutoff": 85,
  "HHsuite_SubjectCoverage_Cutoff": 85,
  "AnnotateIncompleteGenomes_UcfVirus": true,
  "UsingDatabaseIncludingIncompleteRefViruses": false,
  "HMMER_C_EValue_Cutoff": 0.001,
  "HMMER_HitScore_Cutoff": 0,
  "p": 1,
  "Dendrogram": true,
  "Dendrogram_LinkageMethod": "average",
  "Bootstrap": false,
  "N_Bootstrap": 10,
  "Bootstrap_method": "booster",
  "DatabaseAssignmentSimilarityScore_Cutoff": 0.01,
  "N_PairwiseSimilarityScores": 10000,
  "Heatmap_WithDendrogram": true,
  "Heatmap_DendrogramSupport_Cutoff": 0.75,
  "VirusGrouping": true,
  "SimilarityMeasurementScheme": "PG"
}
```
![First pass PL2](/docs/example/HeatmapWithDendrogram.RefVirusGroup=pipeline_1_rna.IncompleteUcfRefGenomes=10.Scheme=PG.Method=average.p=1.0.png)
*First pass PL2 output indicates that our unknown virus has high homology with Amalgaviridae. Numerical results are stored in ./{experiment directory}/Shelves/ClassificationResults.csv. This will be used to inform the design of our more granular second pass.*

7. Construct a second pass VMR with the construct_second_pass_vmr endpoint.
```
{
  "save_path": "./data/",
  "vmr_name": "latest_vmr.csv",
  "save_name": "latest_vmr_secondpass_filter.csv",
  "filter_level": "Family",
  "filter_name": "Amalgaviridae"
}
```

8. Run PL1 in full, on second pass VMR, to build fine grained database of data on specific family (Amalgaviridae). Optional parameters are enabled to increase accuracy.
```
{
  "genbank_email": "name@provider.com",
  "GenomeDescTableFile": "./data/latest_vmr_secondpass_filter.csv",
  "ShelveDir": "./output/pipeline_2_amalgaviridae",
  "Database": null,
  "Database_Header": null,
  "TaxoGrouping_Header": "Taxonomic grouping",
  "TaxoGroupingFile": null,
  "GenomeSeqFile": "./output/pipeline_2_amalgaviridae.gb",
  "ProteinLength_Cutoff": 100,
  "IncludeProteinsFromIncompleteGenomes": true,
  "BLASTp_evalue_Cutoff": 0.001,
  "BLASTp_PercentageIden_Cutoff": 50,
  "BLASTp_QueryCoverage_Cutoff": 75,
  "BLASTp_SubjectCoverage_Cutoff": 75,
  "BLASTp_num_alignments": 1000000,
  "N_CPUs": 24,
  "MUSCLE_GapOpenCost": -3,
  "MUSCLE_GapExtendCost": 0,
  "ProtClustering_MCLInflation": 2,
  "N_AlignmentMerging": 0,
  "PPHMMClustering_MCLInflation_ForAlnMerging": 5,
  "HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging": true,
  "HHsuite_evalue_Cutoff": 0.000001,
  "HHsuite_pvalue_Cutoff": 0.05,
  "HHsuite_QueryCoverage_Cutoff": 85,
  "HHsuite_SubjectCoverage_Cutoff": 85,
  "AnnotateIncompleteGenomes": false,
  "HMMER_C_EValue_Cutoff": 0.001,
  "HMMER_HitScore_Cutoff": 0,
  "RemoveSingletonPPHMMs": false,
  "N_VirusesOfTheClassToIgnore": 1,
  "PPHMMSorting": false,
  "PPHMMClustering_MCLInflation_ForPPHMMSorting": 2,
  "p": 1,
  "Dendrogram": true,
  "Dendrogram_LinkageMethod": "average",
  "Bootstrap": true,
  "N_Bootstrap": 10,
  "Bootstrap_method": "booster",
  "Heatmap": false,
  "Heatmap_VirusOrderScheme": null,
  "Heatmap_WithDendrogram": true,
  "Heatmap_DendrogramFile": null,
  "Heatmap_DendrogramSupport_Cutoff": 0.75,
  "VirusGrouping": true,
  "VirusGroupingFile": null,
  "N_Sampling": 10,
  "SamplingStrategy": "balance_with_repeat",
  "SampleSizePerGroup": 10,
  "SimilarityMeasurementScheme": "PG"
}
```
![Second pass PL1](/docs/example/HeatmapWithDendrogram.IncompleteGenomes=0.Scheme=PG.Method=average.p=1.0.support_cutoff=0.75_pl2.png)
*Second pass PL1 output, showing all members of Amalgaviridae family*

9. Run PL2 to compare unknown virus with second pass PL1 database. Optional parameters are enabled to enhance accuracy.
```
{
  "genbank_email": "name@provider.com",
  "GenomeDescTableFile": "data/latest_vmr_secondpass_filter.csv",
  "ShelveDir_UcfVirus": "output/my_experiment_secondpass",
  "ShelveDirs_RefVirus": "output/pipeline_2_amalgaviridae",
  "GenomeDescTableFile_UcfVirus": "docs/example/my_vmr.csv",
  "Database": null,
  "Database_Header": null,
  "GenomeSeqFile_UcfVirus": "docs/example/my_genbank.gb",
  "GenomeSeqFiles_RefVirus": "output/pipeline_2_amalgaviridae.gb",
  "UseUcfVirusPPHMMs": true,
  "ProteinLength_Cutoff": 100,
  "IncludeProteinsFromIncompleteGenomes": true,
  "BLASTp_evalue_Cutoff": 0.001,
  "BLASTp_PercentageIden_Cutoff": 50,
  "BLASTp_QueryCoverage_Cutoff": 75,
  "BLASTp_SubjectCoverage_Cutoff": 75,
  "BLASTp_num_alignments": 1000000,
  "N_CPUs": 24,
  "MUSCLE_GapOpenCost": -3,
  "MUSCLE_GapExtendCost": 0,
  "ProtClustering_MCLInflation": 2,
  "N_AlignmentMerging": 0,
  "PPHMMClustering_MCLInflation_ForAlnMerging": 5,
  "HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging": true,
  "HHsuite_evalue_Cutoff": 0.000001,
  "HHsuite_pvalue_Cutoff": 0.05,
  "HHsuite_QueryCoverage_Cutoff": 85,
  "HHsuite_SubjectCoverage_Cutoff": 85,
  "AnnotateIncompleteGenomes_UcfVirus": true,
  "UsingDatabaseIncludingIncompleteRefViruses": false,
  "HMMER_C_EValue_Cutoff": 0.001,
  "HMMER_HitScore_Cutoff": 0,
  "p": 1,
  "Dendrogram": true,
  "Dendrogram_LinkageMethod": "average",
  "Bootstrap": true,
  "N_Bootstrap": 10,
  "Bootstrap_method": "booster",
  "DatabaseAssignmentSimilarityScore_Cutoff": 0.01,
  "N_PairwiseSimilarityScores": 10000,
  "Heatmap_WithDendrogram": true,
  "Heatmap_DendrogramSupport_Cutoff": 0.75,
  "VirusGrouping": true,
  "SimilarityMeasurementScheme": "PG"
}
```
![Second pass PL2](/docs/example/HeatmapWithDendrogram.RefVirusGroup=pipeline_2_amalgaviridae.IncompleteUcfRefGenomes=10.Scheme=PG.Method=average.p=1.0_pl2.png)
*Output demonstrates how unclassified virus relates to existing reference viruses. Numerical results are shown in ./{experiment directory}/Shelves/ClassificationResults.{csv}{txt}*
