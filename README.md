```
   _____ _____       __      ___ _______     __      _____
  / ____|  __ \     /\ \    / (_)__   __|    \ \    / /__ \
 | |  __| |__) |   /  \ \  / / _   | |_   _   \ \  / /   ) |
 | | |_ |  _  /   / /\ \ \/ / | |  | | | | |   \ \/ /   / /
 | |__| | | \ \  / ____ \  /  | |  | | |_| |    \  /   / /_
  \_____|_|  \_\/_/    \_\/   |_|  |_|\__, |     \/   |____|
                                       __/ |
                                      |___/
 ```
GRAViTy-V2 is a viral taxonomy application developed by Richard Mayne, Pakorn Aiewsakun, Dann Turner, Evelien Adriaenssens and Peter Simmonds.

Based on original GRAViTy framework https://github.com/PAiewsakun/GRAViTy

Please cite:

Mayne, R., Aiewsakun, P., Turner, D. et al. (2024) GRAViTy-V2: a grounded viral taxonomy application. Preprint. https://www.biorxiv.org/content/10.1101/2024.07.26.605250v1

Aiewsakun, P., Simmonds, P. The genomic underpinnings of eukaryotic virus taxonomy: creating a sequence-based framework for family-level virus classification. Microbiome 6, 38 (2018). https://doi.org/10.1186/s40168-018-0422-7

## Background
GRAViTy (genome relationships applied to virus taxonomy) is a framework for identifying and classifying virus genomes. GRAViTyV2 is a development on the original framework written in Python 3, which builds in manifold optimizations, wider support, machine learning component refinements, an interactive API, Docker container compatibility, automatic acquisition and filtering of virus metadata resource (VMR) files and a different workflow.

The software as provided may either be run from the shell/other command line, or through using a browser-based GUI.

In spite of best efforts to simplify the user experience, GRAViTy-V2 is a complex piece of software: please read the user guide and journal article before use, as this will likely help you to achieve optimal results and reduce software run time.

![GRAViTy-V2 Process Flow](docs/gravity_flow_v2.png "GRAViTy-V2 Process Flow")

## Installation instructions:
This guide is tested on Windows Subsystems Linux 22.04, Ubuntu 22.04 LTS and Ubuntu Server 22.04 LTS. We anticipate that variations in operating system and system architecture (in particular, servers) will necessitate amendments to the process outlined below, hence this guide is advisory only. Always check the terminal output. This guide is not a substite for computing expertise, consequently users who don't understand the installation instructions should defer to their system administrator's guidance.

1. Ensure you have Conda installed: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
1. Create a new Conda environment called `gravity` with python 3.10: ```conda create -n gravity python=3.10```
1. Activate your Conda environment: ```conda activate gravity```
1. Navigate to directory: ```cd {path-to-dir}```
1. Install command line tool requirements: ```bash install-reqs-local.sh```
1. Install Python libraries with pip: ```pip install -r requirements.txt```
1. (Bash only) ```$ echo "alias gravity='conda activate castanet && uvicorn app.api:app --reload --port 8000'" >> ~/.bashrc```, then refresh your terminal with ```$ source ~/.bashrc```.
1. Check that installation of libraries and third party tools was successful. Activate your GRAViTy-V2 Conda environment with ```$ conda activate gravity```, then call the dependency test with ```$ python3 -m app.test.dep_test```. Check the terminal output for success/error messages.
1. Trigger an example run using the dataset included in the repository (in ./data/eval/). ```$ python3 -m app.test.example_run```. Check the terminal output for success/error messages. Output will be saved to ./output/GRAViTyV2_example_run.

## Quick-start: GUI
This guide section covers use of the three GRAViTy-V2 premade pipelines, using the graphical user interface (GUI). These pipelines are adapted for common use cases --- datasets comprising similar viruses, divergent viruses and viruses with long, single ORFs --- and require a minimum of user input to run. This guide specifically examines using the `similar_viruses` pipeline, but the other two premade pipelines function in exactly the same manner.

1. Navigate to your GRAViTy-V2 directory with ```$ cd```
1. Start the GRAViTy-V2 server with ```$ gravity``` (or, if not using a Bash shell/not followed the final step of installation instructions, activate your Conda environment with ```$ conda activate gravity``` then start the server with ```$ uvicorn app.api:app --reload --port 8000```)
1. Visit the graphical user interface in your browser at ```http://127.0.0.1:8000/docs```
1. Trigger an example run using the dataset included in the repository. Expand the `similar_viruses` box, then hit the "try it now" button to make the dialogue box editable.
1. Input the following three fields, which are described below. N.b. you must preserve the formatting of the dialogue! This means that double quote marks, colons, commas, brackets etc. must not be deleted or moved. If you want to run this on your own data, edit the field in double quote marks AFTER the colon.
   * GenomeDescTableFile: This is your VMR-like document, which must be in CSV format and contain at least the columns: [Realm,Order,Family,Genus,Species,Exemplar or additional isolate,Virus name(s),Virus name abbreviation(s),Virus isolate designation,Virus GENBANK accession,Virus REFSEQ accession,Genome coverage,Genome composition,Host source,Baltimore Group,Genetic code table]. Refer to the ICTV VMR (https://ictv.global/vmr) and example file referenced below for format.
   * ExpDir: This is the directory to save your experiment output. If it doesn't already exist, GRAViTy-V2 will create it.
   * GenomeSeqFile: This is the location of your GenBank file, containing all of the sequence data for the viruses you are analysing. If it doesn't already exist, GRAViTy-V2 will attempt to download all of the sequence data according to accesion numbers in your VMR-like document (GenomeDescTableFile), then save it to this file name.
```
{
  "GenomeDescTableFile": "./data/eval/example_vmr.csv",
  "ExpDir": "./output/quick_start_example",
  "GenomeSeqFile": "./output/quick_start_example.gb"
}
```
1. Observe the output in the terminal and check for errors. When the analysis has finished, examine the output in your ./output/quick_start_example directory.

## Quick-start: CLI
CLI readme to follow shortly!

## GRAViTy-V2 Custom Pipelines
GRAViTy-V2 custom pipelines offer total control over a user's workflow. Entrypoints may be found in the GUI under the "Custom pipelines" header. Users can run the full pipeline (`pipeline_i_full`), or from individual stages in the pipeline.

Descriptions for each parameter may be found below. Argument typing (i.e. expected input formats) for each may be found at the bottom of the GUI page, in the `schemas` section.

* GenomeDescTableFile: Path to VMR-like document (CSV)
* ExpDir: Directory to save output data. If doesn't exist, GRAViTy-V2 will create.
* GenomeSeqFile: Path to GenBank file containing sequence data. If doesn't exist, GRAViTy-V2 will attempt to download from NCBI Entrez, using virus accession IDs from input VMR (GenomeDescTableFile)
* TaxoGrouping_Header: Column title in VMR to draw taxon boundaries around in output. Choose from "Family", "Genus" "Taxonomic grouping", where latter is a custom column that users can create in the input CSV document.
* genbank_email: Provide an email address for use in the NCBI Entrez sequence download system.
* ProteinLength_Cutoff: Miminum length (in AA) for a protein sequence to be used in GRAViTy-V2 analysis. Recommended range: 70-100.
* Mash_p_val_cutoff: p value above which Mash hit will be ignored (as part of initial ORF clustering, during PPHMMDb construction). Recommended to not change.
* Mash_sim_score_cutoff: similarity score below which Mash hit will be ignored (as part of initial ORF clustering, during PPHMMDb construction). Recommended range 0.7-0.95.
* ProtClustering_MCLInflation: mcl algorithm inflation score (as part of initial ORF clustering, during PPHMMDb construction). Recommended range 2.0-4.0.
* N_AlignmentMerging: Number of iterations for alignment merging rounds (as part of initial ORF clustering, during PPHMMDb construction). Warning: this can dramatically increase run times and has only a few, highly specific applications. Recommended to not change.
* PPHMMClustering_MCLInflation_ForAlnMerging: mcl algorithm inflation score (as part of alignment merging). Argument is ignored if N_AlignmentMerging == 0. Recommended to not change.
* HHsuite_evalue_Cutoff: [IF SORTING PPHMMS]
* HHsuite_pvalue_Cutoff:
* HHsuite_QueryCoverage_Cutoff:
* HHsuite_SubjectCoverage_Cutoff:
* HMMER_C_EValue_Cutoff: Hmmer E value above which to ignore hits (as part of PPHMM signature generator). Recommended to not change.
* HMMER_HitScore_Cutoff: Hmmer hitscore below which to ignore hits (as part of PPHMM signature generator). Recommended range 0-40 (use will be contextual; 25 is a good starting point, although 0 is default).
* p: Scale factor for similarity matrix. Recommended to not change.
* Dendrogram_LinkageMethod: Select from ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
* Bootstrap: If true, calculate bootstrap support for dendrogram.
* N_Bootstrap: If Bootstrap is true, number of iterations to calculate.
* Bootstrap_method: Method for calculating bootstrap support. Select from ["sumtrees", "booster"]
* VirusGrouping: Perform virus grouping (as part of virus classification). Recommended true.
* SimilarityMeasurementScheme: Scheme for constructing similarity matrix. Choose from: ["P", "G", "L", "PG", "PL", "R", "RG", "PR"], based on whether PPHMM signatures (P), GOM signatures (G), distance correlation (L), PPHMM ratios (R) or combinations of these are used. Recommended "PG".
* Heatmap_DendrogramSupport_Cutoff: Threshold for bootstrap support to be shown on GRAViTY-V2 heatmap dendrograms.
* PphmmNeighbourhoodWeight: Apply weighting to profile scores when all genomes being analysed are highly similar (as part of similarity matrix construction). Recommended to not change.
* PphmmSigScoreThreshold: Disregard protein profiles where signature scores are below (as part of similarity matrix construction). Recommended to not change.
* UseBlast: If true, use BLASTp instead of Mash (as part of initial ORF clustering, during PPHMMDb construction). Recommended when comparing datasets comprising mostly very similar viruses (e.g. all within a single family).
* NThreads: Number of threads for multiprocessing. Choose from either an integer or "auto". Recommended to not change, unless users experience out of memory errors or are working on HPCs.
* ClustAlnScheme: Determines which algorithm is used to do protein alignments (as part of initial ORF clustering, during PPHMMDb construction). Choose from ["local", "global", "auto"], to correspond to Mafft L-INS-i, G-INS-i and auto algorithms, respectively. Recommended to "auto" if unsure, as this is the fastest option.
* RemoveSingletonPPHMMs: If true, disregard scores from PPHMMs shared by only *n* genomes (see N_VirusesOfTheClassToIgnore). Recommended to not change.
* N_VirusesOfTheClassToIgnore: If true and RemoveSingletonPPHMMs true, number of genomes above which must possess a PPHMM for it to be included in the analysis. Recommended to not change.
* PPHMMSorting: If true, sort PPHMMs (as part of PPHMMDb construction). This can dramatically increase run times and will rarely increase performance. Recommended false.
* N_Sampling: Number of mutual information scores sample size. Recommended to not change.
* SamplingStrategy: If N_Sampling > 0, strategy for calculating mutual information.
* SampleSizePerGroup: If N_Sampling > 0, sample size for each group.


## GRAViTy-V2 Workflows
GRAViTyV2 is a resource-intensive application and it would be intractable to construct databases comprising every known viral genome. For this reason, we advise breaking down queries into several stages at differing levels of granularity: we call this the "pass" system.

Imagine you have an RNA virus sequence with unknown taxonomy and you have minimal prior knowledge of its provenance. A "first pass" might compare your unknown sequence against representative members of each family in *Riboviria*. Your results then indicate that your unknown sequence is most similar to the sequence from *Flaviviridae*. A sensible consequent "second pass" would therefore be to re-run GRAViTy-V2, comparing your unknown sequence to a wide variety of genomes *Flaviviridae*, which will result in a much more precise taxonomy.

## Troubleshooting
The following information is designed to help troubleshoot common issues.
1. *GRAViTyV2 can't find a VMR file*. Download an up-to-date VMR via the scrape_vmr endpoint, or construct your own and point to it with the API call parameter "GenomeDescTableFile".
1. *Can't access GRAViTyV2 API in browser*. Ensure program is running in shell/Docker container. Ensure address (localhost/IP of Docker container) and port are correct; if deployed on a server, network or cloud resource, contact your administrator as complex infrastructure such as network loops and proxies are likely to be in place.
1. *GRAViTyV2 is taking a very long time to compute and/or crashes with memory overflow errors*. The software is extremely resource-intensive and incorrect usage will lead to intractable compute jobs. Of particular note, ensure that database size << 1000 genomes and consider disabling optional steps which are extremely expensive (e.g. bootstrapping, mutual information calculation, use of unclassified virus PPHMMs). Ensure you have specified the maximum number of threads you have available for the ```NThreads``` parameter, which will ensure maximum efficiency on BLAST, HMMER, HHSuite and Bootstrap functions. The default values provided for BLASTp, HMMER, HHSuite and Bootstrap functions are for guidance only: thoughtful optimization of these parameters can reduce run time by several orders of magnitude.
1. *How do I find out what all of the parameters for the API functions mean?* Scroll down on the API front-end to the "Schemas" section. Each function description can be expanded to give details of each parameter's type and functionality. High level descriptions of the functions and the rationale for their use are included in the documentation (./docs).
1. *Labels in my output heatmap/dendrogram don't seem to correspond with my input VMR labels* Ensure that there are no commas in your free text fields in the input VMR: GRAViTyV2 works by reading in comma-delimited (CSV) files, which struggle when commas are present in strings. The simplest, non-programmatic way of doing this is to do find-replace all commas with semicolons, via Excel or similar spreadsheet software.

Can't find an answer to your issue? Try emailing us or opening an issue :)

## Guide for contributors
Although forking is encouraged, we will only consider pull requests which address bugs and performance issues. Contributors will please configure pre-commit hooks to match ours, as detailed in the .pre-commit-config.yaml file.

## Disclaimer
The material embodied in this software is provided to you "as-is", “with all faults”, and without warranty of any kind, express, implied or otherwise, including without limitation, any warranty of fitness for a particular purpose, warranty of non-infringement, or warranties of any kind concerning the safety, suitability, lack of viruses, inaccuracies, or other harmful components of this software. There are inherent dangers in the use of any software, and you are solely responsible for determining whether this software is compatible with your equipment and other software installed on your equipment. You are convert_fasta_to_genbankalso solely responsible for the protection of your equipment and backup of your data, and the developers/providers will not be liable for any damages you may suffer in connection with using, modifying, or distributing this software. Without limiting the foregoing, the developers/providers make no warranty that: the software will meet your requirements; the software will be uninterrupted, timely, secure, or error-free; the results that may be obtained from the use of the software will be effective, accurate, or reliable; the quality of the software will meet your expectations; any errors in the software will be identified or corrected.

Software and its documentation made available here could include technical or other mistakes, inaccuracies, or typographical errors. The developers/providers may make changes to the software or documentation made available here may be out of date, and the developers/providers make no commitment to update such materials.

The developers/providers assume no responsibility for errors or omissions in the software or documentation available from here.

In no event shall the developers/providers be liable to you or anyone else for any direct, special, incidental, indirect, or consequential damages of any kind, or any damages whatsoever, including without limitation, loss of data, loss of profit, loss of use, savings or revenue, or the claims of third parties, whether or not the developers/providers have been advised of the possibility of such damages and loss, however caused, and on any theory of liability, arising out of or in connection with the possession, use, or performance of this software.

The use of this software is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system, or networked devices, or loss of data that results from such activities. No advice or information, whether oral or written, obtained by you from the developers/providers shall create any warranty for the software.
