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
GRAViTy V2 adaptation by Mayne, R., Aiewsakun, P., Simmonds., P. et al. (2022)

Based on original GRAViTy software https://github.com/PAiewsakun/GRAViTy

Please cite:

Aiewsakun, P., Simmonds, P. The genomic underpinnings of eukaryotic virus taxonomy: creating a sequence-based framework for family-level virus classification. Microbiome 6, 38 (2018). https://doi.org/10.1186/s40168-018-0422-7

## Background
GRAViTy (genome relationships applied to virus taxonomy) is a framework for identifying and classifying virus genomes. GRAViTyV2 is a development on the original framework written in Python 3, which builds in manifold optimizations, wider support, machine learning component refinements, an interactive API, Docker container compatibility, automatic acquisition and filtering of virus metadata resource (VMR) files and a different workflow.

The software as provided may either be run from the shell/other command line, or via a Docker container. Both methods have their own benefits and detriments; the Docker version will automatically set up your environment and dependencies, both of which need to be set up manually if running from the shell, but the former will also require knowledge of containerisation to access the application programming interface (API), import/export data and make any changes to the code base. If unsure, consult your compute cluster's administrator as to which option is preferable.

We assume for the purposes of this readme that the user will interact with GRAViTyV2's API via its graphical user interface (GUI), either via a localhost or LAN/WAN connection. This being said, the user may also interact with the API via the command line (e.g. with curl) or a third-party tool, such as Postman.

GRAViTy comprises two distinct pipelines. Pipeline 1 (PL1) is for analysing and creating "databases" of information pertaining to these analyses, on "reference" viruses, i.e. genomes for which we have previously generated taxonomic data. Pipeline 2 (PL2) is for identifying and classifying "unclassified" viruses, by doing comparisons with the genomes in databases generated during PL1. Both pipelines generate additional statistics and visualisations, such as heatmaps. Users will please refer to our documentation (./docs) and publications for further details.

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
xxxxxxxxxxxxx

## GRAViTyV2 Workflows
GRAViTyV2 is a resource-intensive application and it would be intractable to construct databases comprising every known viral genome. For this reason, we advise breaking down queries into several stages which should begin with constructing several databases based on significantly restricted datasets using Pipeline I, then comparing your samples against these with Pipeline II to get an approximate indication of taxonomy (hereafter known as doing a FIRST PASS run). First pass databases may be constructed with the scrape_vmr and construct_first_pass_vmr endpoints; we may choose to create several first pass databases where a simple filtering function may be used to ensure that single representative members of each family are included, or otherwise single members of each genus or species are included if upper taxonomic ranks contain only a small number of examples. It is good practice to keep databases under 800 samples, as compute time scales non-linearly with sample size: we may choose therefore, to find a logical way to split first pass databases into multiple portions, such as into RNA, dsDNA and ssDNA. First pass databases take the longest time to compute, but need only be constructed each time a new VMR is issued.

Consequently, a far-more detailed database can be constructed by running Pipeline I again, using a small subset of viral genomes that are known to be similar to the families/genuses indicated by the results of the first pass search. E.g. if the first pass indicates that an unknown sample is likely to belong to *Caulimoviridae*, a SECOND PASS database may be constructed using Pipeline I, to include all examples from this Family. Then, Pipeline II may be run against this second pass database to get a precise taxonomic evaluation.

In summary, an efficient GRAViTy workflow could look like the following:
1. Scrape VMR, then run first pass filter function to reduce number of samples. It may be necessary to repeat this step several times if splitting your first pass databases into smaller segments, e.g. by nucleic acid type.
1. Run Pipeline I several times to construct first pass databases.
1. Run Pipeline II against your first pass databases to get an approximate estimation of taxonomy.
1. Run the second pass filter function to create a second pass database, using a subset of the VMR that contains samples that the first pass database experiments have indicated are genetically similar with your unknown samples.
1. Generate your second pass database by running Pipeline I.
1. Run Pipeline II using your second pass database, to get a precise taxonomic evaluation.

Example workflows are included in the documentation (./docs). N.b. no two workloads should look the same and many run parameters, including the number of passes (more than two may be required), so the user should ensure they understand how GRAViTy works before they begin designing workflows.

## Troubleshooting
The following information is designed to help troubleshoot common issues.
1. *GRAViTyV2 can't find a VMR file*. Download an up-to-date VMR via the scrape_vmr endpoint, or construct your own and point to it with the API call parameter "GenomeDescTableFile".
1. *GRAViTyV2 fails just after a subprocess.Popen() command, such as BLASTp*. Ensure that all command line tools (mcl, hmmer, ncbi-blast+, hhsuite, muscle, booster) are installed, can be called from the command line and are installed in the locations GRAViTyV2 expects (see paths in ./install-reqs.sh). Ensure versions installed are compatible with your operating system and architecture (in particular, your CPU instruction set if using a server).
1. *Can't access GRAViTyV2 API in browser*. Ensure program is running in shell/Docker container. Ensure address (localhost/IP of Docker container) and port are correct; if deployed on a server, network or cloud resource, contact your administrator as complex infrastructure such as network loops and proxies are likely to be in place.
1. *GRAViTyV2 is taking a very long time to compute and/or crashes with memory overflow errors*. The software is extremely resource-intensive and incorrect usage will lead to intractable compute jobs. Of particular note, ensure that database size << 1000 genomes and consider disabling optional steps which are extremely expensive (e.g. bootstrapping, mutual information calculation, use of unclassified virus PPHMMs). Ensure you have specified the maximum number of threads you have available for the ```N_CPUs``` parameter, which will ensure maximum efficiency on BLAST, HMMER, HHSuite and Bootstrap functions. The default values provided for BLASTp, HMMER, HHSuite and Bootstrap functions are for guidance only: thoughtful optimization of these parameters can reduce run time by several orders of magnitude.
1. *How do I specify "nothing" for a parameter in the API, e.g. I don't want to specify a ```Database``` parameter*. The API accepts JSON-like objects: to specify no input for these parameters, use ```null```.
1. *How do I find out what all of the parameters for the API functions mean?* Scroll down on the API front-end to the "Schemas" section. Each function description can be expanded to give details of each parameter's type and functionality. High level descriptions of the functions and the rationale for their use are included in the documentation (./docs).
1. *Labels in my output heatmap/dendrogram don't seem to correspond with my input VMR labels* Ensure that there are no commas in your free text fields in the input VMR: GRAViTyV2 works by reading in comma-delimited (CSV) files, which struggle when commas are present in strings. The simplest, non-programmatic way of doing this is to do find-replace all commas with semicolons, via Excel or similar spreadsheet software.
1. *Muscle is returning error "Segmentation fault"* Ensure muscle has installed properly by trying to run it through the CLI. If the program won't respond at all, check you have the correct version installed (i.e. operating system and bit type) and manually reinstall, see https://drive5.com/muscle/manual/. If you're having trouble, we highly recommend using Conda (https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to automatically install this package separately, with ```conda install -c bioconda muscle=3.8```. Conversely, if the program will initialise but a Segmentation Fault returns mid-job, this usually implies that your machine has run out of memory during calculations.

## Guide for contributors
Although forking is encouraged, we will only consider pull requests which address bugs and performance issues. Contributors will please configure pre-commit hooks to match ours, as detailed in the .pre-commit-config.yaml file.

## Change log
09.06.23
1. Muscle downgrade to better support work involving extremely long sequences
1. Fasta to Genbank conversion fixes
1. Input Accession IDs will now not be looked up on Entrez if the user specifies a pre-computed Genbank file, which reduces chances of accidental Accession ID matches with user-specified naming conventions
1. Enhanced support for Booster bootstrapping
1. Segment concatenation endpoint now arranges sequences by size
1. Figure size and component scaling is now partially automated, which enhances readability of GRAViTy heatmaps when large numbers of samples are used
1. General code tidying and quality of life changes

## Disclaimer
The material embodied in this software is provided to you "as-is", “with all faults”, and without warranty of any kind, express, implied or otherwise, including without limitation, any warranty of fitness for a particular purpose, warranty of non-infringement, or warranties of any kind concerning the safety, suitability, lack of viruses, inaccuracies, or other harmful components of this software. There are inherent dangers in the use of any software, and you are solely responsible for determining whether this software is compatible with your equipment and other software installed on your equipment. You are convert_fasta_to_genbankalso solely responsible for the protection of your equipment and backup of your data, and the developers/providers will not be liable for any damages you may suffer in connection with using, modifying, or distributing this software. Without limiting the foregoing, the developers/providers make no warranty that: the software will meet your requirements; the software will be uninterrupted, timely, secure, or error-free; the results that may be obtained from the use of the software will be effective, accurate, or reliable; the quality of the software will meet your expectations; any errors in the software will be identified or corrected.

Software and its documentation made available here could include technical or other mistakes, inaccuracies, or typographical errors. The developers/providers may make changes to the software or documentation made available here may be out of date, and the developers/providers make no commitment to update such materials.

The developers/providers assume no responsibility for errors or omissions in the software or documentation available from here.

In no event shall the developers/providers be liable to you or anyone else for any direct, special, incidental, indirect, or consequential damages of any kind, or any damages whatsoever, including without limitation, loss of data, loss of profit, loss of use, savings or revenue, or the claims of third parties, whether or not the developers/providers have been advised of the possibility of such damages and loss, however caused, and on any theory of liability, arising out of or in connection with the possession, use, or performance of this software.

The use of this software is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system, or networked devices, or loss of data that results from such activities. No advice or information, whether oral or written, obtained by you from the developers/providers shall create any warranty for the software.
