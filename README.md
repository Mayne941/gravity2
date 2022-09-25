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

## Installation instructions:
This guide is tested on Windows Subsystems Linux 20.04, Ubuntu 20.04 LTS and Ubuntu Server 20.04 LTS. We anticipate that variations in operating system and system architecture (in particular, servers) will necessitate amendments to the process outlined below, hence this guide is advisory only. Please see the Troubleshooting section for further details.

Where parameters appear in curly brackets, replace these with values that correspond to your environment.

### Shell (developed for Ubuntu Linux 20.04 LTS)
1. Navigate to directory: ```cd {path-to-dir}```
1. Install command line tool requirements: ```bash install-reqs.sh```
1. Install Python libraries with pip: ```pip install -r requirements.txt```
1. Source vars and start API: ```source env_vars.sh && python3 -m uvicorn app.api:app --reload```
1. Navigate to http://localhost:8000/docs in browser
1. Follow instructions on Swagger UI

### Docker
1. Download and install Docker, see https://docs.docker.com/get-docker/
1. Build and start Docker container: ```bash build.sh && bash run.sh```
1. Check Docker container has built successfully and is running: ```docker ps```
1. Navigate to API: http://{docker-container-ip}:8000/docs  (container IP may be found by in the Docker Desktop app, or through ```docker inspect {container-id}```, using container ID from previous step)

## Run instructions
xxxxxxxxXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxXXXXXXXXXx
. Instructions for navigating API.

## Troubleshooting
The following information is designed to help troubleshoot common issues.
1. *GRAViTyV2 can't find a VMR file*. Download an up-to-date VMR via the scrape_vmr endpoint, or construct your own and point to it with the API call parameter "GenomeDescTableFile".
1. *GRAViTyV2 fails just after a subprocess.Popen() command, such as BLASTp*. Ensure that all command line tools (mcl, hmmer, ncbi-blast+, hhsuite, muscle, booster) are installed, can be called from the command line and are installed in the locations GRAViTyV2 expects (see paths in ./install-reqs.sh). Ensure versions installed are compatible with your operating system and architecture (in particular, your CPU instruction set if using a server).
1. *Can't access GRAViTyV2 API in browser*. Ensure program is running in shell/Docker container. Ensure address (localhost/IP of Docker container) and port are correct; if deployed on a server, network or cloud resource, contact your administrator as complex infrastructure such as network loops and proxies are likely to be in place.
1. *GRAViTyV2 is taking a very long time to compute and/or crashes with memory overflow errors*. The software is extremely resource-intensive and incorrect usage will lead to intractable compute jobs. Of particular note, ensure that database size << 1000 genomes and consider disabling optional steps which are extremely expensive (e.g. bootstrapping, mutual information calculation, use of unclassified virus PPHMMs). Ensure you have specified the maximum number of threads you have available for the ```N_CPUs``` parameter, which will ensure maximum efficiency on BLAST, HMMER, HHSuite and Bootstrap functions.
1. *How do I specify "nothing" for a parameter in the API, e.g. I don't want to specify a ```Database``` parameter*. The API accepts JSON-like objects: to specify no input for these parameters, use ```null```.

## Guide for contributors
Although forking is encouraged, we will only consider pull requests which address bugs and performance issues. Contributors will please configure pre-commit hooks to match ours, as detailed in the .pre-commit-config.yaml file.

## Disclaimer
The material embodied in this software is provided to you "as-is", “with all faults”, and without warranty of any kind, express, implied or otherwise, including without limitation, any warranty of fitness for a particular purpose, warranty of non-infringement, or warranties of any kind concerning the safety, suitability, lack of viruses, inaccuracies, or other harmful components of this software. There are inherent dangers in the use of any software, and you are solely responsible for determining whether this software is compatible with your equipment and other software installed on your equipment. You are also solely responsible for the protection of your equipment and backup of your data, and the developers/providers will not be liable for any damages you may suffer in connection with using, modifying, or distributing this software. Without limiting the foregoing, the developers/providers make no warranty that: the software will meet your requirements; the software will be uninterrupted, timely, secure, or error-free; the results that may be obtained from the use of the software will be effective, accurate, or reliable; the quality of the software will meet your expectations; any errors in the software will be identified or corrected.

Software and its documentation made available here could include technical or other mistakes, inaccuracies, or typographical errors. The developers/providers may make changes to the software or documentation made available here may be out of date, and the developers/providers make no commitment to update such materials.

The developers/providers assume no responsibility for errors or omissions in the software or documentation available from here.

In no event shall the developers/providers be liable to you or anyone else for any direct, special, incidental, indirect, or consequential damages of any kind, or any damages whatsoever, including without limitation, loss of data, loss of profit, loss of use, savings or revenue, or the claims of third parties, whether or not the developers/providers have been advised of the possibility of such damages and loss, however caused, and on any theory of liability, arising out of or in connection with the possession, use, or performance of this software.

The use of this software is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system, or networked devices, or loss of data that results from such activities. No advice or information, whether oral or written, obtained by you from the developers/providers shall create any warranty for the software.
