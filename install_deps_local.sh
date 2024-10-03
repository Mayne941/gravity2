#!/bin/bash
conda config --add channels conda-forge
conda config --add channels bioconda

# make folders
mkdir output

# mcl installation
conda install -y -c bioconda mcl=22.282

# hmmer installation
conda install -y -c bioconda hmmer=3.4

# blast installation
conda install -y -c bioconda blast=2.16.0

# mash installation
conda install -y -c bioconda mash=2.3

# hhsuite installation
conda install -y -c bioconda hhsuite=3.3.0

# booster installation
wget https://github.com/evolbioinfo/booster/releases/download/v0.1.2/booster_linux64
chmod 777 booster_linux64
