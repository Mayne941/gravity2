#!/bin/bash
conda config --add channels conda-forge
conda config --add channels bioconda

# make folders
mkdir output

# mcl installation
conda install -y -c bioconda mcl

# hmmer installation
conda install -y -c bioconda hmmer

# blast installation
conda install -y -c bioconda blast

# samtools installation
conda install -y -c bioconda samtools

# mash installation
conda install -y -c bioconda mash

# hhsuite installation
conda install -y -c bioconda hhsuite

# booster installation
wget https://github.com/evolbioinfo/booster/releases/download/v0.1.2/booster_linux64
chmod 777 booster_linux64
