#!/bin/bash

# make folders
mkdir output

# mcl installation
apt-get install mcl -y

# hmmer installation
apt-get install hmmer -y

# blast installation. Assumes Linux x64 build!
apt install ncbi-blast+

# RM < TODO Install MAFFT

# RM < TODO install mash

# hhsuite installation # RM < TODO deprecate for conda install as some users have issues
mkdir -p ~/programs/hh-suite && cd ~/programs/hh-suite
git clone https://github.com/soedinglab/hh-suite.git .
mkdir build && cd build
cmake ..
make

# booster installation
wget https://github.com/evolbioinfo/booster/releases/download/v0.1.2/booster_linux64
chmod 777 booster_linux64
