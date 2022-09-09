#!/bin/bash

# mcl installation
apt-get install mcl -y

# hmmer installation
apt-get install hmmer -y

# blast installation. Assumes Linux x64 build!
# apt-get install ncftp -y
# ncftpget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.13.0+-x64-arm-linux.tar.gz
# tar zxvpf ncbi-blast-2.13.0+-x64-arm-linux.tar.gz
# export PATH=$PATH:$HOME/repo/ncbi-blast-2.13.0+/bin
apt install ncbi-blast+

# hhsuite installation
mkdir -p ~/programs/hh-suite && cd ~/programs/hh-suite
git clone https://github.com/soedinglab/hh-suite.git .
mkdir build && cd build
cmake ..
make
export HHLIB="~/programs/hh-suite"
export PATH="$PATH:$HHLIB/bin:$HHLIB/scripts"

# muscle installation. Assumes Linux x64 build!
mkdir -p ~/programs/muscle && cd ~/programs/muscle
wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64
cd src/
make

# booster installation
cd ~/workspace
wget https://github.com/evolbioinfo/booster/releases/download/v0.1.2/booster_linux64
