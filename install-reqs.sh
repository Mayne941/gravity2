#!/bin/bash

apt-get install mcl -y
apt-get install hmmer -y

#blast installation
apt-get install ncftp -y
ncftpget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.13.0+-x64-arm-linux.tar.gz
tar zxvpf ncbi-blast-2.13.0+-x64-arm-linux.tar.gz
export PATH=$PATH:$HOME/repo/ncbi-blast-2.13.0+/bin

#hhsuite installation
mkdir -p ~/programs/hh-suite && cd ~/programs/hh-suite
git clone https://github.com/soedinglab/hh-suite.git .
mkdir build && cd build
cmake ..
make


#muscle installation
wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64

