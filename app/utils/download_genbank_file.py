from Bio import Entrez
import os


def DownloadGenBankFile(GenomeSeqFile, SeqIDLists, email):
    '''Hit GenBank to get data. Needs user to provide email authentication'''
    if not os.path.exists("/".join(GenomeSeqFile.split("/")[:-1])):
        os.makedirs("/".join(GenomeSeqFile.split("/")[:-1]))

    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=", ".join(
        [", ".join(x) for x in SeqIDLists]), rettype="gb", retmode="text")
    with open(GenomeSeqFile, "w") as GenomeSeqFile_handle:
        GenomeSeqFile_handle.write(handle.read())

    handle.close()
