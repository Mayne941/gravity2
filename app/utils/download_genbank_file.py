from Bio import Entrez
import os

from app.utils.error_handlers import raise_gravity_error

def DownloadGenBankFile(GenomeSeqFile, SeqIDLists, email):
    '''Hit GenBank to get data. Needs user to provide email authentication'''
    if not os.path.exists("/".join(GenomeSeqFile.split("/")[:-1])):
        os.makedirs("/".join(GenomeSeqFile.split("/")[:-1]))

    Entrez.email = email
    try:
        handle = Entrez.efetch(db="nucleotide", id=", ".join([", ".join(x) for x in SeqIDLists]), rettype="gb", retmode="text")
    except Exception as e:
        raise_gravity_error(f"Failed to pull Genbank Data from NCBI Entrez with exception: {e}"
                         f"This is usually a temporary problem due to NCBI server down time, try again in a few minutes!")
    print("Writing contents of genbank object to file")
    with open(GenomeSeqFile, "w") as GenomeSeqFile_handle:
        GenomeSeqFile_handle.write(handle.read())

    handle.close()
