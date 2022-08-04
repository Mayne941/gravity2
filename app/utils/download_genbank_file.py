from Bio import Entrez
import os

def DownloadGenBankFile (GenomeSeqFile, SeqIDLists):
	'''Hit GenBank to get data. Needs user to provide email authentication'''
	if not os.path.exists("/".join(GenomeSeqFile.split("/")[:-1])):
		os.makedirs("/".join(GenomeSeqFile.split("/")[:-1]))
	
	Entrez.email = "director@mayneba.com" # Email addr is hard-coded - suggest specific addr used for each compute node
	handle = Entrez.efetch(db = "nucleotide", id=", ".join([", ".join(x) for x in SeqIDLists]), rettype="gb", retmode="text")
	with open(GenomeSeqFile, "w") as GenomeSeqFile_handle:
		GenomeSeqFile_handle.write(handle.read())
		
	handle.close()
