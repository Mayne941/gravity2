import numpy as np
import sys
from .dcor import dcor
from app.utils.stdout_utils import progress_bar, clean_stdout

def GOMSignatureTable_Constructor (PPHMMLocationTable, GOMDB, GOMIDList):
	N_Viruses		= len(PPHMMLocationTable)
	N_GOMs			= len(GOMIDList)
	GOMSignatureTable	= np.empty((N_Viruses,0))
	GOM_i			= 1
	for GOM in GOMIDList:
		GOMSignatureList, Virus_i = [], 1
		for PPHMMLocation in PPHMMLocationTable:
			RelevantPPHMMIndices = np.where(list(map(any, list(zip(list(map(any, GOMDB[GOM].transpose() != 0)), PPHMMLocation != 0)))))[0]
			GOMSignatureList.append(dcor(GOMDB[GOM][:, RelevantPPHMMIndices].T, PPHMMLocation[RelevantPPHMMIndices].reshape(-1,1)))

			progress_bar(f"\033[K GOM construction {GOM} ({GOM_i}/{N_GOMs}): [{'='*int(float(Virus_i)/N_Viruses*20)}] {Virus_i}/{N_Viruses} GOMs \r")
			Virus_i += 1
		
		clean_stdout()
		GOMSignatureTable = np.column_stack((GOMSignatureTable, GOMSignatureList))
		GOM_i = GOM_i+1
	
	return GOMSignatureTable

