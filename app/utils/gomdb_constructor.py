
def GOMDB_Constructor (TaxoGroupingList, PPHMMLocationTable, GOMIDList):
	'''Generate genomic organisation model (GOM) database'''
	GOMDb = {}
	for id in GOMIDList:
		GOMDb[id] = PPHMMLocationTable[TaxoGroupingList == id,:]

	return GOMDb