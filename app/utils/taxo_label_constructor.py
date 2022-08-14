def TaxoLabel_Constructor(SeqIDLists, FamilyList, GenusList, VirusNameList):
	'''Clean & join taxonomy labels in prep for graphing functions'''
	return list(map("_".join,list(zip(["/".join(SeqIDList) if len(SeqIDList)<=3 else "/".join(SeqIDList[0:3])+"/..." for SeqIDList in SeqIDLists],
				[Family.replace(" ", "-") for Family in FamilyList.astype("str")],
				[Genus.replace(" ", "-") for Genus in GenusList.astype("str")],
				[VirusName.replace(" ", "-") for VirusName in VirusNameList.astype("str")],
				))
			))