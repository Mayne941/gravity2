def TaxoLabel_Constructor (SeqIDLists, FamilyList, GenusList, VirusNameList):
	return list(map("_".join,list(zip(["/".join(SeqIDList) if len(SeqIDList)<=3 else "/".join(SeqIDList[0:3])+"/..." for SeqIDList in SeqIDLists],
				[Family.replace(" ", "-") for Family in FamilyList],
				[Genus.replace(" ", "-") for Genus in GenusList],
				[VirusName.replace(" ", "-") for VirusName in VirusNameList],
				))
			))