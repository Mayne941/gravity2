from app.utils.console_messages import section_header
from app.utils.generate_fnames import generate_file_names

import numpy as np
from collections import Counter
import os
import re
import pickle


class ReadGenomeDescTable:
    '''Parse, transform and load input VMR to GRAViTy compatible format, store as shelve'''

    def __init__(self,
                 payload,
                 GenomeDescTableFile,
                 GenomeSeqFile,
                 ExpDir,
                 RefreshGenbank=False
                 ) -> None:
        self.payload = payload
        self.GenomeDescTableFile = GenomeDescTableFile
        self.fnames = generate_file_names(payload,ExpDir)
        self.GenomeSeqFile = GenomeSeqFile
        if not os.path.exists(self.fnames["OutputDir"]):
            os.makedirs(self.fnames["OutputDir"])
        self.refresh_genbank = RefreshGenbank
        self.BaltimoreList, self.OrderList, self.FamilyList, \
            self.SubFamList, self.GenusList, self.VirusNameList, \
            self.SeqIDLists, self.SeqStatusList, self.TaxoGroupingList, \
            self.TranslTableList, self.DatabaseList, self.VirusIndexList = [], [], [], [],\
            [], [], [], [],\
            [], [], [], []

    def open_table(self) -> None:
        '''Open VMR file and read in relevant data to volatile'''
        print("- Read the GenomeDesc table")
        with open(self.GenomeDescTableFile, "r") as GenomeDescTable_txt:
            header_raw = next(GenomeDescTable_txt)
            header = header_raw.replace("\n", "").replace(
                "\r", "").replace("\t", "").split(",")
            Baltimore_i = header.index("Baltimore Group")
            Order_i = header.index("Order")
            Family_i = header.index("Family")
            Subfamily_i = header.index("Subfamily")
            Genus_i = header.index("Genus")
            VirusName_i = header.index("Virus name(s)")
            SeqID_i = header.index("Virus GENBANK accession")
            SeqStatus_i = header.index("Genome coverage")
            TranslTable_i = header.index("Genetic code table")

            if self.payload['Database_Header'] != None:
                Database_i = header.index(self.payload['Database_Header'])

            if self.payload['TaxoGrouping_Header'] != None:
                TaxoGrouping_i = header.index(self.payload['TaxoGrouping_Header'])
            else:
                TaxoGrouping_i = header.index("Family")

            for Virus_i, Line in enumerate(GenomeDescTable_txt):
                Line = Line.replace("\n", "").replace(
                    "\r", "").replace("\t", "").split(",")
                if not os.path.isfile(self.GenomeSeqFile) or self.refresh_genbank:
                    try:
                        '''Regular exp to extract accession numbers'''
                        # print(f"*** You DIDN'T specify an input Genbank file, OR we're doing an end-to-end test, so I'm going to try to extract sequence IDs from the VMR")
                        SeqIDList = re.findall(
                            r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{4}[0-9]{6,8}|[A-Z]{2}_[0-9]{6}", Line[SeqID_i])
                    except IndexError as ex:
                        raise SystemExit(
                            f"Error parsing VMR: exception {ex}"
                            f"On line {Line}"
                            f"Accession ID could not be extracted."
                            f"GRAViTy terminated."
                        )
                else:
                    '''If we specify the genbank file, don't extract accession IDS and consequently download'''
                    # print(f"*** You DID specify an input Genbank file, so I'm NOT going to try to extract sequence IDs from the VMR")
                    SeqIDList = []

                if SeqIDList != [] or Line[SeqID_i] != "":
                    '''Ignore record without sequences'''
                    self.VirusIndexList.append(Virus_i)
                    self.BaltimoreList.append(Line[Baltimore_i])
                    self.OrderList.append(Line[Order_i])
                    self.FamilyList.append(Line[Family_i])
                    self.SubFamList.append(Line[Subfamily_i])
                    self.GenusList.append(Line[Genus_i])
                    self.VirusNameList.append(re.sub(r"^\/|\/$", "", re.sub(r"[\/ ]{2,}", "/", re.sub(
                        r"[^\w^ ^\.^\-]+", "/", re.sub(r"[ ]{2,}", " ", Line[VirusName_i])))))
                    self.TaxoGroupingList.append(Line[TaxoGrouping_i])
                    self.SeqStatusList.append(Line[SeqStatus_i])

                    if SeqIDList != []:
                        self.SeqIDLists.append(SeqIDList)
                    else:
                        self.SeqIDLists.append([Line[SeqID_i]])

                    try:
                        self.TranslTableList.append(int(Line[TranslTable_i]))
                    except ValueError:
                        print(
                            f"Genetic code is not specified. GRAViTy will use the standard code for {self.VirusNameList[-1]}")
                        self.TranslTableList.append(1)
                    if self.payload['Database_Header'] != None:
                        self.DatabaseList.append(Line[Database_i])


    def update_desc_table(self) -> dict:
        '''Create dictionary in GRAViTy structure for saving to persistent storage'''
        master_data = {}
        if self.payload['AnnotateIncompleteGenomes']:
            '''If including incomplete, take all seqs'''
            IncludedGenomes_IndexList = [SeqStatus_i for SeqStatus_i, SeqStatus in enumerate(
                self.SeqStatusList)]
        else:
            '''Else only take coding complete seqs'''
            IncludedGenomes_IndexList = [SeqStatus_i for SeqStatus_i, SeqStatus in enumerate(
                self.SeqStatusList) if "Complete" in SeqStatus]

        if self.payload['Database'] != None:
            '''Filter to user-provided database if present'''
            master_data["DatabaseList"] = self.DatabaseList[IncludedGenomes_IndexList]

        '''Create master data. Arrays are LINKED.'''
        self.BaltimoreList = master_data["BaltimoreList"] = self.BaltimoreList[IncludedGenomes_IndexList]
        self.OrderList = master_data["OrderList"] = self.OrderList[IncludedGenomes_IndexList]
        self.FamilyList = master_data["FamilyList"] = self.FamilyList[IncludedGenomes_IndexList]
        self.SubFamList = master_data["SubFamList"] = self.SubFamList[IncludedGenomes_IndexList]
        self.GenusList = master_data["GenusList"] = self.GenusList[IncludedGenomes_IndexList]
        self.VirusNameList = master_data["VirusNameList"] = self.VirusNameList[IncludedGenomes_IndexList]
        self.SeqIDLists = master_data["SeqIDLists"] = self.SeqIDLists[IncludedGenomes_IndexList]
        self.SeqStatusList = master_data["SeqStatusList"] = self.SeqStatusList[IncludedGenomes_IndexList]
        self.TranslTableList = master_data["TranslTableList"] = self.TranslTableList[IncludedGenomes_IndexList]
        self.TaxoGroupingList = master_data["TaxoGroupingList"] = self.TaxoGroupingList[IncludedGenomes_IndexList]
        self.DatabaseList = master_data["DatabaseList"] = self.DatabaseList
        return master_data

    def save_desc_table(self, table):
        '''Save dictionary in GRAViTy structure to persistent storage'''
        pickle.dump(table, open(self.fnames["ReadGenomeDescTablePickle"], "wb"))

    def entrypoint(self) -> None:
        '''RM < TODO DOCSTRING'''
        section_header("Read the GenomeDesc table")

        '''Open & parse input VMR (txt) file'''
        self.open_table()

        if self.payload['TaxoGroupingFile'] != None:
            '''PL1 has option to specify taxonomic grouping level.'''
            self.TaxoGroupingList = []
            with open(self.payload['TaxoGroupingFile'], "r") as TaxoGrouping_txt:
                for TaxoGrouping in TaxoGrouping_txt:
                    TaxoGrouping = TaxoGrouping.split("\r\n")[0].split("\n")[0]
                    self.TaxoGroupingList.append(TaxoGrouping)
            self.TaxoGroupingList = np.array(self.TaxoGroupingList)
            self.TaxoGroupingList = self.TaxoGroupingList[self.VirusIndexList]

        '''Generate rows for output'''
        all_desc_table = {}
        self.BaltimoreList = all_desc_table["BaltimoreList"] = np.array(
            self.BaltimoreList)
        self.OrderList = all_desc_table["OrderList"] = np.array(self.OrderList)
        self.FamilyList = all_desc_table["FamilyList"] = np.array(
            self.FamilyList)
        self.SubFamList = all_desc_table["SubFamList"] = np.array(
            self.SubFamList)
        self.GenusList = all_desc_table["GenusList"] = np.array(self.GenusList)
        self.VirusNameList = all_desc_table["VirusNameList"] = np.array(
            self.VirusNameList)

        '''Ensure SeqIDLists will be a h list (i.e. not a v list)'''
        self.SeqIDLists = np.array(self.SeqIDLists, dtype="object")

        SeqIDFlatList = [
            SeqID for SeqIDList in self.SeqIDLists for SeqID in SeqIDList]
        if len(SeqIDFlatList) != len(set(SeqIDFlatList)):
            '''Check if there exist duplicated accession numbers, fail IF TRUE'''
            print("The following accession numbers appear more than once: ")
            duplicates = [SeqID for SeqID, count in Counter(
                SeqIDFlatList).items() if count > 1]
            raise SystemExit(
                f"Error parsing VMR."
                f"The following accession numbers appear more than once: {duplicates}"
                f"GRAViTy terminated."
            )

        '''Update master data. Arrays are LINKED.'''
        all_desc_table["SeqIDLists"] = self.SeqIDLists
        self.SeqStatusList = all_desc_table["SeqStatusList"] = np.array(
            self.SeqStatusList)
        self.TaxoGroupingList = all_desc_table["TaxoGroupingList"] = np.array(
            self.TaxoGroupingList)
        self.TranslTableList = all_desc_table["TranslTableList"] = np.array(
            self.TranslTableList)
        self.DatabaseList = all_desc_table["DatabaseList"] = np.array(
            self.DatabaseList)

        '''Create ReadGenomeDescTable "all genomes' db'''
        print("- Save variables to ReadGenomeDescTable pickle")
        all_desc_table = self.update_desc_table()
        self.save_desc_table(
            all_desc_table)
