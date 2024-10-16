from app.utils.console_messages import section_header
from app.utils.generate_fnames import generate_file_names
from app.utils.error_handlers import raise_gravity_error, raise_gravity_warning

import pandas as pd
import numpy as np
from collections import Counter
import os
import re
import pickle
import csv


class ReadGenomeDescTable:
    '''Parse, transform and load input VMR to GRAViTy compatible format, store as shelve'''
    def __init__(self,
                 payload,
                 GenomeDescTableFile,
                 GenomeSeqFile,
                 ExpDir,
                 RefreshGenbank=False,
                 is_secondpass=False
                 ) -> None:
        self.payload = payload
        self.GenomeDescTableFile = GenomeDescTableFile
        self.fnames = generate_file_names(payload,ExpDir)
        self.GenomeSeqFile = GenomeSeqFile
        if not os.path.exists(self.fnames["OutputDir"]):
            os.makedirs(self.fnames["OutputDir"])
        self.refresh_genbank = RefreshGenbank
        self.is_secondpass = is_secondpass
        self.no_acc_cnt = 1
        self.BaltimoreList, self.OrderList, self.FamilyList, \
            self.SubFamList, self.GenusList, self.VirusNameList, \
            self.SeqIDLists, self.SeqStatusList, self.TaxoGroupingList, \
            self.TranslTableList, self.DatabaseList, self.VirusIndexList,\
            self.transl_table_errors = [], [], [], [],\
                [], [], [], [], [], [], [], [], []

    def open_table(self) -> None:
        '''Open VMR file and read in relevant data to volatile'''
        print("- Read the GenomeDesc table")
        try:
            df = pd.read_csv(self.GenomeDescTableFile, index_col=0)
        except:
            raise_gravity_error(f"Failed to read your input VMR-like document (GenomeDescTableFile). Check it's a valid CSV.")
        df = df.fillna("")
        crap_pandas = [i for i in df.columns if "Unnamed" in i]
        for i in crap_pandas:
            df = df.drop(columns=i)

        try:
            df["Virus name(s)"] = df.apply(lambda x: self.get_names(x), axis=1)
            df["Genetic code table"] = df.apply(lambda x: self.transl_table_check(x), axis=1)
            df["Virus GENBANK accession"] = df.apply(lambda x: self.get_accession(x), axis=1)
        except KeyError as e:
            raise_gravity_error(f"Your input VMR-like document isn't structured appropriately, see documentation for instructions.\nException: {e}")

        self.VirusIndexList = [i for i in range(1, df.shape[0] + 1)]
        try: # TODO < TIDY THIS
            self.BaltimoreList = df["Baltimore Group"].tolist()
        except KeyError:
            self.BaltimoreList = ["" for i in df["Virus GENBANK accession"].tolist()]
        try:
            self.OrderList = df["Order"].tolist()
        except KeyError:
            self.OrderList = ["" for i in df["Virus GENBANK accession"].tolist()]
        self.FamilyList = df["Family"].tolist()
        self.SubFamList = df["Subfamily"].tolist()
        self.GenusList = df["Genus"].tolist()
        self.VirusNameList = df["Virus name(s)"].tolist()
        try:
            self.TaxoGroupingList = df[self.payload['TaxoGrouping_Header']].tolist()
        except KeyError:
            raise_gravity_error(f"You specified to group your input VMR by the column header: {self.payload['TaxoGrouping_Header']}"
                                f"but GRAViTy couldn't find this column in the file! Please check and amend your TaxoGrouping_Header parameter.")
        try:
            self.SeqStatusList = df["Genome coverage"].tolist()
        except KeyError:
            self.SeqStatusList = [1 for i in df["Virus GENBANK accession"].tolist()]
        self.SeqIDLists = [i.split(", ") for i in df["Virus GENBANK accession"].tolist()]
        self.TranslTableList = df["Genetic code table"].tolist()

    def get_names(self, row):
        # try: # RM < TODO Incorporate commented code for curtailing long names?
        #     if not self.is_secondpass:
        #         if not row["Virus name(s)"]:
        #             row["Virus name(s)"] == ""
        #         return re.sub(r"^\/|\/$", "", re.sub(r"[\/ ]{2,}", "/", re.sub(
        #             r"[^\w^ ^\.^\-]+", "/", re.sub(r"[ ]{2,}", " ", row["Virus name(s)"]))))
        #     else:
        #         '''If PL2, append virus description to accession ID'''
        #         part1 = re.sub(r"^\/|\/$", "", re.sub(r"[\/ ]{2,}", "/", re.sub(r"[^\w^ ^\.^\-]+", "/", re.sub(r"[ ]{2,}", " ", row["Virus name(s)"]))))
        #         if len(row["Virus isolate designation"]) == 0:
        #             part2 = ""
        #         part2 = re.sub(r"[^\w\s]", "", row["Virus isolate designation"][0:39]).replace(" ", "_")
        #         if len(row["Virus isolate designation"]) > 40:
        #             part2 += "..."
        #         return f'{part1}_{part2}'
        try:
            if self.is_secondpass:
                if not row["Virus name(s)"]:
                    row["Virus name(s)"] == ""
                virus_name = re.sub(r"^\/|\/$", "", re.sub(r"[\/ ]{2,}", "/", re.sub(r"[^\w^ ^\.^\-]+", "/", re.sub(r"[ ]{2,}", " ", row["Virus name(s)"])))).replace(" ","_")
                accession = ", ".join(re.findall(r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{4}[0-9]{6,8}|[A-Z]{2}_[0-9]{6}|SRR[0-9]{7,8}", row["Virus GENBANK accession"]))
                if len(accession.split(",")) > 1:
                    accession = accession.split(",")[0]
                return f'{accession}_{row["Family"]}_{row["Genus"]}_{virus_name}'
            else:
                if not row["Virus name(s)"]:
                    row["Virus name(s)"] == ""
                return re.sub(r"^\/|\/$", "", re.sub(r"[\/ ]{2,}", "/", re.sub(
                    r"[^\w^ ^\.^\-]+", "/", re.sub(r"[ ]{2,}", " ", row["Virus name(s)"]))))

        except Exception as ex:

            raise_gravity_error(f"At least one line in your input CSV (genome desc table) has empty fields (or fields causing another error): {ex}")

    def transl_table_check(self, row):
        try:
            if row["Genetic code table"] != 1:
                if not "complete" in str(row["Genome coverage"]).lower():
                    self.transl_table_errors.append(row["Virus name(s)"])
                return int(1)
            else:
                return int(row["Genetic code table"])
        except:
            raise_gravity_warning(f"No 'Genetic code table' column found or invalid value therein for {row['Virus name(s)']}, so setting as default translation table")
            return 1

    def get_accession(self, row):
        try:
            # if "SRR" in row["Virus GENBANK accession"]:
            #     th = ", ".join()
            th = ", ".join(re.findall(r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{4}[0-9]{6,8}|[A-Z]{2}_[0-9]{6}|SRR[0-9]{7,8}", row["Virus GENBANK accession"]))
            if th == "":
                raise
            return th
        except:
            self.no_acc_cnt += 1
            raise_gravity_warning(f"No accession number for {row['Virus name(s)']}: this might cause GRAViTy to error later on if you're not providing your own GenBank file!")
            return f'{row["Virus GENBANK accession"]}' if not row["Virus GENBANK accession"] == "" else f'No_data_{self.no_acc_cnt - 1}'

    def update_desc_table(self) -> dict:
        '''Create dictionary in GRAViTy structure for saving to persistent storage'''
        master_data = {}
        if self.payload['AnnotateIncompleteGenomes']:
            '''If including incomplete, take all seqs'''
            IncludedGenomes_IndexList = [SeqStatus_i for SeqStatus_i, _ in enumerate(
                self.SeqStatusList)]
        else:
            '''Else only take coding complete seqs'''
            IncludedGenomes_IndexList = [SeqStatus_i for SeqStatus_i, _ in enumerate(
                self.SeqStatusList)]# if "Complete" in SeqStatus] # RM < TODO Removed as this can break a user workflow - guidance on using CCS is probably preferable

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
        duplicates = [SeqID for SeqID, count in Counter(
            SeqIDFlatList).items() if count > 1 and not SeqID == ""]
        if len(duplicates) > 0:
            '''Check if there exist duplicated accession numbers, fail IF TRUE'''
            print("The following accession numbers appear more than once: ")
            raise_gravity_error(
                f"Error parsing VMR."
                f"The following accession numbers appear more than once and GRAViTy only supports unique entries:\n{duplicates}"
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

        if self.no_acc_cnt > 1:
            raise_gravity_warning(f"{self.no_acc_cnt} genomes specified in your desc file had no accession IDs.")

        '''Create ReadGenomeDescTable "all genomes' db'''
        print("- Save variables to ReadGenomeDescTable pickle")
        all_desc_table = self.update_desc_table()
        self.save_desc_table(
            all_desc_table)
