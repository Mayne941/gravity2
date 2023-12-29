from sklearn.feature_selection import mutual_info_classif
import numpy as np
import os
import pickle

from app.utils.ordered_set import OrderedSet
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_pickle
from app.utils.stdout_utils import progress_bar, clean_stdout


class MutualInformationCalculator:
    def __init__(self,
                 ShelveDir,
                 IncludeIncompleteGenomes=False,
                 VirusGroupingFile=None,
                 N_Sampling=100,
                 SamplingStrategy="balance_with_repeat",
                 SampleSizePerGroup=10,
                 ):
        self.ShelveDir = ShelveDir
        self.VariableShelveDir = ShelveDir+"/output"
        self.IncludeIncompleteGenomes = IncludeIncompleteGenomes
        self.VirusGroupingFile = VirusGroupingFile
        self.N_Sampling = N_Sampling
        self.SamplingStrategy = SamplingStrategy
        self.SampleSizePerGroup = SampleSizePerGroup
        self.genomes, self.ref_annotations, self.pphmmdb_construction = {}, {}, {}

    def mkdirs(self) -> str:
        '''1/3: Return all dirs for db storage and retrieval'''
        MutualInformationScoreDir = self.VariableShelveDir+"/MutualInformationScore"
        if not os.path.exists(MutualInformationScoreDir):
            os.makedirs(MutualInformationScoreDir)

        return MutualInformationScoreDir

    def read_virus_grouping_file(self) -> dict:
        '''2/3: Read virus grouping file to dict'''
        VirusGroupingTable = []
        with open(self.VirusGroupingFile, "r") as VirusGrouping_txt:
            VirusGroupingSchemeList = next(VirusGrouping_txt)
            VirusGroupingSchemeList = VirusGroupingSchemeList.split("\r\n")[0].split("\n")[
                0].split("\t")
            for Line in VirusGrouping_txt:
                Line = Line.split("\r\n")[0].split("\n")[0].split("\t")
                VirusGroupingTable.append(Line)

        VirusGroupingTable = np.array(VirusGroupingTable).T
        VirusGroupingDict = {VirusGroupingScheme: np.array(
            VirusGroupingList) for VirusGroupingScheme, VirusGroupingList in zip(VirusGroupingSchemeList, VirusGroupingTable)}
        VirusGroupingDict = {"Overvall": self.genomes["TaxoGroupingList"]}
        return VirusGroupingDict

    def calculate_mis(self, PPHMMDesc, VirusGroupingDict, MutualInformationScoreDir) -> None:
        '''3/3: Calculate mutual info scores (sklearn), output to np to .txt, then save pickle'''
        ResultDict, N_VirusGroupingSchemes, VirusGroupingScheme_i = {}, len(
            VirusGroupingDict), 1
        for VirusGroupingScheme, VirusGroupingList in VirusGroupingDict.items():
            '''Compute mutual information between PPHMM scores and virus classification'''
            ResultDict[VirusGroupingScheme] = {}
            VirusGroupingIDList = [VirusGroupingID for VirusGroupingID in OrderedSet(
                VirusGroupingList) if VirusGroupingID != "-"]
            IngroupVirus_IndexList = np.where([(VirusGroupingID != "-" and not VirusGroupingID.startswith(
                "O_")) for VirusGroupingID in VirusGroupingList.astype("str")])[0]
            IngroupVirus_PPHMM_IndexList = np.where(np.sum(
                self.ref_annotations["PPHMMSignatureTable"][IngroupVirus_IndexList] != 0, axis=0) != 0)[0]
            IngroupVirus_PPHMMDesc = PPHMMDesc[IngroupVirus_PPHMM_IndexList]
            PPHMMSignatureTable_Subset = self.ref_annotations["PPHMMSignatureTable"][:,
                                                                                     IngroupVirus_PPHMM_IndexList]

            SampledMutualInformationTable = []
            for Sampling_i in range(self.N_Sampling):
                if self.SamplingStrategy == None:
                    SampledVirus_IndexList = np.where(
                        [VirusGroupingID != "-" for VirusGroupingID in VirusGroupingList])[0]

                elif self.SamplingStrategy == "balance_without_repeat":
                    SampledVirus_IndexList = sum([list(np.random.permutation(np.where(VirusGroupingList == VirusGroupingID)[0])[
                                                 :self.SampleSizePerGroup]) for VirusGroupingID in VirusGroupingIDList], [])

                elif self.SamplingStrategy == "balance_with_repeat":
                    SampledVirus_IndexList = sum([list(np.random.choice(a=np.where(VirusGroupingList == VirusGroupingID)[
                                                 0], size=self.SampleSizePerGroup, replace=True)) for VirusGroupingID in VirusGroupingIDList], [])

                '''Call SK-Learn mutual info classifier. This component scales badly but no replacement currently exists.'''
                SampledMutualInformationTable.append(mutual_info_classif(
                    PPHMMSignatureTable_Subset[SampledVirus_IndexList], VirusGroupingList[SampledVirus_IndexList]))

                progress_bar(
                    f"\033[K Virus grouping scheme = %{VirusGroupingScheme} ({VirusGroupingScheme_i}/{N_VirusGroupingSchemes}): [{'='*int(float(Sampling_i+1)/self.N_Sampling*20)}] {Sampling_i+1}/{self.N_Sampling} samplings \r")
            clean_stdout()

            SampledMutualInformationTable = np.array(
                SampledMutualInformationTable)
            AverageMutualInformationScoreList = np.mean(
                SampledMutualInformationTable, axis=0)

            VirusGroupingScheme_i += 1

            '''Sort the PPHMMs according to mutual information scores'''
            PPHMMOrder = list(list([zip(*sorted([(AverageMutualInformationScore_i, AverageMutualInformationScore) for AverageMutualInformationScore_i, AverageMutualInformationScore in enumerate(AverageMutualInformationScoreList)],
                                                key=lambda x: x[1],
                                                reverse=True,
                                                )
                                        )][0]
                                   )[0])

            ResultDict[VirusGroupingScheme]["AverageMutualInformationScoreList"] = [
                x for _, x in sorted(zip(PPHMMOrder, AverageMutualInformationScoreList))]
            ResultDict[VirusGroupingScheme]["PPHMMDesc"] = [
                x for _, x in sorted(zip(PPHMMOrder, IngroupVirus_PPHMMDesc))]
            ResultDict[VirusGroupingScheme]["PPHMMSignatureTable"] = [
                x for _, x in sorted(zip(PPHMMOrder, PPHMMSignatureTable_Subset))]

            '''Save results to numpy stack and .txt'''
            Dat = np.column_stack((list(map(", ".join, self.genomes["SeqIDLists"])),
                                   self.genomes["VirusNameList"],
                                   self.genomes["BaltimoreList"],
                                   self.genomes["OrderList"],
                                   self.genomes["FamilyList"],
                                   self.genomes["SubFamList"],
                                   self.genomes["GenusList"],
                                   self.genomes["TaxoGroupingList"],
                                   VirusGroupingList,
                                   PPHMMSignatureTable_Subset[:, PPHMMOrder],
                                   ))

            Dat = np.vstack((["Sequence ID", "Virus name", "Baltimore classification group", "Order", "Family", "Subfamily", "Genus", "TaxoGrouping", "Virus grouping"]+IngroupVirus_PPHMMDesc[PPHMMOrder].tolist(),
                             Dat,
                             ["Average mutual information score"]+[""]*8 +
                             AverageMutualInformationScoreList[PPHMMOrder].astype(
                                 str).tolist(),
                             ))

            MutualInformationScoreFile = f"{MutualInformationScoreDir}/MIScore.Scheme={VirusGroupingScheme}.txt"
            np.savetxt(fname=MutualInformationScoreFile,
                       X=Dat,
                       fmt='%s',
                       delimiter="\t")

            with open(MutualInformationScoreFile, "a") as MutualInformationScore_txt:
                MutualInformationScore_txt.write(
                    f"N_Sampling: {self.N_Sampling}\n Sampling strategy: {self.SamplingStrategy}\n Sample size per virus group: {self.SampleSizePerGroup}\n")

        '''Save all results to pickle'''
        if self.IncludeIncompleteGenomes == True:
            VariableShelveFile = self.VariableShelveDir + \
                "/MutualInformationCalculator.AllGenomes.p"
        elif self.IncludeIncompleteGenomes == False:
            VariableShelveFile = self.VariableShelveDir + \
                "/MutualInformationCalculator.CompleteGenomes.p"

        pickle.dump(ResultDict, open(VariableShelveFile, "wb"))

    def main(self):
        '''Compute mutual information scores'''
        section_header("#Compute mutual information score")

        '''1/3: Make fpaths'''
        MutualInformationScoreDir = self.mkdirs()

        '''2/3: Retrieve variables from VMR, ref annotations, PPHMMDB, virus grouping file'''
        self.genomes = retrieve_genome_vars(
            self.VariableShelveDir)
        self.ref_annotations = retrieve_pickle(self.VariableShelveDir) # TODO ref
        if "PPHMMSignatureTable_coo" in self.ref_annotations.keys():
            self.ref_annotations["PPHMMSignatureTable_coo"] = self.ref_annotations["PPHMMSignatureTable_coo"].toarray(
            )
        if "PPHMMLocationTable_coo" in self.ref_annotations.keys():
            self.ref_annotations["PPHMMLocationTable_coo"] = self.ref_annotations["PPHMMLocationTable_coo"].toarray(
            )

        self.pphmmdb_construction = retrieve_pickle(self.VariableShelveDir) # TODO pphmmdb
        PPHMMDesc = [
            "PPHMM|"+ClusterDesc for ClusterDesc in self.pphmmdb_construction["ClusterDescList"].astype("str")]
        PPHMMDesc = np.array(PPHMMDesc)

        if self.VirusGroupingFile == None:
            VirusGroupingDict = {"Overvall": self.genomes["TaxoGroupingList"]}
        else:
            VirusGroupingDict = self.read_virus_grouping_file()

        '''3/3: Calculate mutual information score between PPHMM and virus classification, save results'''
        self.calculate_mis(PPHMMDesc, VirusGroupingDict,
                           MutualInformationScoreDir)
