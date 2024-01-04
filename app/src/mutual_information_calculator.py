from sklearn.feature_selection import mutual_info_classif
from alive_progress import alive_it
import numpy as np
import pickle
import pandas as pd

from app.utils.ordered_set import OrderedSet
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_pickle
from app.utils.stdout_utils import progress_msg
from app.utils.generate_fnames import generate_file_names
from app.utils.mkdirs import mkdir_mi_scorer

class MutualInformationCalculator:
    def __init__(self,
                 payload,
                 ExpDir,
                 ):
        self.payload = payload,
        if type(self.payload) == tuple: self.payload = self.payload[0] # RM < TODO Fix this weirdness
        self.fnames = generate_file_names(payload,ExpDir)
        self.genomes = retrieve_genome_vars(self.fnames['ReadGenomeDescTablePickle'])
        self.pphmmdb_construction = retrieve_pickle(self.fnames['PphmmdbPickle'])
        self.ref_annotations = retrieve_pickle(self.fnames['RefAnnotatorPickle'])
        mkdir_mi_scorer(self.fnames)

    def get_pphmm_clusters(self):
        '''1/3: Retrieve PPHMM clusters and descriptions from PPHMMDB Constructor pickle'''
        return np.array(["PPHMM|"+ClusterDesc for ClusterDesc in self.pphmmdb_construction["ClusterDescList"].astype("str")])

    def read_virus_grouping_file(self) -> dict:
        '''2/3: (Opt) Read virus grouping file to dict'''
        progress_msg("Loading virus groupings from GRAViTy grouping file")
        VirusGroupingTable = []
        with open(self.fnames['VirusGroupingFile'], "r") as VirusGrouping_txt:
            VirusGroupingSchemeList = next(VirusGrouping_txt)
            VirusGroupingSchemeList = VirusGroupingSchemeList.split("\r\n")[0].split("\n")[
                0].split("\t")
            for Line in VirusGrouping_txt:
                Line = Line.split("\r\n")[0].split("\n")[0].split("\t")
                VirusGroupingTable.append(Line)
        '''Group at GRAViTy-assigned family level'''
        return {"overall": np.array([i[1] for i in VirusGroupingTable if len(i) > 1])}

    def calculate_mis(self, PPHMMDesc, VirusGroupingDict) -> None:
        '''3/3: Calculate mutual info scores (sklearn), output to np to .txt, then save pickle'''
        progress_msg("Calculating Mutual Information Score (this can take a while...)")
        ResultDict = {}
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
            for _ in alive_it(range(self.payload['N_Sampling'])):
                if self.payload['SamplingStrategy'] == None:
                    SampledVirus_IndexList = np.where(
                        [VirusGroupingID != "-" for VirusGroupingID in VirusGroupingList])[0]

                elif self.payload['SamplingStrategy'] == "balance_without_repeat":
                    SampledVirus_IndexList = sum([list(np.random.permutation(np.where(VirusGroupingList == VirusGroupingID)[0])[
                                                 :self.payload['SampleSizePerGroup']]) for VirusGroupingID in VirusGroupingIDList], [])

                elif self.payload['SamplingStrategy'] == "balance_with_repeat":
                    SampledVirus_IndexList = sum([list(np.random.choice(a=np.where(VirusGroupingList == VirusGroupingID)[
                                                 0], size=self.payload['SampleSizePerGroup'], replace=True)) for VirusGroupingID in VirusGroupingIDList], [])

                '''Call SK-Learn mutual info classifier. This component scales badly but no replacement currently exists.'''
                SampledMutualInformationTable.append(mutual_info_classif(
                    PPHMMSignatureTable_Subset[SampledVirus_IndexList], VirusGroupingList[SampledVirus_IndexList]))

            SampledMutualInformationTable = np.array(
                SampledMutualInformationTable)
            AverageMutualInformationScoreList = np.mean(
                SampledMutualInformationTable, axis=0)


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
            '''Save'''
            df = pd.DataFrame(Dat)
            df.columns = df.iloc[0]
            df = df.iloc[1:]
            df.to_csv(self.fnames['MutualInformationScoreFile'])

        '''Save all results to pickle'''
        pickle.dump(ResultDict, open(self.fnames['MiScorePickle'], "wb"))

    def main(self):
        '''Compute mutual information scores'''
        section_header("Compute mutual information score")

        '''1/3: Load PPHMM data'''
        progress_msg("Loading & grouping PPHMM data")
        PPHMMDesc = self.get_pphmm_clusters()

        '''2/3: Group PPHMMs'''
        if self.payload["VirusGrouping"]:
            VirusGroupingDict = self.read_virus_grouping_file()
        else:
            VirusGroupingDict = {"overall": self.genomes["TaxoGroupingList"]}

        '''3/3: Calculate mutual information score between PPHMM and virus classification, save results'''
        self.calculate_mis(PPHMMDesc, VirusGroupingDict)
