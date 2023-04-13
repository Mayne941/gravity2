import numpy as np
from ete3 import Tree
from sklearn.svm import SVC

from app.utils.ordered_set import OrderedSet
from app.utils.stdout_utils import progress_bar, clean_stdout


def PairwiseSimilarityScore_Cutoff_Dict_Constructor(SimMat, TaxoGroupingList, N_PairwiseSimilarityScores) -> dict:
    '''Accessory to classify (6): Use machine learning classifier to construct similarity score cutoff dict'''
    N_Viruses = len(SimMat)
    PairwiseSimilarityScore_Cutoff_Dict = {TaxoGrouping: {"PairwiseSimilarityScore_InterClass_List": [
    ], "PairwiseSimilarityScore_IntraClass_List": []} for TaxoGrouping in OrderedSet(TaxoGroupingList)}
    for Virus_i in range(N_Viruses):
        for Virus_j in range(Virus_i, N_Viruses):
            Class_Virus_i = TaxoGroupingList[Virus_i]
            Class_Virus_j = TaxoGroupingList[Virus_j]
            if Class_Virus_i == Class_Virus_j:
                PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_i]["PairwiseSimilarityScore_IntraClass_List"].append(
                    SimMat[Virus_i][Virus_j])
            else:
                PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_i]["PairwiseSimilarityScore_InterClass_List"].append(
                    SimMat[Virus_i][Virus_j])
                PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_j]["PairwiseSimilarityScore_InterClass_List"].append(
                    SimMat[Virus_i][Virus_j])

    for TaxoGrouping in OrderedSet(TaxoGroupingList):
        PairwiseSimilarityScore_IntraClass_List = PairwiseSimilarityScore_Cutoff_Dict[
            TaxoGrouping]["PairwiseSimilarityScore_IntraClass_List"]
        PairwiseSimilarityScore_InterClass_List = PairwiseSimilarityScore_Cutoff_Dict[
            TaxoGrouping]["PairwiseSimilarityScore_InterClass_List"]

        if len(PairwiseSimilarityScore_IntraClass_List) > N_PairwiseSimilarityScores:
            PairwiseSimilarityScore_IntraClass_List = np.random.choice(
                PairwiseSimilarityScore_IntraClass_List, N_PairwiseSimilarityScores, replace=False).tolist()
            PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping][
                "PairwiseSimilarityScore_IntraClass_List"] = PairwiseSimilarityScore_IntraClass_List

        if len(PairwiseSimilarityScore_InterClass_List) > N_PairwiseSimilarityScores:
            PairwiseSimilarityScore_InterClass_List = np.random.choice(
                PairwiseSimilarityScore_InterClass_List, N_PairwiseSimilarityScores, replace=False).tolist()
            PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping][
                "PairwiseSimilarityScore_InterClass_List"] = PairwiseSimilarityScore_InterClass_List

        PairwiseSimilarityList = np.array(
            PairwiseSimilarityScore_IntraClass_List+PairwiseSimilarityScore_InterClass_List).reshape(-1, 1)
        Labels = np.array(["Intra"]*len(PairwiseSimilarityScore_IntraClass_List) +
                          ["Inter"]*len(PairwiseSimilarityScore_InterClass_List))

        # RM <<<<<<<<<<<< SVC -- NEEDS REPLACING, TRAINING AND TUNING
        clf = SVC(kernel="linear", class_weight="balanced")
        if np.unique(Labels).shape[0] == 1:
            # RM < Set control to outside of main cluster if only 1 detected
            Labels[-1] = "Extra"
        clf			.fit(PairwiseSimilarityList, Labels)
        PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["CutOff"] = - \
            clf.intercept_[0]/clf.coef_[0][0]

    return PairwiseSimilarityScore_Cutoff_Dict


def TaxonomicAssignmentProposerAndEvaluator(SimMat_UcfVirusesVSRefViruses, TaxoGrouping, VirusDendrogram, TaxoLabelList_RefVirus, TaxoLabelList_UcfVirus, PairwiseSimilarityScore_Cutoff_Dict):
    '''Accessory to classify (6): Propose a taxonomic class to each unclassified viruses, using the 1-nearest neighbour algorithm'''
    MaxSimScoreList = np.max(SimMat_UcfVirusesVSRefViruses, axis=1)
    TaxoOfMaxSimScoreList = TaxoGrouping[np.argmax(
        SimMat_UcfVirusesVSRefViruses, axis=1)]

    vd_temp = VirusDendrogram ####################
    print("Evaluating taxonomic assignments")
    '''Rename the reference virus leaves in the dendrogram to class label'''
    VirusDendrogram = Tree(VirusDendrogram)
    for LeafNode in VirusDendrogram.get_leaves():
        if LeafNode.name in TaxoLabelList_RefVirus:
            LeafNode.name = TaxoGrouping[np.where(
                np.array(TaxoLabelList_RefVirus) == LeafNode.name)[0]][0]

    TaxoAssignmentList, PhyloStatList, N_UcfViruses = [
    ], [], len(SimMat_UcfVirusesVSRefViruses)
    for UcfVirus_i in range(N_UcfViruses):
        CandidateTaxoAssignment = TaxoOfMaxSimScoreList[UcfVirus_i]
        MaxSimScore = MaxSimScoreList[UcfVirus_i]

        original_label = TaxoLabelList_UcfVirus[UcfVirus_i][-8:]
        '''1st criterion: see if the similarity score between the unclassified virus and the best match is greater than the similarity cutoff of the proposed class'''
        if MaxSimScore > PairwiseSimilarityScore_Cutoff_Dict[CandidateTaxoAssignment]["CutOff"]:
            '''if pass the 1st criterion Prune the dendrogram so that it contains only reference sequences and the virus of interest (VoI)'''
            VirusDendrogram_tmp = VirusDendrogram.copy()
            try:
                for q in np.delete(TaxoLabelList_UcfVirus, UcfVirus_i):
                    VirusDendrogram_tmp.get_leaves_by_name(q)[0].delete()
            except Exception as ex:
                print(
                    f"Error on Taxonomic Assignment: could not find dendrogram node to delete: {ex}")

            try:
                '''Get the leaf names (taxonomic assignments) of the sister clade of the unclassified virus'''
                VoI = TaxoLabelList_UcfVirus[UcfVirus_i]
                VoISisterClade = VirusDendrogram_tmp.search_nodes(name=VoI)[
                    0].get_sisters()[0]
                VoISisterClade_Leafnames = OrderedSet(
                    VoISisterClade.get_leaf_names())
                if VoISisterClade.is_leaf():
                    VoISisterClade_Leafnames1 = OrderedSet(
                        VoISisterClade.get_leaf_names())
                    VoISisterClade_Leafnames2 = []
                else:
                    VoISisterClade_Leafnames1 = OrderedSet(
                        VoISisterClade.get_children()[0].get_leaf_names())
                    VoISisterClade_Leafnames2 = OrderedSet(
                        VoISisterClade.get_children()[1].get_leaf_names())

                '''Get the leaf names (taxonomic assignments) of the immediate out group of the unclassified virus'''
                ImmediateAncestorNode = VirusDendrogram_tmp.search_nodes(name=VoI)[
                    0].up
                if len(ImmediateAncestorNode.get_sisters()) != 0:
                    ImmediateOutCladeLeafnames = OrderedSet(
                        ImmediateAncestorNode.get_sisters()[0].get_leaf_names())
                else:
                    ImmediateOutCladeLeafnames = []

                '''2nd criterion: see if the candidate taxonomic assignment is supported by the dendrogram...'''
                if len(VoISisterClade_Leafnames) == 1 and [CandidateTaxoAssignment] == VoISisterClade_Leafnames:
                    '''If the sister clade consists entirely of the candidate class...'''
                    if VoISisterClade_Leafnames == ImmediateOutCladeLeafnames:
                        '''and the out group is the same, then assign VoI to the candidate class'''
                        TaxoAssignmentList.append(CandidateTaxoAssignment)
                        '''Embedded within a clade of a single class '''
                        PhyloStatList.append("1")

                    else:
                        '''... but the out group isn't, having a sister relationship with the proposed class and they are similar enough (1st criterion)'''
                        TaxoAssignmentList.append(CandidateTaxoAssignment)
                        PhyloStatList.append("2")

                elif len(ImmediateOutCladeLeafnames) == 1 and [CandidateTaxoAssignment] == ImmediateOutCladeLeafnames:
                    '''If the immediate outgroup consists entirely of the candidate class...'''
                    if(VoISisterClade_Leafnames1 == [CandidateTaxoAssignment] or VoISisterClade_Leafnames2 == [CandidateTaxoAssignment]):
                        '''.. and the VoI is sandwished between 2 branches of the candidate class...'''
                        TaxoAssignmentList.append(CandidateTaxoAssignment)
                        '''Sandwiched between 2 branches of the candidate class'''
                        PhyloStatList.append("3")
                    else:
                        '''else the candidate class is accepted on the ground that it has a paraphyletic relationship with the candidate class (just inside)'''
                        TaxoAssignmentList.append(CandidateTaxoAssignment)
                        '''Having a paraphyletic relationship with the candidate class (just inside)'''
                        PhyloStatList.append("4")

                elif (VoISisterClade_Leafnames1 == [CandidateTaxoAssignment] or VoISisterClade_Leafnames2 == [CandidateTaxoAssignment]):
                    '''If one of the two branches in the sister clade consists entirely of the candidate class...'''
                    TaxoAssignmentList.append(CandidateTaxoAssignment)
                    '''Having a paraphyletic relationship with the candidate class (just outside)'''
                    PhyloStatList.append("5")
                else:
                    TaxoAssignmentList.append(
                        f"Unclassified - NSBD: {original_label}")
                    '''The candidate class is not supported by the dendrogram'''
                    PhyloStatList.append("6")
            except:
                continue
        else:
            '''The unclassified virus isn't similar enough to the members of the candidate class'''
            TaxoAssignmentList.append(f"Unclassified - NS: {original_label}")
            PhyloStatList.append("NA")

        progress_bar(
            f"\033[K Taxonomic assignment evaluation: [{'='*int(float(UcfVirus_i+1)/N_UcfViruses*20)}] {UcfVirus_i+1}/{N_UcfViruses} viruses \r")
    clean_stdout()

    return MaxSimScoreList, TaxoOfMaxSimScoreList, TaxoAssignmentList, PhyloStatList
