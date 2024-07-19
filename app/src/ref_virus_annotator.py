from Bio import Phylo, AlignIO
from io import StringIO
from collections import Counter
from copy import copy
from scipy.sparse import coo_matrix
from alive_progress import alive_it
import pandas as pd
import numpy as np
import os
import re
import string
import random
import glob
import pickle

from app.utils.line_count import LineCount
from app.utils.ordered_set import OrderedSet
from app.utils.dist_mat_to_tree import DistMat2Tree
from app.utils.gomdb_constructor import GOMDB_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_pickle
from app.utils.error_handlers import raise_gravity_error
from app.utils.shell_cmds import shell
from app.utils.mkdirs import mkdir_ref_annotator
from app.utils.stdout_utils import progress_msg
from app.utils.generate_fnames import generate_file_names
#
from app.utils.parallel_sig_generator import PPHMMSignatureTable_Constructor
# from app.utils.pphmm_signature_table_constructor import PPHMMSignatureTable_Constructor_DEPRECATED as PPHMMSignatureTable_Constructor #########

class RefVirusAnnotator:
    def __init__(self,
                 payload,
                 GenomeSeqFile,
                 ExpDir,
                 ):
        '''Params'''
        self.payload = payload
        self.fnames = generate_file_names(payload,ExpDir)
        self.GenomeSeqFile = GenomeSeqFile
        self.genomes = retrieve_genome_vars(self.fnames["ReadGenomeDescTablePickle"])
        mkdir_ref_annotator(self.fnames, self.payload['PPHMMSorting'])
        '''Placehodler dirs & objects'''
        self.PPHMMSignatureTable, self.PPHMMLocationTable, self.NaivePPHMMLocationTable = [], [], []

    def remove_singleton_pphmms(self):
        '''Remove singleton PPHMMs from the PPHMM databases'''
        progress_msg(f"- Determining and removing PPHMMs shared by less than {self.payload['N_VirusesOfTheClassToIgnore']} genomes")

        # TODO < RM BETTER SINGLETON INDEX FINDER
        sig_table = pd.DataFrame(self.PPHMMSignatureTable)
        sig_table_bin = sig_table.applymap(lambda x: 1 if x > 0 else 0)
        pphmm_usage_counts = sig_table_bin.apply(lambda x: sum(x)).tolist()
        SingletonPPHMM_IndexList = []
        for idx, pphmm_count in enumerate(pphmm_usage_counts):
            if pphmm_count == 0 or pphmm_count < self.payload['N_VirusesOfTheClassToIgnore']:
                SingletonPPHMM_IndexList.append(idx)
        #

        '''Identify and remove non informative singleton PPHMMs from the database'''
        # N_VirusesPerClass = Counter(self.genomes["TaxoGroupingList"])
        # CandSingletonPPHMM_IndexList = np.where(
        #     np.sum(self.PPHMMSignatureTable > 0, axis=0) <= 1)[0]
        # InfoSingletonPPHMM_IndexList = [PPHMM_i for PresentPPHMM_IndexList in [np.where(PPHMMSignature > 0)[0] for PPHMMSignature in self.PPHMMSignatureTable] if set(
        #     PresentPPHMM_IndexList).issubset(CandSingletonPPHMM_IndexList) for PPHMM_i in PresentPPHMM_IndexList]
        # CandSingletonPPHMM_IndexList = [
        #     PPHMM_i for PPHMM_i in CandSingletonPPHMM_IndexList if PPHMM_i not in InfoSingletonPPHMM_IndexList]
        # SingletonPPHMM_IndexList = []
        # for PPHMM_i in CandSingletonPPHMM_IndexList:
        #     VirusWithTheSingletonPPHMM_i = np.where(
        #         self.PPHMMSignatureTable[:, PPHMM_i] != 0)[0]
        #     if len(VirusWithTheSingletonPPHMM_i) == 0:
        #         '''If zero hits''' # TODO Why would there be zero hits?
        #         SingletonPPHMM_IndexList.append(PPHMM_i)
        #     else:
        #         if N_VirusesPerClass[self.genomes["TaxoGroupingList"][VirusWithTheSingletonPPHMM_i[0]]] > self.payload['N_VirusesOfTheClassToIgnore']:
        #             SingletonPPHMM_IndexList.append(PPHMM_i)

        SelectedPPHMM_IndexList = np.array(
            list(range(self.PPHMMSignatureTable.shape[1])))
        SelectedPPHMM_IndexList = np.delete(
            arr=SelectedPPHMM_IndexList, obj=SingletonPPHMM_IndexList)

        '''Add ".ToBeDeleted" tag to each alignment and PPHMM file'''
        for f in glob.iglob(os.path.join(self.fnames['ClustersDir'], '*.fasta')):
            os.rename(f, f + '.ToBeDeleted')

        for f in glob.iglob(os.path.join(self.fnames['HMMER_PPHMMDir'], '*.hmm')):
            os.rename(f, f + '.ToBeDeleted')

        '''Rename and re-number informative alignment and PPHMM files'''
        PPHMM_j = 0
        for PPHMM_i in SelectedPPHMM_IndexList:
            ClusterFile_j = f"{self.fnames['ClustersDir']}/Cluster_{PPHMM_j}.fasta"
            ClusterFile_i = f"{ClusterFile_j}.ToBeDeleted"
            shell(f"mv {ClusterFile_i} {ClusterFile_j}")

            HMMER_PPHMMFile_j = f"{self.fnames['HMMER_PPHMMDir']}/PPHMM_{PPHMM_j}.hmm"
            shell(f"mv {self.fnames['HMMER_PPHMMDir']}/PPHMM_{PPHMM_i}.hmm.ToBeDeleted {HMMER_PPHMMFile_j}")

            '''Change the PPHMM name annotation'''
            with open(HMMER_PPHMMFile_j, "r+") as HMMER_PPHMM_txt:
                Contents = HMMER_PPHMM_txt.readlines()
                Contents[1] = f"NAME  Cluster_{PPHMM_j}\n"
                Contents = "".join(Contents)
                '''Put cursor at the beginning of the file'''
                HMMER_PPHMM_txt.seek(0)
                '''Write the contents, everything after the cursor'''
                HMMER_PPHMM_txt.write(Contents)
                HMMER_PPHMM_txt.truncate()

            PPHMM_j += 1

        '''Delete singleton clusters and PPHMMs (those with ".ToBeDeleted" tag) from the database'''
        for f in glob.glob(os.path.join(self.fnames['ClustersDir'], "*.ToBeDeleted")):
            os.remove(f)

        for f in glob.glob(os.path.join(self.fnames['HMMER_PPHMMDir'], "*.ToBeDeleted")):
            os.remove(f)

        '''Make a database of the remaining PPHMMs'''
        shell(f"find {self.fnames['HMMER_PPHMMDir']} -name '*.hmm' -exec cat {{}} \; > {self.fnames['HMMER_PPHMMDb']}")
        shell(f"hmmpress -f {self.fnames['HMMER_PPHMMDb']}")

        'Remake the summary file of the new PPHMM database'
        HMMER_PPHMMDbSummaryFile = f"{self.fnames['HMMER_PPHMMDb']}_Summary.txt"
        Contents, PPHMM_i = [], 0
        with open(HMMER_PPHMMDbSummaryFile, "r") as HMMER_PPHMMDbSummary_txt:
            Contents.append(next(HMMER_PPHMMDbSummary_txt))
            for Line in HMMER_PPHMMDbSummary_txt:
                Line = Line.split("\t")
                if int(Line[0].split("_")[-1]) in SelectedPPHMM_IndexList:
                    Line[0] = "Cluster_%s" % PPHMM_i
                    Contents.append("\t".join(Line))
                    PPHMM_i = PPHMM_i + 1

        Contents = "".join(Contents)
        with open(HMMER_PPHMMDbSummaryFile, "w") as HMMER_PPHMMDbSummary_txt:
            HMMER_PPHMMDbSummary_txt.write(Contents)

        '''Remove singleton PPHMMs from PPHMMSignatureTable and PPHMMLocationTable'''
        self.PPHMMSignatureTable = np.delete(
            arr=self.PPHMMSignatureTable, obj=SingletonPPHMM_IndexList, axis=1)
        self.PPHMMLocationTable = np.delete(
            arr=self.PPHMMLocationTable, obj=SingletonPPHMM_IndexList, axis=1)
        self.NaivePPHMMLocationTable = np.delete( # RM < TODO Fix 30/04/24
            arr=self.NaivePPHMMLocationTable, obj=SingletonPPHMM_IndexList, axis=1)
        '''Reorganise the cluster meta data from PPHMMDBConstruction output'''
        parameters = retrieve_pickle(self.fnames['PphmmdbPickle'])

        '''Remove singleton PPHMMs' associated meta data, overwrite PPHMDBConstruction.p'''
        updated_parameters = {}
        updated_parameters["ClusterIDList"] = [
            f"Cluster_{Cluster_i}" for Cluster_i in range(len(SelectedPPHMM_IndexList))]
        updated_parameters["ClusterDescList"] = parameters["ClusterDescList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterSizeList"] = parameters["ClusterSizeList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterProtSeqIDList"] = parameters["ClusterProtSeqIDList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterSizeByTaxoGroupingList"] = parameters["ClusterSizeByTaxoGroupingList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterSizeByProtList"] = parameters["ClusterSizeByProtList"][SelectedPPHMM_IndexList]
        pickle.dump(updated_parameters, open(self.fnames['PphmmdbPickle'], "wb"))

    def sort_pphmm(self):
        '''Sort PPHMMs by similarity order, via AVA comparison. Overwrite db.'''
        progress_msg("- Sorting PPHMMs in order of similarity and rewriting PPHMM database. This may take a while.")

        '''Determine the PPHMM order by 'virus profile' similarity'''
        TraitValueTable = copy(self.PPHMMSignatureTable.transpose())
        N_PPHMMs = len(TraitValueTable)
        TaxoLabelList = list(range(N_PPHMMs))

        TraitValueTable[TraitValueTable !=
                        0], SimMat_traits, PPHMM_i = 1, [], 1
        for TraitValues in alive_it(TraitValueTable):
            TraitValuesMat = np.tile(TraitValues, (N_PPHMMs, 1))
            SimMat_traits	.append(np.sum(np.minimum(TraitValuesMat, TraitValueTable),
                                  axis=1)/np.sum(np.maximum(TraitValuesMat, TraitValueTable), axis=1))
            PPHMM_i += 1

        '''Constructe a PPHMM dendrogram, and extract the PPHMM order'''
        SimMat_traits = np.array(SimMat_traits)
        TreeNewick_traits = DistMat2Tree(DistMat=1 - SimMat_traits,
                                         LeafList=TaxoLabelList,
                                         Dendrogram_LinkageMethod="average")

        TreeNewick_traits = Phylo.read(StringIO(TreeNewick_traits), "newick")
        TreeNewick_traits.ladderize(reverse=True)
        PPHMMOrder_ByTree = np.array(
            [int(Clade.name) for Clade in TreeNewick_traits.get_terminals()])

        '''Cluster PPHMMs by PPHMM similiarty, make dbs'''
        N_PPHMMs = self.PPHMMSignatureTable.shape[1]
        progress_msg("\t- Clustering PPHMMs by similarity.")
        for PPHMM_i in alive_it(range(N_PPHMMs)):
            AlnClusterFile = f"{self.fnames['ClustersDir']}/Cluster_{PPHMM_i}.fasta"
            Aln = AlignIO.read(AlnClusterFile, "fasta")
            N_Seqs = len(Aln)
            L_Seqs = Aln.get_alignment_length()
            HHsuite_PPHMMFile = f"{self.fnames['HHsuite_PPHMMDir']}/PPHMM_{PPHMM_i}.hhm"
            shell(f"hhmake -i {AlnClusterFile} -o {HHsuite_PPHMMFile} -seq {N_Seqs+1} -name Cluster_{PPHMM_i} -id 100 -M 50 -v 0")

        '''Rebuild the HHsuite PPHMM database'''
        wdir = f'{"/".join(self.fnames["HHsuite_PPHMMDB"].split("/")[:-2])}/HHSuiteDB/'
        fname = f"{wdir}/mycluster"
        shell(f"mkdir {wdir}")
        shell(f"rm {self.fnames['HHsuite_PPHMMDB']}_hhm.ffdata {self.fnames['HHsuite_PPHMMDB']}_hhm.ffindex",
                    "PPHMMDB Contruction: rebuild hhsite db, remove existing db")
        shell(f"ffindex_build -s {fname}_msa.ffdata {fname}_msa.index {self.fnames['ClustersDir']}")
        shell(f"ffindex_apply {fname}_msa.ffdata {fname}_msa.index -i {fname}_a3m.ffindex -d {fname}_a3m.ffdata -- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0")
        shell(f"rm {fname}_msa.ffdata {fname}_msa.index")
        shell(f"ffindex_apply {fname}_a3m.ffdata {fname}_a3m.ffindex -i {fname}_hhm.ffindex -d {fname}_hhm.ffdata -- hhmake -i stdin -o stdout -v 0")
        shell(f"cstranslate -f -x 0.3 -c 4 -I a3m -i {fname}_a3m -o {fname}_cs219")
        shell(f"sort -k3 -n -r {fname}_cs219.ffindex | cut -f1 > sorting.dat") # Parameterise sorting file
        shell(f"ffindex_order sorting.dat {fname}_hhm.ff{{data,index}} {fname}_hhm_ordered.ffdata {fname}_hhm_ordered.ffindex")
        shell(f"mv {fname}_hhm_ordered.ffindex {fname}_hhm.ffindex")
        shell(f"mv {fname}_hhm_ordered.ffdata {fname}_hhm.ffdata")
        shell(f"ffindex_order sorting.dat {fname}_a3m.ffdata {fname}_a3m.ffindex {fname}_a3m_ordered.ffdata {fname}_a3m_ordered.ffindex")
        shell(f"mv {fname}_a3m_ordered.ffindex {fname}_a3m.ffindex")
        shell(f"mv {fname}_a3m_ordered.ffdata {fname}_a3m.ffdata")

        '''Determine PPHMM-PPHMM similarity (AVA hhsearch)'''
        hhsearchDir = f'{self.fnames["HHsuiteDir"]}/hhsearch_{"".join(random.choice(string.ascii_uppercase + string.digits)for _ in range(10))}'
        os.makedirs(hhsearchDir)
        hhsearchOutFile = f"{hhsearchDir}/hhsearch.stdout.hhr"
        N_PPHMMs = LineCount(f"{fname}_hhm.ffindex")

        SeenPair, SeenPair_i, PPHMMSimScoreCondensedMat = {}, 0, []
        hit_match = re.compile(
            r"[A-Z0-9]+\|[A-Z0-9]+.[0-9]{1}\s[a-zA-Z0-9_ ]{0,10}")

        progress_msg("\t- Regenerating protein profile scores.")
        for PPHMM_i in alive_it(range(0, N_PPHMMs)):
            HHsuite_PPHMMFile = f"{self.fnames['HHsuite_PPHMMDir']}/PPHMM_{PPHMM_i}.hhm"
            shell(f"hhsearch -maxmem 9.0 -i {HHsuite_PPHMMFile} -d {fname} -o {hhsearchOutFile} -e {self.payload['HHsuite_evalue_Cutoff']} -E {self.payload['HHsuite_evalue_Cutoff']} -cov {self.payload['HHsuite_SubjectCoverage_Cutoff']} -z 1 -b 1 -id 100 -global -v 0 -cpu {self.payload['N_CPUs']}")
            # breakpoint()

            '''Open output from hhsearch'''
            with open(hhsearchOutFile, 'r') as hhsearchOut_txt:
                Content = hhsearchOut_txt.readlines()
                QueryLength = int(Content[1].split()[1])
                '''Iterate over best hits and extract data line by line, append to sim score matrix'''
                for Line in Content[9:]:
                    if Line == "\n":
                        break
                    else:
                        '''Filter variable length Hit name'''
                        Line = re.sub(hit_match, "", Line)
                        Line = Line.replace("(", " ").replace(")", " ").split()
                        try:
                            PPHMM_j = int(Line[0])
                            evalue = float(Line[2])
                            pvalue = float(Line[3])
                            PPHMMSimScore = float(Line[4])
                            Col = float(Line[6])
                            SubjectLength = int(Line[-1])
                        except ValueError as ex:
                            raise_gravity_error(f"Failed to extract PPHMM data from HHsearch output (ref virus annotator, sort PPHMMs function). This happens when HHsearch outputs malformed text files, check the output for the entry listed in this exception: {ex}")
                        qcovs = Col/QueryLength*100
                        scovs = Col/SubjectLength*100
                        if (evalue <= self.payload['HHsuite_evalue_Cutoff'] and pvalue <= self.payload['HHsuite_pvalue_Cutoff'] and qcovs >= self.payload['HHsuite_QueryCoverage_Cutoff'] and scovs >= self.payload['HHsuite_SubjectCoverage_Cutoff']):
                            Pair = ", ".join(
                                sorted(map(str, [PPHMM_i, PPHMM_j])))
                            if Pair in SeenPair:
                                '''If the pair has already been seen...'''
                                if PPHMMSimScore > PPHMMSimScoreCondensedMat[SeenPair[Pair]][2]:
                                    '''and if the new PPHMMSimScore is higher...'''
                                    PPHMMSimScoreCondensedMat[SeenPair[Pair]
                                                              ][2] = PPHMMSimScore
                            else:
                                SeenPair[Pair] = SeenPair_i
                                PPHMMSimScoreCondensedMat.append(
                                    [PPHMM_i, PPHMM_j, PPHMMSimScore])
                                SeenPair_i += 1

            shell(f"rm {hhsearchOutFile}")

        '''Structure and make similarity score matrix'''
        PPHMMSimScoreCondensedMat = np.array(PPHMMSimScoreCondensedMat)
        PPHMMSimScoreCondensedMat = np.array(
            [PPHMMSimScorePair for PPHMMSimScorePair in PPHMMSimScoreCondensedMat if PPHMMSimScorePair[0] < PPHMMSimScorePair[1]])
        PPHMMSimScoreCondensedMatFile = f"{hhsearchDir}/PPHMMSimScoreCondensedMat.txt"
        np.savetxt(fname=PPHMMSimScoreCondensedMatFile,
                   X=PPHMMSimScoreCondensedMat,
                   fmt='%s',
                   delimiter="\t",
                   header="PPHMM_i\tPPHMM_j\tPPHMMSimScore")

        '''Cluster PPHMMs based on hhsearch scores, using the MCL algorithm'''
        PPHMM_ClusterFile = f"{hhsearchDir}/PPHMM_Clusters.txt"
        shell(f"mcl {PPHMMSimScoreCondensedMatFile} --abc -o {PPHMM_ClusterFile} -I 2") # RM < Hard code for now
        # shell(f"mcl {PPHMMSimScoreCondensedMatFile} --abc -o {PPHMM_ClusterFile} -I {self.payload['PPHMMClustering_MCLInflation']}")

        SeenPPHMMIDList = []
        with open(PPHMM_ClusterFile, 'r') as PPHMM_Cluster_txt:
            for Cluster in PPHMM_Cluster_txt.readlines():
                SeenPPHMMIDList.extend(Cluster.split("\n")[0].split("\t"))

        with open(PPHMM_ClusterFile, 'a') as PPHMM_Cluster_txt:
            PPHMM_Cluster_txt.write("\n".join(
                list(set([str(float(x)) for x in range(N_PPHMMs)])-set(SeenPPHMMIDList))))

        with open(PPHMM_ClusterFile, 'r') as PPHMM_Cluster_txt:
            PPHMMOrder_ByMCL = PPHMM_Cluster_txt.readlines()

        PPHMMOrder_ByMCL = [list(map(int, list(map(float, Cluster.split(
            "\n")[0].split("\t"))))) for Cluster in PPHMMOrder_ByMCL]

        '''Delete the hhsuite shelve directory and database'''
        shell(f"rm -rf {self.fnames['HHsuiteDir']}")

        '''Determine the final PPHMM order'''
        PPHMMOrder, PPHMMClusterSeparation_IndexList = [], []
        PPHMMOrder_ByTree_tmp = copy(PPHMMOrder_ByTree)
        PPHMMOrder_ByMCL_tmp = copy(PPHMMOrder_ByMCL)
        while len(PPHMMOrder_ByTree_tmp) != 0:
            PPHMM_i = PPHMMOrder_ByTree_tmp[0]
            Cluster_i = np.where(
                [PPHMM_i in cluster for cluster in PPHMMOrder_ByMCL_tmp])[0][0]
            Cluster = PPHMMOrder_ByMCL_tmp.pop(Cluster_i)
            i, PPHMM_i = list(zip(
                *[(i, PPHMM_i) for i, PPHMM_i in enumerate(PPHMMOrder_ByTree_tmp) if PPHMM_i in Cluster]))
            PPHMMOrder.extend(PPHMM_i)
            PPHMMClusterSeparation_IndexList.append(len(PPHMMOrder))
            PPHMMOrder_ByTree_tmp = np.delete(PPHMMOrder_ByTree_tmp, i)

        '''Reorganise the PPHMM database'''
        '''Add ".tmp" tag to each alignment and PPHMM file'''
        for f in glob.iglob(os.path.join(self.fnames['ClustersDir'], '*.fasta')):
            os.rename(f, f + '.tmp')

        for f in glob.iglob(os.path.join(self.fnames['HMMER_PPHMMDir'], '*.hmm')):
            os.rename(f, f + '.tmp')

        '''Rename and re-number alignment and PPHMM files'''
        PPHMM_j = 0
        for PPHMM_i in PPHMMOrder:
            ClusterFile_j = f"{self.fnames['ClustersDir']}/Cluster_{PPHMM_j}.fasta"
            ClusterFile_i = f"{ClusterFile_j}.tmp"
            shell(f"mv {ClusterFile_i} {ClusterFile_j}")

            HMMER_PPHMMFile_j = f"{self.fnames['HMMER_PPHMMDir']}/PPHMM_{PPHMM_j}.hmm"
            shell(f"mv {self.fnames['HMMER_PPHMMDir']}/PPHMM_{PPHMM_i}.hmm.tmp {HMMER_PPHMMFile_j}")

            '''Change the PPHMM name annotation'''
            with open(HMMER_PPHMMFile_j, "r+") as HMMER_PPHMM_txt:
                Contents = HMMER_PPHMM_txt.readlines()
                Contents[1] = f"NAME  Cluster_{PPHMM_j}\n"
                Contents = "".join(Contents)
                HMMER_PPHMM_txt.seek(0)
                HMMER_PPHMM_txt.write(Contents)
                HMMER_PPHMM_txt.truncate()

            PPHMM_j += 1

        '''Make a new (ordered) PPHMM database'''
        progress_msg("\t - Building ordered PPHMM DB")
        shell(f"find {self.fnames['HMMER_PPHMMDir']} -name '*.hmm' -exec cat {{}} \; > {self.fnames['HMMER_PPHMMDb']}")
        shell(f"hmmpress -f {self.fnames['HMMER_PPHMMDb']}")

        '''Remake the summary file of the new PPHMM database'''
        HMMER_PPHMMDbSummaryFile = f"{self.fnames['HMMER_PPHMMDb']}_Summary.txt"
        with open(HMMER_PPHMMDbSummaryFile, "r") as HMMER_PPHMMDbSummary_txt:
            Contents = HMMER_PPHMMDbSummary_txt.readlines()
            Header = Contents.pop(0)[2:]
            Contents = list(np.array(Contents)[PPHMMOrder])
            PPHMM_i = 0
            for Content in Contents:
                Content = Content.split("\t")
                Content[0] = f"Cluster_{PPHMM_i}"
                Contents[PPHMM_i] = "\t".join(Content)
                PPHMM_i += 1
            Contents = [Header] + Contents

        Contents = "".join(Contents)
        with open(HMMER_PPHMMDbSummaryFile, "w") as HMMER_PPHMMDbSummary_txt:
            HMMER_PPHMMDbSummary_txt.write(Contents)

        progress_msg("\t- Sorting new PPHMMSignatureTable and PPHMMLocationTable")
        self.PPHMMSignatureTable = self.PPHMMSignatureTable[:, PPHMMOrder]
        self.PPHMMLocationTable = self.PPHMMLocationTable[:, PPHMMOrder]

        '''Load cluster meta data from PPHMMDBConstruction pickle, reorganise'''
        parameters = retrieve_pickle(self.fnames['PphmmdbPickle'])

        '''Sort PPHMMs' associated meta data, overwrite PPHMMDBConstruction pickle'''
        updated_parameters = {}
        updated_parameters["ClusterIDList"] = [
            f"Cluster_{Cluster_i}" for Cluster_i in range(len(PPHMMOrder))]
        updated_parameters["ClusterDescList"] = parameters["ClusterDescList"][PPHMMOrder]
        updated_parameters["ClusterSizeList"] = parameters["ClusterSizeList"][PPHMMOrder]
        updated_parameters["ClusterProtSeqIDList"] = parameters["ClusterProtSeqIDList"][PPHMMOrder]
        updated_parameters["ClusterSizeByTaxoGroupingList"] = parameters["ClusterSizeByTaxoGroupingList"][PPHMMOrder]
        updated_parameters["ClusterSizeByProtList"] = parameters["ClusterSizeByProtList"][PPHMMOrder]
        pickle.dump(updated_parameters, open(self.fnames['PphmmdbPickle'], "wb"))

    def save_sig_tables(self, GOMIDList, GOMSignatureTable, GOMDB):
        '''8/8: Load in PPHMMDB, to extract cluster list...'''
        parameters = retrieve_pickle(self.fnames['PphmmdbPickle'])

        '''Combine with original data from VMR, save as CSV for shared PPHMM heatmaps'''
        PPHMMDesc = [
            f"PPHMM|{ClusterDesc}" for ClusterDesc in parameters["ClusterDescList"].astype("str")]
        header = ["Virus name"] + PPHMMDesc
        table_data = np.column_stack((self.genomes["VirusNameList"],
                                      self.PPHMMSignatureTable))
        out_df = pd.DataFrame(table_data, index=None, columns=header)
        out_df.to_csv(self.fnames["PphmmAndGomSigs"])

        '''Save'''
        parameters["PPHMMSignatureTable"] = self.PPHMMSignatureTable
        parameters["GOMSignatureTable"] = GOMSignatureTable
        parameters["PPHMMLocationTable"] = self.PPHMMLocationTable
        parameters["NaivePPHMMLocationTable"] = self.NaivePPHMMLocationTable
        parameters["GOMIDList"] = GOMIDList
        parameters["PPHMMSignatureTable_coo"] = coo_matrix(
            self.PPHMMSignatureTable)
        parameters["PPHMMLocationTable_coo"] = coo_matrix(
            self.PPHMMLocationTable)
        parameters["GOMDB_coo"] = {GOMID: coo_matrix(
            GOM) for GOMID, GOM in GOMDB.items()}

        pickle.dump(parameters, open(self.fnames['RefAnnotatorPickle'], "wb"))

    def main(self):
        '''
        Generate PPHMM signature table, PPHMM location table, GOM database, and GOM signature table for
        reference viruses, using their own PPHMM database. When RemoveSingletonPPHMMs == TRUE, singleton
        PPHMMs will be removed if the virus that has it belongs to a group that contains more than
        "N_VirusesOfTheClassToIgnore" members; i.e. groups with less than 'N_VirusesOfTheClassToIgnore'
        members will be ignored.
        '''
        section_header(
            "Generate PPHMM signature & loc tables, GOM database & GOM sig table")

        '''3/8 : Generate PPHMMSignatureTable and PPHMMLocationTable'''
        self.PPHMMSignatureTable, self.PPHMMLocationTable, self.NaivePPHMMLocationTable = PPHMMSignatureTable_Constructor(
                    self.genomes,
                    self.payload,
                    self.fnames,
                    self.GenomeSeqFile,
                    self.fnames['HMMER_PPHMMDb']
        )

        if self.payload['RemoveSingletonPPHMMs']:
            '''4/8 : (OPT) Remove singleton PPHMMs'''
            self.remove_singleton_pphmms()

        if self.payload['PPHMMSorting'] == True:
            '''5/8 : (OPT) Sort PPHMMs'''
            self.sort_pphmm()

        '''6/8 : Make GOM database'''
        GOMIDList = OrderedSet([TaxoGrouping for TaxoGrouping in self.genomes["TaxoGroupingList"].astype(
            'str') if not TaxoGrouping.startswith(("_", "*"))])
        GOMDB = GOMDB_Constructor(
            self.genomes["TaxoGroupingList"], self.PPHMMLocationTable, GOMIDList)

        progress_msg("- Generate GOM signature table")
        '''7/8 : Make GOM signature table'''
        GOMSignatureTable = GOMSignatureTable_Constructor(
            self.PPHMMLocationTable, GOMDB, GOMIDList)

        '''8/8 : Save PPHMMSignatureTable and GOMSignatureTable'''
        self.save_sig_tables(GOMIDList, GOMSignatureTable, GOMDB)
