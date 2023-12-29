from Bio import Phylo, AlignIO, SeqIO
from io import StringIO
from collections import Counter
from copy import copy
from scipy.sparse import coo_matrix
from alive_progress import alive_it
import pandas as pd
import numpy as np
import subprocess
import os
import re
import string
import random
import glob
import pickle
import time

from app.utils.line_count import LineCount
from app.utils.ordered_set import OrderedSet
from app.utils.dist_mat_to_tree import DistMat2Tree
from app.utils.gomdb_constructor import GOMDB_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars
from app.utils.orf_identifier import find_orfs
from app.utils.error_handlers import raise_gravity_error
from app.utils.shell_cmds import shell
from app.utils.mkdirs import mkdir_ref_annotator
from app.utils.stdout_utils import progress_msg, warning_msg

class RefVirusAnnotator:
    def __init__(self,
                 GenomeSeqFile,
                 ShelveDir,
                 ProteinLength_Cutoff=100,
                 IncludeIncompleteGenomes=False,
                 HMMER_N_CPUs=20,
                 HMMER_C_EValue_Cutoff=1E-3,
                 HMMER_HitScore_Cutoff=0,
                 RemoveSingletonPPHMMs=False,
                 N_VirusesOfTheClassToIgnore=1,
                 PPHMMSorting=False,
                 HHsuite_evalue_Cutoff=1E-3,
                 HHsuite_pvalue_Cutoff=0.05,
                 HHsuite_N_CPUs=10,
                 HHsuite_QueryCoverage_Cutoff=75,
                 HHsuite_SubjectCoverage_Cutoff=75,
                 PPHMMClustering_MCLInflation=2,
                 ):
        '''Params'''
        self.GenomeSeqFile = GenomeSeqFile
        self.ExpDir = ShelveDir
        self.ExpOutputDir = f"{self.ExpDir}/output"
        self.pphmmdb_pickle = f"{self.ExpOutputDir}/PPHMMDBConstruction.p"
        self.IncludeIncompleteGenomes = IncludeIncompleteGenomes
        self.ProteinLength_Cutoff = ProteinLength_Cutoff
        self.HMMER_N_CPUs = HMMER_N_CPUs
        self.HMMER_C_EValue_Cutoff = HMMER_C_EValue_Cutoff
        self.HMMER_HitScore_Cutoff = HMMER_HitScore_Cutoff
        self.RemoveSingletonPPHMMs = RemoveSingletonPPHMMs
        self.N_VirusesOfTheClassToIgnore = N_VirusesOfTheClassToIgnore
        self.PPHMMSorting = PPHMMSorting
        self.HHsuite_evalue_Cutoff = HHsuite_evalue_Cutoff
        self.HHsuite_pvalue_Cutoff = HHsuite_pvalue_Cutoff
        self.HHsuite_N_CPUs = HHsuite_N_CPUs
        self.HHsuite_QueryCoverage_Cutoff = HHsuite_QueryCoverage_Cutoff
        self.HHsuite_SubjectCoverage_Cutoff = HHsuite_SubjectCoverage_Cutoff
        self.PPHMMClustering_MCLInflation = PPHMMClustering_MCLInflation
        '''Placehodler dirs & objects'''
        self.genomes = {}
        self.PPHMMSignatureTable, self.PPHMMLocationTable = [], []

    def PPHMMSignatureTable_Constructor(self, HMMER_hmmscanDir, HMMER_PPHMMDb):
        '''2/8: Create PPHMM signature and location tables'''
        progress_msg("- Generating PPHMM signature table and PPHMM location table")
        '''Load GenBank record'''
        # RM << TODO WHY DON'T WE USE THE SIG TABLE CONSTRUCTOR FROM UTILS??
        GenBankDict = SeqIO.index(self.GenomeSeqFile, "fasta" if os.path.splitext(self.GenomeSeqFile)[1] in [".fas", ".fst", ".fasta"] else "gb")
        Records_dict = {}
        for i in GenBankDict.items():
            ### RM < TEST: BACK TRANSCRIBE IF RNA SEQS USED
            if "u" in str(i[1].seq).lower():
                i[1].seq = i[1].seq.back_transcribe()
            Records_dict[i[0].split(".")[0]] = i[1]

        Records_dict = {k.split(".")[0]: v for k, v in Records_dict.items()}

        '''Specify PPHMMQueryFile and PPHMMScanOutFile'''
        PPHMMDB_Summary = f"{HMMER_PPHMMDb}_Summary.txt"
        N_PPHMMs = LineCount(PPHMMDB_Summary)-1
        PPHMMQueryFile = f"{HMMER_hmmscanDir}/QProtSeqs.fasta"
        PPHMMScanOutFile = f"{HMMER_hmmscanDir}/PPHMMScanOut.txt"
        PPHMMSignatureTable = np.empty((0, N_PPHMMs))
        PPHMMLocMiddleBestHitTable = np.empty((0, N_PPHMMs))

        for SeqIDList, TranslTable, BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping in alive_it(zip(self.genomes["SeqIDLists"], self.genomes["TranslTableList"], self.genomes["BaltimoreList"], self.genomes["OrderList"], self.genomes["FamilyList"], self.genomes["SubFamList"], self.genomes["GenusList"], self.genomes["VirusNameList"], self.genomes["TaxoGroupingList"])):
            GenBankSeqList, GenBankIDList, GenBankDescList = [], [], []
            for SeqID in SeqIDList:
                GenBankRecord = Records_dict[SeqID]
                GenBankSeqList.append(GenBankRecord.seq)
                GenBankIDList.append(GenBankRecord.id)
                GenBankDescList.append(GenBankRecord.description)

            '''sort lists by sequence/segment lengths'''
            _, GenBankSeqList, GenBankIDList, GenBankDescList = list(zip(
                *sorted(zip(list(map(len, list(map(str, GenBankSeqList)))), GenBankSeqList, GenBankIDList, GenBankDescList), reverse=True)))

            '''Get each orf for a genome'''
            ProtList, ProtIDList = [], []
            for GenBankSeq, GenBankID in zip(GenBankSeqList, GenBankIDList):
                prot, prot_id, _ = find_orfs(GenBankID, GenBankSeq, TranslTable, self.ProteinLength_Cutoff,
                                             taxonomy_annots=[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping])
                ProtList += prot
                ProtIDList += prot_id

            with open(PPHMMQueryFile, "w") as f:
                '''Write each translated ORF to the PPHMM Query File'''
                for i in range(len(ProtIDList)):
                    f.write(f">{ProtIDList[i]}\n{str(ProtList[i].seq)}\n")

            '''HMMSCAN query file vs PPHMM database'''
            shell(f"hmmscan --cpu {self.HMMER_N_CPUs} --noali --nobias --domtblout {PPHMMScanOutFile} {HMMER_PPHMMDb} {PPHMMQueryFile}",
                                "Ref virus annotator: signature table constructor, HMMScan call")

            PPHMMIDList, PPHMMScoreList, FeatureFrameBestHitList, FeatureLocFromBestHitList, \
                FeatureLocToBestHitList, FeatureDescList = [], [], [], [], [], []

            with open(PPHMMScanOutFile, "r") as PPHMMScanOut_txt:
                '''Open HMMSCAN results, collate scores for ORFs vs each PPHMM in database'''
                for Line in PPHMMScanOut_txt:
                    # RM < TODO Break out into fn as nesting here is upsetting me
                    if Line[0] == "#":
                        '''If header'''
                        continue
                    Line = Line.split()
                    '''Concatenate the cluster description back'''
                    Line[22] = " ".join(Line[22:])
                    Line = Line[:23]
                    C_EValue = float(Line[11])
                    HitScore = float(Line[7])
                    OriAASeqlen = float(len(GenBankSeq))/3

                    if C_EValue >= self.HMMER_C_EValue_Cutoff or HitScore <= self.HMMER_HitScore_Cutoff:
                        '''Threshold at user-set values'''
                        continue
                    '''Determine the frame and the location of the hit'''
                    iden = int(Line[0].split('_')[-1])
                    HitFrom = int(Line[17])
                    HitTo = int(Line[18])
                    HitMid = float(HitFrom+HitTo)/2
                    Frame = int(np.ceil(HitMid/OriAASeqlen)) if np.ceil(HitMid/OriAASeqlen) <= 3 else int(-(np.ceil(HitMid/OriAASeqlen)-3))
                    LocFrom = int(HitFrom % OriAASeqlen)
                    if LocFrom == 0:
                        '''If the hit occurs preciously from the end of the sequence'''
                        LocFrom = int(OriAASeqlen)
                    LocTo = int(HitTo % OriAASeqlen)
                    if LocTo == 0:
                        '''If the hit occurs preciously to the end of the sequence'''
                        LocTo = int(OriAASeqlen)
                    if LocTo < LocFrom:
                        '''The hit (falsely) spans across sequences of different frames'''
                        if np.ceil(HitFrom/OriAASeqlen) <= 3:
                            HitFrom_Frame = int(
                                np.ceil(HitFrom/OriAASeqlen))
                        else:
                            HitFrom_Frame = int(
                                -(np.ceil(HitFrom/OriAASeqlen)-3))

                        if np.ceil(HitTo/OriAASeqlen) <= 3:
                            HitTo_Frame = int(
                                np.ceil(HitTo/OriAASeqlen))
                        else:
                            HitTo_Frame = int(-(np.ceil(HitTo/OriAASeqlen)-3))

                        if Frame == HitFrom_Frame:
                            LocTo = int(OriAASeqlen)
                        elif Frame == HitTo_Frame:
                            LocFrom = int(1)
                        elif HitFrom_Frame != Frame and Frame != HitTo_Frame:
                            LocFrom = int(1)
                            LocTo = int(OriAASeqlen)
                        else:
                            warning_msg(
                                "Something is wrong with this PPHMMDB hit location determination")

                    if iden not in PPHMMIDList:
                        '''If new hit'''
                        Best_C_EValue = C_EValue
                        PPHMMIDList.append(iden)
                        PPHMMScoreList.append(HitScore)
                        FeatureDescList.append(
                            Line[22].split('|')[0])
                        FeatureFrameBestHitList.append(Frame)
                        FeatureLocFromBestHitList.append(LocFrom*3)
                        FeatureLocToBestHitList.append(LocTo*3)

                    elif iden in PPHMMIDList and C_EValue < Best_C_EValue:
                        '''Not new hit but score better than last'''
                        Best_C_EValue = C_EValue
                        FeatureFrameBestHitList[-1] = Frame
                        FeatureLocFromBestHitList[-1] = LocFrom*3
                        FeatureLocToBestHitList[-1] = LocTo*3

                    else:
                        '''Not new hit and score not as good as last'''
                        continue

            '''Absolute coordinate with orientation info encoded into it: +ve if the gene is present on the (+)strand, otherwise -ve'''
            FeatureLocMiddleBestHitList = np.zeros(N_PPHMMs)
            FeatureLocMiddleBestHitList[PPHMMIDList] = np.mean(np.array([FeatureLocFromBestHitList, FeatureLocToBestHitList]), axis=0)*(
                np.array(FeatureFrameBestHitList)/abs(np.array(FeatureFrameBestHitList)))
            PPHMMLocMiddleBestHitTable = np.vstack(
                (PPHMMLocMiddleBestHitTable, FeatureLocMiddleBestHitList))
            FeatureValueList = np.zeros(N_PPHMMs)
            FeatureValueList[PPHMMIDList] = PPHMMScoreList
            PPHMMSignatureTable = np.vstack(
                (PPHMMSignatureTable, FeatureValueList))

        '''Delete HMMER_hmmscanDir'''
        shell(f"rm -rf {HMMER_hmmscanDir}")
        return PPHMMSignatureTable, PPHMMLocMiddleBestHitTable

    def remove_singleton_pphmms(self, ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDb):
        '''Remove singleton PPHMMs from the PPHMM databases'''
        progress_msg("\tDetermining and removing PPHMMs ")

        '''Identify and remove non informative singleton PPHMMs from the database'''
        N_VirusesPerClass = Counter(self.genomes["TaxoGroupingList"])
        CandSingletonPPHMM_IndexList = np.where(
            np.sum(self.PPHMMSignatureTable > 0, axis=0) <= 1)[0]
        InfoSingletonPPHMM_IndexList = [PPHMM_i for PresentPPHMM_IndexList in [np.where(PPHMMSignature > 0)[0] for PPHMMSignature in self.PPHMMSignatureTable] if set(
            PresentPPHMM_IndexList).issubset(CandSingletonPPHMM_IndexList) for PPHMM_i in PresentPPHMM_IndexList]
        CandSingletonPPHMM_IndexList = [
            PPHMM_i for PPHMM_i in CandSingletonPPHMM_IndexList if PPHMM_i not in InfoSingletonPPHMM_IndexList]
        SingletonPPHMM_IndexList = []
        for PPHMM_i in CandSingletonPPHMM_IndexList:
            VirusWithTheSingletonPPHMM_i = np.where(
                self.PPHMMSignatureTable[:, PPHMM_i] != 0)[0]
            if len(VirusWithTheSingletonPPHMM_i) == 0:
                SingletonPPHMM_IndexList.append(PPHMM_i)
            else:
                if N_VirusesPerClass[self.genomes["TaxoGroupingList"][VirusWithTheSingletonPPHMM_i[0]]] > self.N_VirusesOfTheClassToIgnore:
                    SingletonPPHMM_IndexList.append(PPHMM_i)

        SelectedPPHMM_IndexList = np.array(
            list(range(self.PPHMMSignatureTable.shape[1])))
        SelectedPPHMM_IndexList = np.delete(
            arr=SelectedPPHMM_IndexList, obj=SingletonPPHMM_IndexList)

        '''Add ".ToBeDeleted" tag to each alignment and PPHMM file'''
        for f in glob.iglob(os.path.join(ClustersDir, '*.fasta')):
            os.rename(f, f + '.ToBeDeleted')

        for f in glob.iglob(os.path.join(HMMER_PPHMMDir, '*.hmm')):
            os.rename(f, f + '.ToBeDeleted')

        '''Rename and re-number informative alignment and PPHMM files'''
        PPHMM_j = 0
        for PPHMM_i in SelectedPPHMM_IndexList:
            ClusterFile_i = f"{ClustersDir}/Cluster_{PPHMM_i}.fasta.ToBeDeleted"
            ClusterFile_j = f"{ClustersDir}/Cluster_{PPHMM_j}.fasta"
            shell(f"mv {ClusterFile_i} {ClusterFile_j}")

            HMMER_PPHMMFile_j = f"{HMMER_PPHMMDir}/PPHMM_{PPHMM_j}.hmm"
            shell(f"mv {HMMER_PPHMMDir}/PPHMM_{PPHMM_i}.hmm.ToBeDeleted {HMMER_PPHMMFile_j}")

            '''Change the PPHMM name annotation'''
            with open(HMMER_PPHMMFile_j, "r+") as HMMER_PPHMM_txt:
                Contents = HMMER_PPHMM_txt.readlines()
                Contents[1] = "NAME  Cluster_%s\n" % PPHMM_j
                Contents = "".join(Contents)
                '''Put cursor at the beginning of the file'''
                HMMER_PPHMM_txt.seek(0)
                '''Write the contents, everything after the cursor'''
                HMMER_PPHMM_txt.write(Contents)
                HMMER_PPHMM_txt.truncate()

            PPHMM_j += 1

        '''Delete singleton clusters and PPHMMs (those with ".ToBeDeleted" tag) from the database'''
        for f in glob.glob(os.path.join(ClustersDir, "*.ToBeDeleted")):
            os.remove(f)

        for f in glob.glob(os.path.join(HMMER_PPHMMDir, "*.ToBeDeleted")):
            os.remove(f)

        '''Make a database of the remaining PPHMMs'''
        shell(f"find {HMMER_PPHMMDir} -name '*.hmm' -exec cat {{}} \; > {HMMER_PPHMMDb}")
        shell(f"hmmpress -f {HMMER_PPHMMDb}")

        'Remake the summary file of the new PPHMM database'
        HMMER_PPHMMDbSummaryFile = f"{HMMER_PPHMMDb}_Summary.txt"
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

        '''Reorganise the cluster meta data from PPHMMDBConstruction output'''
        parameters = pickle.load(open(self.pphmmdb_pickle, "rb"))

        '''Remove singleton PPHMMs' associated meta data, overwrite PPHMDBConstruction.p'''
        updated_parameters = {}
        updated_parameters["ClusterIDList"] = [
            f"Cluster_{Cluster_i}" for Cluster_i in range(len(SelectedPPHMM_IndexList))]
        updated_parameters["ClusterDescList"] = parameters["ClusterDescList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterSizeList"] = parameters["ClusterSizeList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterProtSeqIDList"] = parameters["ClusterProtSeqIDList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterSizeByTaxoGroupingList"] = parameters["ClusterSizeByTaxoGroupingList"][SelectedPPHMM_IndexList]
        updated_parameters["ClusterSizeByProtList"] = parameters["ClusterSizeByProtList"][SelectedPPHMM_IndexList]
        pickle.dump(updated_parameters, open(self.pphmmdb_pickle, "wb"))

    def sort_pphmm(self, ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDb, HHsuite_PPHMMDir, HHsuite_PPHMMDB, HHsuiteDir):
        '''Sort PPHMMs by similarity order, via AVA comparison. Overwrite db.'''
        progress_msg("Sorting PPHMMs")

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
        for PPHMM_i in alive_it(range(N_PPHMMs)):
            AlnClusterFile = f"{ClustersDir}/Cluster_{PPHMM_i}.fasta"
            Aln = AlignIO.read(AlnClusterFile, "fasta")
            N_Seqs = len(Aln)
            L_Seqs = Aln.get_alignment_length()
            HHsuite_PPHMMFile = f"{HHsuite_PPHMMDir}/PPHMM_{PPHMM_i}.hhm"
            shell(f"hhmake -i {AlnClusterFile} -o {HHsuite_PPHMMFile} -seq {N_Seqs+1} -name Cluster_{PPHMM_i} -id 100 -M 50 -v 0")

        '''Rebuild the HHsuite PPHMM database'''
        wdir = f'{"/".join(HHsuite_PPHMMDB.split("/")[:-2])}/HHSuiteDB/'
        fname = f"{wdir}/mycluster"
        shell(f"mkdir {wdir}")
        shell(f"rm {HHsuite_PPHMMDB}_hhm.ffdata {HHsuite_PPHMMDB}_hhm.ffindex",
                    "PPHMMDB Contruction: rebuild hhsite db, remove existing db")
        shell(f"ffindex_build -s {fname}_msa.ffdata {fname}_msa.index {ClustersDir}")
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
        hhsearchDir = f'{HHsuiteDir}/hhsearch_{"".join(random.choice(string.ascii_uppercase + string.digits)for _ in range(10))}'
        os.makedirs(hhsearchDir)
        hhsearchOutFile = f"{hhsearchDir}/hhsearch.stdout.hhr"
        N_PPHMMs = LineCount(f"{fname}_hhm.ffindex")

        SeenPair, SeenPair_i, PPHMMSimScoreCondensedMat = {}, 0, []
        hit_match = re.compile(
            r"[A-Z0-9]+\|[A-Z0-9]+.[0-9]{1}\s[a-zA-Z0-9_ ]{0,10}")
        for PPHMM_i in range(0, N_PPHMMs):
            HHsuite_PPHMMFile = f"{HHsuite_PPHMMDir}/PPHMM_{PPHMM_i}.hhm"
            shell(f"hhsearch -maxmem 9.0 -i {HHsuite_PPHMMFile} -d {fname} -o {hhsearchOutFile} -e {self.HHsuite_evalue_Cutoff} -E {self.HHsuite_evalue_Cutoff} -cov {self.HHsuite_SubjectCoverage_Cutoff} -z 1 -b 1 -id 100 -global -v 0 -cpu {self.HHsuite_N_CPUs}")

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
                        if (evalue <= self.HHsuite_evalue_Cutoff and pvalue <= self.HHsuite_pvalue_Cutoff and qcovs >= self.HHsuite_QueryCoverage_Cutoff and scovs >= self.HHsuite_SubjectCoverage_Cutoff):
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
        shell(f"mcl {PPHMMSimScoreCondensedMatFile} --abc -o {PPHMM_ClusterFile} -I {self.PPHMMClustering_MCLInflation}")

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
        shell(f"rm -rf {HHsuiteDir}")

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
        for f in glob.iglob(os.path.join(ClustersDir, '*.fasta')):
            os.rename(f, f + '.tmp')

        for f in glob.iglob(os.path.join(HMMER_PPHMMDir, '*.hmm')):
            os.rename(f, f + '.tmp')

        '''Rename and re-number alignment and PPHMM files'''
        PPHMM_j = 0
        for PPHMM_i in PPHMMOrder:
            ClusterFile_i = f"{ClustersDir}/Cluster_{PPHMM_i}.fasta.tmp"
            ClusterFile_j = f"{ClustersDir}/Cluster_{PPHMM_j}.fasta"
            shell(f"mv {ClusterFile_i} {ClusterFile_j}")

            HMMER_PPHMMFile_j = f"{HMMER_PPHMMDir}/PPHMM_{PPHMM_j}.hmm"
            shell(f"mv {HMMER_PPHMMDir}/PPHMM_{PPHMM_i}.hmm.tmp {HMMER_PPHMMFile_j}")

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
        progress_msg("Building ordered PPHMM DB")
        shell(f"find {HMMER_PPHMMDir} -name '*.hmm' -exec cat {{}} \; > {HMMER_PPHMMDb}")
        shell(f"hmmpress -f {HMMER_PPHMMDb}")

        '''Remake the summary file of the new PPHMM database'''
        HMMER_PPHMMDbSummaryFile = f"{HMMER_PPHMMDb}_Summary.txt"
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

        progress_msg("\tSort PPHMMSignatureTable and PPHMMLocationTable")
        self.PPHMMSignatureTable = self.PPHMMSignatureTable[:, PPHMMOrder]
        self.PPHMMLocationTable = self.PPHMMLocationTable[:, PPHMMOrder]

        '''Load cluster meta data from PPHMMDBConstruction pickle, reorganise'''
        parameters = pickle.load(open(self.pphmmdb_pickle, "rb"))

        '''Sort PPHMMs' associated meta data, overwrite PPHMMDBConstruction pickle'''
        updated_parameters = {}
        updated_parameters["ClusterIDList"] = [
            f"Cluster_{Cluster_i}" for Cluster_i in range(len(PPHMMOrder))]
        updated_parameters["ClusterDescList"] = parameters["ClusterDescList"][PPHMMOrder]
        updated_parameters["ClusterSizeList"] = parameters["ClusterSizeList"][PPHMMOrder]
        updated_parameters["ClusterProtSeqIDList"] = parameters["ClusterProtSeqIDList"][PPHMMOrder]
        updated_parameters["ClusterSizeByTaxoGroupingList"] = parameters["ClusterSizeByTaxoGroupingList"][PPHMMOrder]
        updated_parameters["ClusterSizeByProtList"] = parameters["ClusterSizeByProtList"][PPHMMOrder]
        pickle.dump(updated_parameters, open(self.pphmmdb_pickle, "wb"))

    def save_sig_tables(self, GOMIDList, GOMSignatureTable, GOMDB):
        '''8/8: Load in PPHMMDB, to extract cluster list...'''
        parameters = pickle.load(open(self.pphmmdb_pickle, "rb"))

        '''Combine with original data from VMR'''
        PPHMMDesc = [
            f"PPHMM|{ClusterDesc}" for ClusterDesc in parameters["ClusterDescList"].astype("str")]
        GOMDesc = [f"GOM|{TaxoGrouping}" for TaxoGrouping in GOMIDList]
        header = ["Virus name", "Accession", "Order", "Family", "Subfamily", "Genus", "Class"] + PPHMMDesc + GOMDesc
        table_data = np.column_stack((self.genomes["VirusNameList"],
                                      list(map(", ".join,
                                               self.genomes["SeqIDLists"])),
                                      self.genomes["OrderList"],
                                      self.genomes["FamilyList"],
                                      self.genomes["SubFamList"],
                                      self.genomes["GenusList"],
                                      self.genomes["TaxoGroupingList"],
                                      self.PPHMMSignatureTable,
                                      GOMSignatureTable))
        out_df = pd.DataFrame(table_data, index=None, columns=header)
        out_df.to_csv(f"{self.ExpOutputDir}/PPHMMandGOMsignatures.csv")

        '''Save'''
        parameters["PPHMMSignatureTable"] = self.PPHMMSignatureTable
        parameters["GOMSignatureTable"] = GOMSignatureTable
        parameters["PPHMMLocationTable"] = self.PPHMMLocationTable
        parameters["GOMIDList"] = GOMIDList
        parameters["PPHMMSignatureTable_coo"] = coo_matrix(
            self.PPHMMSignatureTable)
        parameters["PPHMMLocationTable_coo"] = coo_matrix(
            self.PPHMMLocationTable)
        parameters["GOMDB_coo"] = {GOMID: coo_matrix(
            GOM) for GOMID, GOM in GOMDB.items()}

        pickle.dump(parameters, open(f"{self.ExpOutputDir}/RefVirusAnnotator.p", "wb"))

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

        '''1/8 : Make directories'''
        HMMER_hmmscanDir, ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDb, HHsuite_PPHMMDB, HHsuite_PPHMMDir, HHsuite_PPHMMDB, HHsuiteDir = mkdir_ref_annotator(self.ExpDir, self.RemoveSingletonPPHMMs, self.PPHMMSorting)

        '''2/8 : Retrieve variables'''
        self.genomes = retrieve_genome_vars(
            self.ExpOutputDir)

        '''3/8 : Generate PPHMMSignatureTable and PPHMMLocationTable'''
        self.PPHMMSignatureTable, self.PPHMMLocationTable = self.PPHMMSignatureTable_Constructor(
            HMMER_hmmscanDir, HMMER_PPHMMDb)

        if self.RemoveSingletonPPHMMs:
            '''4/8 : (OPT) Remove singleton PPHMMs'''
            self.remove_singleton_pphmms(
                ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDb)

        if self.PPHMMSorting == True:
            '''5/8 : (OPT) Sort PPHMMs'''
            self.sort_pphmm(ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDb,
                            HHsuite_PPHMMDir, HHsuite_PPHMMDB, HHsuiteDir)

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
