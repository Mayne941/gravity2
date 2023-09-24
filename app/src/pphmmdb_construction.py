from ete3 import Tree
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from alive_progress import alive_it
import numpy as np
import os
import operator
import random
import string
import re
import pickle
import pandas as pd

from app.utils.line_count import LineCount
from app.utils.dist_mat_to_tree import DistMat2Tree
from app.utils.download_genbank_file import DownloadGenBankFile
from app.utils.console_messages import section_header
from app.utils.orf_identifier import get_orf_trasl_table, no_orf_match
from app.utils.stdout_utils import clean_stdout, progress_bar, error_handler, progress_msg, warning_msg
from app.utils.retrieve_pickle import retrieve_genome_vars
from app.utils.shell_cmds import shell
from app.utils.mkdirs import mkdir_pphmmdbc
from app.utils.hhr_parse import hhparse
from app.utils.error_handlers import raise_gravity_error
from app.utils.orf_identifier import find_orfs

class PPHMMDBConstruction:
    def __init__(self,
                 genbank_email,
                 GenomeSeqFile,
                 ShelveDir,
                 ProteinLength_Cutoff=100,
                 IncludeIncompleteGenomes=True,
                 BLASTp_evalue_Cutoff=1E-3,
                 BLASTp_PercentageIden_Cutoff=50,
                 BLASTp_QueryCoverage_Cutoff=75,
                 BLASTp_SubjectCoverage_Cutoff=75,
                 BLASTp_num_alignments=1000000,
                 BLASTp_N_CPUs=20,
                 MUSCLE_GapOpenCost=-3.0,
                 MUSCLE_GapExtendCost=-0.0,
                 ProtClustering_MCLInflation=2,
                 N_AlignmentMerging=0,
                 HHsuite_evalue_Cutoff=1E-6,
                 HHsuite_pvalue_Cutoff=0.05,
                 HHsuite_N_CPUs=10,
                 HHsuite_QueryCoverage_Cutoff=85,
                 HHsuite_SubjectCoverage_Cutoff=85,
                 PPHMMClustering_MCLInflation=5,
                 HMMER_PPHMMDb_ForEachRoundOfPPHMMMerging=True,
                 ) -> None:
        self.GenomeSeqFile = GenomeSeqFile
        self.ShelveDir = ShelveDir
        self.ProteinLength_Cutoff = ProteinLength_Cutoff
        self.IncludeIncompleteGenomes = IncludeIncompleteGenomes
        self.BLASTp_evalue_Cutoff = BLASTp_evalue_Cutoff
        self.BLASTp_PercentageIden_Cutoff = BLASTp_PercentageIden_Cutoff
        self.BLASTp_QueryCoverage_Cutoff = BLASTp_QueryCoverage_Cutoff
        self.BLASTp_SubjectCoverage_Cutoff = BLASTp_SubjectCoverage_Cutoff
        self.BLASTp_num_alignments = BLASTp_num_alignments
        self.BLASTp_N_CPUs = BLASTp_N_CPUs
        self.MUSCLE_GapOpenCost = MUSCLE_GapOpenCost
        self.MUSCLE_GapExtendCost = MUSCLE_GapExtendCost
        self.ProtClustering_MCLInflation = ProtClustering_MCLInflation
        self.N_AlignmentMerging = N_AlignmentMerging
        self.HHsuite_evalue_Cutoff = HHsuite_evalue_Cutoff
        self.HHsuite_pvalue_Cutoff = HHsuite_pvalue_Cutoff
        self.HHsuite_N_CPUs = HHsuite_N_CPUs
        self.HHsuite_QueryCoverage_Cutoff = HHsuite_QueryCoverage_Cutoff
        self.HHsuite_SubjectCoverage_Cutoff = HHsuite_SubjectCoverage_Cutoff
        self.PPHMMClustering_MCLInflation = PPHMMClustering_MCLInflation
        self.HMMER_PPHMMDb_ForEachRoundOfPPHMMMerging = HMMER_PPHMMDb_ForEachRoundOfPPHMMMerging
        self.genbank_email = genbank_email
        self.orf_tranl_table = get_orf_trasl_table()
        self.orf_no_match = no_orf_match()
        self.BLASTp_outfmt = '"6 qseqid sseqid pident qcovs qlen slen evalue bitscore"'

    def get_genbank(self, genomes):
        '''3/10: Check if GenBank files exist; dl if needed. Index & transform.'''
        if not os.path.isfile(self.GenomeSeqFile):
            progress_msg("Download GenBank file\n GenomeSeqFile doesn't exist. GRAViTy is downloading the GenBank file(s).")
            DownloadGenBankFile(self.GenomeSeqFile,
                                genomes["SeqIDLists"], self.genbank_email)

        progress_msg("Reading GenBank file")
        GenBankDict = SeqIO.index(self.GenomeSeqFile, "genbank")
        return {k.split(".")[0]: v for k, v in GenBankDict.items()}

    def sequence_extraction(self, genomes, GenBankDict):
        '''4/10: Extract protein sequences from VMR using GenBank data; manually annotate ORFs if needed'''
        progress_msg(
            f"- Extract/predict protein sequences from virus genomes, excluding proteins with lengthes <{self.ProteinLength_Cutoff} aa")
        ProtList, ProtIDList, Virus_i = [
        ], [], 1

        for SeqIDList, TranslTable, BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping in alive_it(zip(genomes["SeqIDLists"], genomes["TranslTableList"], genomes["BaltimoreList"], genomes["OrderList"], genomes["FamilyList"], genomes["SubFamList"], genomes["GenusList"], genomes["VirusNameList"], genomes["TaxoGroupingList"])):
            for SeqID in SeqIDList:
                '''Sometimes an Acc ID doesn't have a matching record (usually when multiple seqs for 1 virus)... - skip if true'''
                try:
                    GenBankRecord = GenBankDict[SeqID]
                except KeyError as ex:
                    raise SystemExit(f"{ex}\n"
                                        f"ERROR: I couldn't extract a sequence ID from the input Genbank file.\n"
                                        f"This usually happens when you've tried to re-run the experiment with an old .gb file.\n"
                                        f"Try deleting your input .gb file ({self.GenomeSeqFile}), then starting again.")
                GenBankID = GenBankRecord.name
                GenBankFeatures = GenBankRecord.features

                '''Extract protein sequences'''
                ContainProtAnnotation = 0
                for Feature in GenBankFeatures:
                    if(Feature.type == 'CDS' and "protein_id" in Feature.qualifiers and "translation" in Feature.qualifiers):
                        ContainProtAnnotation = 1
                        try:
                            ProtName = Feature.qualifiers["product"][0]
                        except KeyError:
                            try:
                                ProtName = Feature.qualifiers["gene"][0]
                            except KeyError:
                                try:
                                    ProtName = Feature.qualifiers["note"][0]
                                except KeyError:
                                    ProtName = "Hypothetical protein"
                        ProtID = Feature.qualifiers["protein_id"][0]
                        ProtSeq = Feature.qualifiers["translation"][0]
                        if len(ProtSeq) >= self.ProteinLength_Cutoff:
                            ProtRecord = SeqRecord(Seq(ProtSeq),
                                                   id=GenBankID+"|"+ProtID,
                                                   name=GenBankID+"|"+ProtID,
                                                   description=ProtName,
                                                   annotations={'taxonomy': [BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping]})
                            ProtList.append(ProtRecord)
                            ProtIDList.append(GenBankID+"|"+ProtID)

                '''If the genome isn't annotated with any ORFs, find some'''
                if ContainProtAnnotation == 0:
                    prots, prot_ids = find_orfs(GenBankID, GenBankRecord.seq, TranslTable, self.ProteinLength_Cutoff,
                                                    taxonomy_annots=[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping])
                    ProtList += prots
                    ProtIDList += prot_ids
            Virus_i += 1
        clean_stdout()
        return ProtList, np.array(ProtIDList)

    def make_blastp_db(self, BLASTSubjectFile, ProtList):
        '''5/10: Make BLAST DB'''
        progress_msg("Building BLASTp DB")
        with open(BLASTSubjectFile, "w") as BLASTSubject_txt:
            SeqIO.write(ProtList, BLASTSubject_txt, "fasta")
        shell(f"makeblastdb -in {BLASTSubjectFile} -dbtype prot", "PPHMMDB Construction: make BLASTp db")

    def blastp_analysis(self, ProtList, BLASTQueryFile, BLASTSubjectFile, BLASTOutputFile, BLASTBitScoreFile):
        '''6/10: Perform ALL-VERSUS-ALL BLASTp analysis'''
        progress_msg("Performing ALL-VERSUS-ALL BLASTp analysis")
        BitScoreMat, SeenPair, SeenPair_i, N_ProtSeqs = [
            ], {}, 0, len(ProtList)
        for ProtSeq_i in alive_it(range(N_ProtSeqs)):
            '''BLAST query fasta file'''
            BLASTQuery = ProtList[ProtSeq_i]
            with open(BLASTQueryFile, "w") as BLASTQuery_txt: _ = SeqIO.write(BLASTQuery, BLASTQuery_txt, "fasta")

            '''Perform BLASTp, load output to dataframe'''
            shell(f'blastp -query {BLASTQueryFile} -db {BLASTSubjectFile} -out {BLASTOutputFile} -evalue {self.BLASTp_evalue_Cutoff} -outfmt {self.BLASTp_outfmt} -num_alignments {self.BLASTp_num_alignments} -num_threads {self.BLASTp_N_CPUs}',
                  f"PPHMMDB Construction: BLASTp (protein ID = {ProtList[ProtSeq_i].id})")

            try:
                blast_df = pd.read_csv(BLASTOutputFile, sep="\t", names=["qseqid", "sseqid", "pident", "qcovs", "qlen", "slen", "evalue", "bitscore"])
            except Exception as e:
                raise FileNotFoundError(f"Could not open BLASTp output with exception: {e}")

            if blast_df.empty:
                '''If no hits in sequence, continue'''
                continue

            '''Load BLASTp results, collate bit scores for matches'''
            blast_iter = blast_df.to_dict(orient="records")
            for i in blast_iter:
                pair = ", ".join(sorted([i["qseqid"], i["sseqid"]]))
                if ((i["qseqid"] != i["sseqid"]) and (i["pident"] >= self.BLASTp_PercentageIden_Cutoff)
                    and (i["qcovs"] >= self.BLASTp_QueryCoverage_Cutoff)
                    and ((i["qcovs"]*i["qlen"]/i["slen"]) >= self.BLASTp_SubjectCoverage_Cutoff)):
                    '''Query must: not match subject, have identity > thresh, have query coverage > thresh and query coverage normalised to subject length > thresh'''
                    if pair in SeenPair:
                        '''If the pair has already been seen...'''
                        if i["bitscore"] > BitScoreMat[SeenPair[pair]][2]:
                            '''...and if the new bitscore is higher'''
                            BitScoreMat[SeenPair[pair]][2] = i["bitscore"]
                    else:
                        SeenPair[pair] = SeenPair_i
                        BitScoreMat.append(
                            [pair.split(", ")[0], pair.split(", ")[1], i["bitscore"]])
                        SeenPair_i += 1
        clean_stdout()

        BitScoreMat = np.array(BitScoreMat)
        '''Save protein-protein similarity scores (BLASTp bit scores)'''
        np.savetxt(fname=BLASTBitScoreFile,
                   X=BitScoreMat,
                   fmt='%s',
                   delimiter="\t",
                   header="SeqID_I\tSeqID_II\tBit score"
                   )

    def mcl_clustering(self, ProtIDList, BLASTBitScoreFile, BLASTProtClusterFile):
        '''7/10: Use Mcl to do clustering on BLASTp bit scores'''
        progress_msg("- Doing protein sequence clustering based on BLASTp bit scores, using the MCL algorithm")
        shell(f"mcl {BLASTBitScoreFile} --abc -o {BLASTProtClusterFile} -I {self.ProtClustering_MCLInflation}",
                    "PPHMMDB Contruction: mcl clustering, call mcl")

        '''For each cluster found, pull out seen proteins by ID'''
        SeenProtIDList = []
        with open(BLASTProtClusterFile, 'r') as BLASTProtCluster_txt:
            for Cluster in BLASTProtCluster_txt.readlines():
                SeenProtIDList.extend(Cluster.split(
                    "\r\n")[0].split("\n")[0].split("\t"))

        '''Append seen protein IDs to clusters'''
        with open(BLASTProtClusterFile, 'a') as BLASTProtCluster_txt:
            BLASTProtCluster_txt.write(
                "\n".join(list(set(ProtIDList)-set(SeenProtIDList))))

    def make_alignments(self, ProtList, ProtIDList, BLASTProtClusterFile, ClustersDir):
        '''8/10: Do protein alignments with muscle align, make cluster alignment annotations'''
        progress_msg("- Make protein alignments")
        _, Cluster_i, Cluster_MetaDataDict = LineCount(
            BLASTProtClusterFile)+1, 0, {}
        with open(BLASTProtClusterFile, 'r') as BLASTProtCluster_txt:
            for Cluster in alive_it(BLASTProtCluster_txt.readlines()):
                HitList, TaxoLists, DescList, Cluster = [], [], [], Cluster.split("\n")[
                    0].split("\t")
                for ProtID in Cluster:
                    HitList.append(
                        ProtList[np.where(ProtIDList == ProtID)[0][0]])
                    TaxoLists.append(HitList[-1].annotations['taxonomy'])
                    DescList.append(HitList[-1].description.replace(", ", " ").replace(",", " ").replace(": ", "_").replace(
                        ":", "_").replace("; ", " ").replace(";", " ").replace(" (", "/").replace("(", "/").replace(")", ""))

                '''Cluster file'''
                AlnClusterFile = f"{ClustersDir}/Cluster_{Cluster_i}.fasta"
                with open(AlnClusterFile, "w") as UnAlnClusterTXT:
                    p = SeqIO.write(HitList, UnAlnClusterTXT, "fasta")

                temp_aln_fname = f"{ClustersDir}/temp.fasta"
                '''Align cluster using muscle'''
                shell(f"mafft --thread {self.HHsuite_N_CPUs} --localpair --maxiterate 1000  {AlnClusterFile} > {temp_aln_fname}",
                        "PPHMMDB Construction: make alignments, main muscle call")

                shell(f"rm {AlnClusterFile} && mv {temp_aln_fname} {AlnClusterFile}",
                      "PPHMMDB COnstruction: move temp mafft file")
                '''Cluster annotations'''
                Cluster_MetaDataDict[Cluster_i] = {"Cluster": Cluster,
                                                "DescList": DescList,
                                                "TaxoLists": TaxoLists,
                                                "AlignmentLength": AlignIO.read(AlnClusterFile, "fasta").get_alignment_length()
                                                }
                Cluster_i += 1
        clean_stdout()
        return Cluster_MetaDataDict

    def pphmm_and_merge_alignments(self, Cluster_MetaDataDict, ClustersDir, HHsuite_PPHMMDir,
                                   HMMER_PPHMMDir, HMMER_PPHMMDbDir, HHsuite_PPHMMDB, HHsuiteDir,
                                   VariableShelveDir):
        '''9/10 (OPTIONAL): Make HHSuite PPHMMs and DB, merge & make protein alignmens'''
        if self.N_AlignmentMerging > 0:
            progress_msg(
                f"- Merge protein alignments, {self.N_AlignmentMerging} rounds of merging")
        elif self.N_AlignmentMerging < 0:
            progress_msg("- Merge protein alignments until exhausted")

        AlignmentMerging_i_round = 0
        wdir = f'{"/".join(HHsuite_PPHMMDB.split("/")[:-2])}/HHSuiteDB/'
        fname = f"{wdir}/mycluster"
        shell(f"mkdir {wdir}")

        '''Make HHsuite PPHMMs from protein alignments'''
        for Cluster_i in alive_it(range(len(Cluster_MetaDataDict))):
            AlnClusterFile = f"{ClustersDir}/Cluster_{Cluster_i}.fasta"
            HHsuite_PPHMMFile = f"{HHsuite_PPHMMDir}/PPHMM_{Cluster_i}.hmm"
            shell(f"hhmake -i {AlnClusterFile} -o {HHsuite_PPHMMFile} -seq {len(Cluster_MetaDataDict[Cluster_i]['Cluster'])+1} -name Cluster_{Cluster_i} -id 100 -M 50 -v 0",
                                 "PPHMMDB Construction: merge alignments, hhmake")

        '''Make PPHMM DB'''
        AlignmentMerging_i_round = self.rebuild_hhsuite_db(
            HHsuite_PPHMMDB, HHsuite_PPHMMDir, AlignmentMerging_i_round, ClustersDir, fname, first=True)

        '''Merge protein alignments'''
        while True:
            if self.HMMER_PPHMMDb_ForEachRoundOfPPHMMMerging == True:
                progress_msg(
                    f"\t\t - Building HMMER PPHMM DB...")
                _ = self.Make_HMMER_PPHMM_DB(HMMER_PPHMMDir=HMMER_PPHMMDir,
                                             HMMER_PPHMMDb=f"{HMMER_PPHMMDbDir}/HMMER_PPHMMDb_{AlignmentMerging_i_round}",
                                             ClustersDir=ClustersDir,
                                             Cluster_MetaDataDict=Cluster_MetaDataDict,
                                             VariableShelveFile=f"{VariableShelveDir}/PPHMMDBConstruction.p")

                shell(f"find {HMMER_PPHMMDir} -type f -name '*.hmm' -delete",
                        "PPHMMDB Construction: merge alignments, hhsearch")

            '''Inter-PPHMM similarity scoring'''
            progress_msg(f"\t - Round {AlignmentMerging_i_round + 1} - Determine PPHMM-PPHMM similarity scores (ALL-VERSUS-ALL hhsearch)")
            hhsearchDir = f"{HHsuiteDir}/hhsearch_{''.join(random.choice(string.ascii_uppercase + string.digits)for _ in range(10))}"
            os.makedirs(hhsearchDir)
            hhsearchOutFile = f"{hhsearchDir}/hhsearch.txt"

            '''Build hhm DBs'''
            progress_msg("\t\t - Building HMM Databases...")
            SeenPair, SeenPair_i, PPHMMSimScoreCondensedMat = {}, 0, []
            N_PPHMMs = LineCount(f"{fname}_hhm.ffindex")
            # for PPHMM_i in alive_it(range(0, N_PPHMMs)):
            for PPHMM_i in range(0, N_PPHMMs):
                HHsuite_PPHMMFile = f'{HHsuite_PPHMMDir}/PPHMM_{PPHMM_i}.hmm'
                hit_match = re.compile(
                    r"[A-Z0-9]+\|[A-Z0-9]+.[0-9]{1}\s[a-zA-Z0-9_ ]{0,10}")
                # RM < TODO PARAMETERISE MAXMEM
                shell(f"hhsearch -maxmem 9.0 -i {HHsuite_PPHMMFile} -d {fname} -o {hhsearchOutFile} -e {self.HHsuite_evalue_Cutoff} -E {self.HHsuite_evalue_Cutoff} -cov {self.HHsuite_SubjectCoverage_Cutoff} -z 1 -b 1 -id 100 -global -v 0 -cpu {self.HHsuite_N_CPUs}",
                                     "PPHMMDB Construction: merge alignments, hhsearch")
                with open(hhsearchOutFile, 'r') as hhsearchOut_txt:
                    Content = hhsearchOut_txt.readlines()
                    QueryLength = int(Content[1].split()[1])
                    for line in Content[9:]:
                        if line == "\n":
                            break
                        else:
                            '''Filter variable length Hit name'''
                            line_fil = re.sub(hit_match, "", line)
                            line_fil = re.sub(r" ~", "-", line_fil)
                            line_fil = line_fil.replace("(", " ").replace(")", " ").split()
                            try:
                                PPHMM_j = int(line_fil[0])
                                evalue = float(line_fil[3])
                                pvalue = float(line_fil[4])
                                PPHMMSimScore = float(line_fil[5])
                                Col = float(line_fil[7])
                                # SubjectLength = int(line_fil[-1])
                            except ValueError as ex:
                                warning_msg(f"Error parsing HHSuite output. Attempting to correct...")
                                try:
                                    line_fil = line.split()
                                    PPHMM_j = int(line_fil[0])
                                    evalue = float(line_fil[4])
                                    pvalue = float(line_fil[5])
                                    PPHMMSimScore = float(line_fil[6])
                                    Col = float(line_fil[8])
                                    print(f"...corrected successfully!")
                                except:
                                    raise_gravity_error(f"Failed to extract PPHMM data from HHsearch output (pphmmdb construction, aln merging function). This happens when HHsearch outputs malformed text files, check the output for the entry listed in this exception: {ex}")
                            qcovs = Col/QueryLength*100
                            #scovs = Col/SubjectLength*100
                            if ( evalue <= self.HHsuite_evalue_Cutoff and # RM < TODO Why would we need these when in HHsuite search terms?
                                 pvalue <= self.HHsuite_pvalue_Cutoff and # RM < TODO ... so check output
                                 qcovs >= self.HHsuite_QueryCoverage_Cutoff):
                                Pair = ", ".join(
                                    sorted(map(str, [PPHMM_i, PPHMM_j])))
                                if Pair in SeenPair:
                                    '''If the pair has already been seen...'''
                                    if PPHMMSimScore > PPHMMSimScoreCondensedMat[SeenPair[Pair]][2]:
                                        '''...and if the new PPHMMSimScore is higher...'''
                                        PPHMMSimScoreCondensedMat[SeenPair[Pair]
                                                                  ][2] = PPHMMSimScore
                                else:
                                    SeenPair[Pair] = SeenPair_i
                                    PPHMMSimScoreCondensedMat.append(
                                        [PPHMM_i, PPHMM_j - 1 , PPHMMSimScore]) # -1 added for zero indedxing
                                    SeenPair_i = SeenPair_i+1

                shell(f"rm {hhsearchOutFile}", "PPHMMDB Construction: merge alignments, remove hhsearch file")
            clean_stdout()

            '''Structure similarity score matrix for saving'''
            PPHMMSimScoreCondensedMat = np.array(PPHMMSimScoreCondensedMat)
            if PPHMMSimScoreCondensedMat.shape[0] == 0:
                raise ValueError(f"Failed on PPHMM Alignment Merging step: no alignments could be merged. Check your thresholds or disable this feature.")

            PPHMMSimScoreMat = np.zeros((N_PPHMMs, N_PPHMMs))
            PPHMMSimScoreMat[list(map(int, PPHMMSimScoreCondensedMat[:, 0])), list(map(
                int, PPHMMSimScoreCondensedMat[:, 1]))] = list(map(float, PPHMMSimScoreCondensedMat[:, 2]))
            PPHMMSimScoreMat[list(map(int, PPHMMSimScoreCondensedMat[:, 1])), list(map(
                int, PPHMMSimScoreCondensedMat[:, 0]))] = list(map(float, PPHMMSimScoreCondensedMat[:, 2]))
            PPHMMSimScoreCondensedMat = np.array(
                [PPHMMSimScorePair for PPHMMSimScorePair in PPHMMSimScoreCondensedMat if PPHMMSimScorePair[0] < PPHMMSimScorePair[1]])
            PPHMMSimScoreCondensedMatFile = f"{hhsearchDir}/PPHMMSimScoreCondensedMat.txt"
            np.savetxt(fname=PPHMMSimScoreCondensedMatFile,
                       X=PPHMMSimScoreCondensedMat,
                       fmt='%s',
                       delimiter="\t",
                       header="PPHMM_i\tPPHMM_j\tPPHMMSimScore")

            '''Cluster PPHMMs using Muscle'''
            progress_msg(
                "\t\t - Cluster PPHMMs based on hhsearch scores, using the MCL algorithm")
            PPHMMClustersFile = f"{hhsearchDir}/PPHMMClusters.txt"
            shell(f"mcl {PPHMMSimScoreCondensedMatFile} --abc -o {PPHMMClustersFile} -I {self.PPHMMClustering_MCLInflation}",
                                 "PPHMMDB Construction: merge alignments, mcl call")

            SeenProtIDList = []
            with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
                for Cluster in PPHMMClusters_txt.readlines():
                    SeenProtIDList.extend(Cluster.split("\n")[0].split("\t"))

            with open(PPHMMClustersFile, 'a') as PPHMMClusters_txt:
                PPHMMClusters_txt.write("\n".join(list(
                    set(map(str, list(map(float, list(range(0, N_PPHMMs))))))-set(SeenProtIDList))))

            progress_msg("\t\t - Check if there are alignments to be merged")
            with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
                N_PPHMMs_AfterMerging = len(PPHMMClusters_txt.readlines())

            if N_PPHMMs_AfterMerging == N_PPHMMs or AlignmentMerging_i_round == self.N_AlignmentMerging:
                progress_msg(
                    "\t\t\t - No alignments to be merged. Stop alignment merging process")
                shell(f"rm -rf {hhsearchDir}", "PPHMMDB Construction: merge alignments, rm hhsearch dir")
                break
            else:
                progress_msg(
                    f"\t\t\t - Merge {N_PPHMMs} alignments to make {N_PPHMMs_AfterMerging} alignments")

            '''Merge alignments, remake PPHMMs'''
            progress_msg("\t\t - Merge protein alignments and remake HHsuite PPHMMs")
            SelfSimScoreList = PPHMMSimScoreMat.diagonal()
            PPHMMDissimScoreMat = 1 - np.nan_to_num(np.transpose(PPHMMSimScoreMat**2 / SelfSimScoreList) / SelfSimScoreList)
            PPHMMDissimScoreMat[PPHMMDissimScoreMat < 0] = 0
            AfterMergingPPHMM_IndexList, AfterMergingPPHMM_i = [], 1

            with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
                for PPHMMCluster in alive_it(PPHMMClusters_txt.readlines()):
                    PPHMMCluster = list(
                        map(int, list(map(float, PPHMMCluster.split("\n")[0].split("\t")))))
                    AfterMergingPPHMM_IndexList.append(min(PPHMMCluster))
                    if len(PPHMMCluster) == 1:
                        pass

                    elif len(PPHMMCluster) >= 2:
                        PPHMMDissimScoreMat_Subset = PPHMMDissimScoreMat[PPHMMCluster][:, PPHMMCluster]
                        PPHMMTreeNewick = DistMat2Tree(DistMat=PPHMMDissimScoreMat_Subset,
                                                    LeafList=PPHMMCluster,
                                                    Dendrogram_LinkageMethod="average")
                        PPHMMTreeNewick = Tree(PPHMMTreeNewick)
                        _ = PPHMMTreeNewick.ladderize()
                        PPHMMTreeNewick = PPHMMTreeNewick.write(format=9)

                        while True:
                            '''Iterate over tree looking for pairs'''
                            m = re.search(r"\((\d+),(\d+)\)", PPHMMTreeNewick)
                            mafft_temp_fname = f"{self.ShelveDir}/temp.fasta"
                            if not m:
                                shell(f"mafft --thread {self.HHsuite_N_CPUs} --localpair --maxiterate 1000  {ClusterFile_i} > {mafft_temp_fname}",
                                        calling_fn="PPHMMDB Construction: merge alignments, mafft refine")
                                shell(f"rm {ClusterFile_i} && mv {mafft_temp_fname} {ClusterFile_i}",
                                        "PPHMMDB Construction: merge alignments, rm mafft intermediates")
                                break

                            PPHMM_i, PPHMM_j = sorted(
                                [int(m.group(1)), int(m.group(2))])
                            PPHMMTreeNewick = re.sub(
                                r"\((\d+),(\d+)\)", str(PPHMM_i), PPHMMTreeNewick, count=1)

                            ClusterFile_i = f"{ClustersDir}/Cluster_{PPHMM_i}.fasta"
                            ClusterFile_j = f"{ClustersDir}/Cluster_{PPHMM_j}.fasta"
                            HHsuite_PPHMMFile_j = f"{HHsuite_PPHMMDir}/PPHMM_{PPHMM_j}.hmm"

                            # RM < TODO BP HERE AND REPLACe WITH MAFFT
                            shell(f"muscle -profile -in1 {ClusterFile_i} -in2 {ClusterFile_j} -out {mafft_temp_fname} -gapopen {self.MUSCLE_GapOpenCost} -gapextend {self.MUSCLE_GapExtendCost}",
                                                 "PPHMMDB Construction: merge alignments, muscle profile")
                            shell(f"rm {ClusterFile_i} {ClusterFile_j} {HHsuite_PPHMMFile_j} && mv {mafft_temp_fname} {ClusterFile_i}",
                                    "PPHMMDB Construction: merge alignments, rm mafft intermediates")

                            Cluster_MetaDataDict[PPHMM_i]["Cluster"] = Cluster_MetaDataDict[PPHMM_i]["Cluster"] + \
                                Cluster_MetaDataDict[PPHMM_j]["Cluster"]
                            Cluster_MetaDataDict[PPHMM_i]["DescList"] = Cluster_MetaDataDict[PPHMM_i]["DescList"] + \
                                Cluster_MetaDataDict[PPHMM_j]["DescList"]
                            Cluster_MetaDataDict[PPHMM_i]["TaxoLists"] = Cluster_MetaDataDict[PPHMM_i]["TaxoLists"] + \
                                Cluster_MetaDataDict[PPHMM_j]["TaxoLists"]
                            del Cluster_MetaDataDict[PPHMM_j]

                        HHsuite_PPHMMFile_i = f"{HHsuite_PPHMMDir}/PPHMM_{PPHMM_i}.hmm"
                        Cluster_MetaDataDict[PPHMM_i]["AlignmentLength"] = AlignIO.read(
                            ClusterFile_i, "fasta").get_alignment_length()
                        shell(f"hhmake -i {ClusterFile_i} -o {HHsuite_PPHMMFile_i} -v 0 -seq {len(Cluster_MetaDataDict[PPHMM_i]['Cluster'])+1} -name Cluster_{PPHMM_i} -id 100 -M 50",
                                             "PPHMMDB Construction: merge alignments, hhmake")
                    AfterMergingPPHMM_i += 1

            '''Rename protein alignments and their associated PPHMMs'''
            AfterMergingPPHMM_IndexList, AfterMergingPPHMM_i = sorted(
                AfterMergingPPHMM_IndexList), 0
            for PPHMM_i in AfterMergingPPHMM_IndexList:
                Cluster_MetaDataDict[AfterMergingPPHMM_i] = Cluster_MetaDataDict.pop(
                    PPHMM_i)

                ClusterFile_i = f"{ClustersDir}/Cluster_{PPHMM_i}.fasta"
                ClusterFile_j = f"{ClustersDir}/Cluster_{AfterMergingPPHMM_i}.fasta"
                shell(f"mv {ClusterFile_i} {ClusterFile_j}",
                            "PPHMMDB Construction: merge alignments, move cluster files")

                HHsuite_PPHMMFile_i = f"{HHsuite_PPHMMDir}/PPHMM_{PPHMM_i}.hmm"
                HHsuite_PPHMMFile_j = f"{HHsuite_PPHMMDir}/PPHMM_{AfterMergingPPHMM_i}.hmm"
                shell(f"mv {HHsuite_PPHMMFile_i} {HHsuite_PPHMMFile_j}",
                            "PPHMMDB Construction: merge alignments, move hhsuite files")

                with open(HHsuite_PPHMMFile_j, "r+") as HHsuite_PPHMM_txt:
                    contents = HHsuite_PPHMM_txt.readlines()
                    contents[1] = "NAME  Cluster_%s\n" % AfterMergingPPHMM_i
                    contents = "".join(contents)
                    HHsuite_PPHMM_txt.seek(0)
                    HHsuite_PPHMM_txt.write(contents)
                    HHsuite_PPHMM_txt.truncate()
                AfterMergingPPHMM_i += 1

            '''Rebuild the HHsuite PPHMM database'''
            AlignmentMerging_i_round = self.rebuild_hhsuite_db(
                HHsuite_PPHMMDB, HHsuite_PPHMMDir, AlignmentMerging_i_round, ClustersDir, fname)
            shell(f"rm -rf {hhsearchDir}",
                        "PPHMMDB Contruction: rebuild hhsite db, remove hhsuite intermediates")

        progress_msg("\tAlignment merging is done.")
        '''Delete dir'''
        shell(f"rm -rf {HHsuiteDir}",
              "PPHMMDB Construction: merge alignments, remove hhsuite dir")

    def rebuild_hhsuite_db(self, HHsuite_PPHMMDB, HHsuite_PPHMMDir, AlignmentMerging_i_round, ClustersDir, fname, first=False):
        '''Rebuild the HHsuite PPHMM database'''
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

        if not first:
            AlignmentMerging_i_round += 1

        return AlignmentMerging_i_round


    def Make_HMMER_PPHMM_DB(self, HMMER_PPHMMDir, HMMER_PPHMMDb, ClustersDir, Cluster_MetaDataDict, VariableShelveFile):
        '''10/10: Build HMMER DB of PPHMMs, create and save summary file'''
        progress_msg("- Make HMMER PPHMMDB and its summary file")
        ClusterSizeList, ClusterSizeByTaxoGroupingList, ClusterSizeByProtList, ClusterTaxoList, ClusterProtSeqIDList, ClusterDescList, N_PPHMMs, pphhmmdb_construction_out = [
            ], [], [], [], [], [], len(Cluster_MetaDataDict), {}

        for PPHMM_i in alive_it(range(N_PPHMMs)):
            '''Iterate over each cluster, build annotations'''
            Cluster = Cluster_MetaDataDict[PPHMM_i]["Cluster"]
            DescList = Cluster_MetaDataDict[PPHMM_i]["DescList"]
            TaxoLists = Cluster_MetaDataDict[PPHMM_i]["TaxoLists"]
            ClusterSizeList.append(len(Cluster))
            ClusterSizeByTaxoGroupingList.append(", ".join([f"{TaxoGrouping}: {N_ProtSeqsInTheTaxoGroup}" for TaxoGrouping, N_ProtSeqsInTheTaxoGroup in sorted(list(Counter(list(zip(*TaxoLists))[-1]).items()),
                                                                                                                                                               key=operator.itemgetter(1),
                                                                                                                                                               reverse=True
                                                                                                                                                               )])
                                                 )
            ClusterSizeByProtList.append(", ".join([f"{Desc}: {N_ProtSeqsWithDesc}" for Desc, N_ProtSeqsWithDesc in sorted(list(Counter(DescList).items()),
                                                                                                                           key=operator.itemgetter(1),
                                                                                                                           reverse=True
                                                                                                                           )])
                                         )
            ClusterTaxo_UniqueBaltimoreGroup, ClusterTaxo_UniqueOrder, ClusterTaxo_UniqueFamily, \
                ClusterTaxo_UniqueSubFam, ClusterTaxo_UniqueGenus, _, _ = [
                    list(set(TaxoList)) for TaxoList in zip(*TaxoLists)]
            ClusterTaxo = []

            for ClusterTaxo_UniqueTaxoLabel in [ClusterTaxo_UniqueBaltimoreGroup, ClusterTaxo_UniqueOrder, ClusterTaxo_UniqueFamily, ClusterTaxo_UniqueSubFam, ClusterTaxo_UniqueGenus]:
                if len(ClusterTaxo_UniqueTaxoLabel) == 1:
                    if ClusterTaxo_UniqueTaxoLabel[0] not in ["", "Unassigned", "unassigned"]:
                        ClusterTaxo = ClusterTaxo + ClusterTaxo_UniqueTaxoLabel
                else:
                    ClusterTaxo = ClusterTaxo + \
                        ["/".join(ClusterTaxo_UniqueTaxoLabel)]
                    break

            ClusterTaxoList.append("; ".join(ClusterTaxo))
            ClusterProtSeqIDList.append(", ".join(Cluster))

            ClusterDescCount = sorted(
                list(Counter(DescList).items()), key=operator.itemgetter(1), reverse=True)
            for ClusterDesc in ClusterDescCount:
                ClusterDesc = ClusterDesc[0] # RM < TODO check in clustedesc.lower() after changing name of hypothetical clusters
                if (("Hypothetical" not in ClusterDesc) and ("hypothetical" not in ClusterDesc) and ("~" not in ClusterDesc)):
                    break

            ClusterDescList.append(f"{ClusterDesc}|{ClusterTaxoList[-1]}")

            '''Make a PPHMM using HMMER hmmbuild'''
            AlnClusterFile = f"{ClustersDir}/Cluster_{PPHMM_i}.fasta"
            HMMER_PPHMMFile = f"{HMMER_PPHMMDir}/PPHMM_{PPHMM_i}.hmm"
            shell(f"hmmbuild {HMMER_PPHMMFile} {AlnClusterFile}", "PPHMMDB Construction: make HMMER PPHMMDB, hmmbuild")

            '''Modify the DESC line in the HMM file to ClusterDesc|ClusterTaxo'''
            with open(HMMER_PPHMMFile, "r+") as HMMER_PPHMM_txt:
                contents = HMMER_PPHMM_txt.readlines()
                contents.insert(2, f"DESC  {ClusterDescList[-1]}\n")
                contents = "".join(contents)
                HMMER_PPHMM_txt.seek(0)
                HMMER_PPHMM_txt.write(contents)
                HMMER_PPHMM_txt.truncate()
        clean_stdout()

        '''Make a HMMER HMM DB with hmmpress'''
        shell(f"find {HMMER_PPHMMDir} -name '*.hmm' -exec cat {{}} \; > {HMMER_PPHMMDb}", "PPHMMDB Construction: make HMMER PPHMMDB, find db")
        shell(f"hmmpress {HMMER_PPHMMDb}", "PPHMMDB Construction: make HMMER PPHMMDB, HMMPress")

        '''Make a PPHMMDBSummary file'''
        ClusterIDList = [
            f"Cluster_{Cluster_i}" for Cluster_i in range(N_PPHMMs)]
        np.savetxt(fname=HMMER_PPHMMDb+"_Summary.txt",
                   X=np.column_stack((ClusterIDList,
                                       ClusterDescList,
                                       ClusterSizeList,
                                       ClusterProtSeqIDList,
                                       ClusterSizeByTaxoGroupingList,
                                       ClusterSizeByProtList)),
                   fmt='%s',
                   delimiter="\t",
                   header="# Cluster ID\tCluster desc\tSequence number\tProtein ID\tNumber of sequences by class\tNumber of sequences by protein")

        '''Prepare output and save as pickle file'''
        pphhmmdb_construction_out["ClusterIDList"] = np.array(ClusterIDList)
        pphhmmdb_construction_out["ClusterDescList"] = np.array(
            ClusterDescList)
        pphhmmdb_construction_out["ClusterSizeList"] = np.array(
            ClusterSizeList)
        pphhmmdb_construction_out["ClusterProtSeqIDList"] = np.array(
            ClusterProtSeqIDList)
        pphhmmdb_construction_out["ClusterSizeByTaxoGroupingList"] = np.array(
            ClusterSizeByTaxoGroupingList)
        pphhmmdb_construction_out["ClusterSizeByProtList"] = np.array(
            ClusterSizeByProtList)
        pickle.dump(pphhmmdb_construction_out, open(VariableShelveFile, "wb"))

    def main(self):
        '''Entrypoint to PPHMMDB construction functions'''
        section_header(
            "#Build a database of virus protein profile hidden Markov models (PPHMMs) #")

        '''1/10: Build db dirs'''
        BLASTQueryFile, BLASTSubjectFile, BLASTOutputFile, BLASTBitScoreFile, \
            BLASTProtClusterFile, ClustersDir, HMMER_PPHMMDir, HMMER_PPHMMDbDir, \
            HMMER_PPHMMDb, VariableShelveDir, HHsuiteDir, HHsuite_PPHMMDir, HHsuite_PPHMMDB = mkdir_pphmmdbc(self.ShelveDir)

        '''2/10: Retrieve Variables'''
        genomes = retrieve_genome_vars(
            VariableShelveDir, self.IncludeIncompleteGenomes)

        '''3/10: Get GenBank files if not exists'''
        GenBankDict = self.get_genbank(genomes)

        '''4/10: Sequence extraction'''
        ProtList, ProtIDList = self.sequence_extraction(genomes, GenBankDict)

        '''5/10: Make BLAST DB'''
        self.make_blastp_db(BLASTSubjectFile, ProtList)

        '''6/10: Do BLASTp analysis, save output'''
        self.blastp_analysis(ProtList, BLASTQueryFile,
                             BLASTSubjectFile, BLASTOutputFile, BLASTBitScoreFile)

        '''7/10: Cluster using Muscle'''
        self.mcl_clustering(ProtIDList, BLASTBitScoreFile,
                            BLASTProtClusterFile)

        '''8/10: Make Alignments'''
        Cluster_MetaDataDict = self.make_alignments(
            ProtList, ProtIDList, BLASTProtClusterFile, ClustersDir)

        '''9/10: Make PPHMMs, DB and Merge Alignments'''
        if self.N_AlignmentMerging != 0:
            self.pphmm_and_merge_alignments(Cluster_MetaDataDict, ClustersDir, HHsuite_PPHMMDir,
                                            HMMER_PPHMMDir, HMMER_PPHMMDbDir, HHsuite_PPHMMDB, HHsuiteDir, VariableShelveDir)

        '''10/10: Make HMMER DB, save to persistent storage'''
        self.Make_HMMER_PPHMM_DB(HMMER_PPHMMDir, HMMER_PPHMMDb, ClustersDir,
                                 Cluster_MetaDataDict, VariableShelveDir+"/PPHMMDBConstruction.p")
