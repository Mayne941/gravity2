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
from app.utils.orf_identifier import get_orf_trasl_table, no_orf_match, find_orfs
from app.utils.stdout_utils import clean_stdout, progress_msg, warning_msg
from app.utils.retrieve_pickle import retrieve_genome_vars
from app.utils.shell_cmds import shell
from app.utils.mkdirs import mkdir_pphmmdbc
from app.utils.error_handlers import raise_gravity_error
from app.utils.taxo_label_constructor import TaxoLabel_Constructor
from app.utils.generate_fnames import generate_file_names
from app.utils.blast import blastp_analysis

class PPHMMDBConstruction:
    def __init__(self,
                 payload,
                 GenomeSeqFile,
                 ExpDir,
                 ) -> None:
        self.payload = payload
        self.fnames = generate_file_names(payload,ExpDir)
        self.GenomeSeqFile = GenomeSeqFile
        self.orf_tranl_table = get_orf_trasl_table()
        self.orf_no_match = no_orf_match()
        self.genomes = retrieve_genome_vars(self.fnames["ReadGenomeDescTablePickle"])
        mkdir_pphmmdbc(self.fnames)

    def get_genbank(self):
        '''3/10: Check if GenBank files exist; dl if needed. Index & transform.'''
        if not os.path.isfile(self.GenomeSeqFile):
            progress_msg("Download GenBank file\n GenomeSeqFile doesn't exist. GRAViTy is downloading the GenBank file(s).")
            DownloadGenBankFile(self.GenomeSeqFile, # RM < TODO We should give users option to provide their own gb file if not on genbank
                                self.genomes["SeqIDLists"], self.payload["genbank_email"])

        progress_msg("Reading GenBank file")
        GenBankDict = SeqIO.index(self.GenomeSeqFile, "genbank")
        CleanGenbankDict = {} # TODO remove redundancy with sig table constructor
        for i in GenBankDict.items():
            ### RM < TODO TEST: BACK TRANSCRIBE IF RNA SEQS USED
            if "u" in str(i[1].seq).lower():
                i[1].seq = i[1].seq.back_transcribe()
            CleanGenbankDict[i[0].split(".")[0]] = i[1]

        return {k.split(".")[0]: v for k, v in CleanGenbankDict.items()}

    def sequence_extraction(self, GenBankDict):
        '''4/10: Extract protein sequences from VMR using GenBank data; manually annotate ORFs if needed'''
        progress_msg(
            f"- Extract/predict protein sequences from virus genomes, excluding proteins with lengthes <{self.payload['ProteinLength_Cutoff']} aa")
        ProtList, ProtIDList, Virus_i = [
        ], [], 1
        raw_seqs = {}

        for SeqIDList, TranslTable, BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping in alive_it(zip(self.genomes["SeqIDLists"], self.genomes["TranslTableList"], self.genomes["BaltimoreList"], self.genomes["OrderList"], self.genomes["FamilyList"], self.genomes["SubFamList"], self.genomes["GenusList"], self.genomes["VirusNameList"], self.genomes["TaxoGroupingList"])):

            for SeqID in SeqIDList:
                '''Sometimes an Acc ID doesn't have a matching record (usually when multiple seqs for 1 virus)... - skip if true'''
                try:
                    GenBankRecord = GenBankDict[SeqID]
                except KeyError as ex:
                    raise_gravity_error(f"{ex}\n"
                                        f"ERROR: I couldn't extract a sequence ID from the input Genbank file.\n"
                                        f"This usually happens when you've tried to re-run the experiment with an old .gb file.\n"
                                        f"Try deleting your input .gb file ({self.GenomeSeqFile}), then starting again.")

                GenBankID = GenBankRecord.name
                GenBankFeatures = GenBankRecord.features

                '''Extract protein sequences'''
                ContainProtAnnotation = False
                for Feature in GenBankFeatures:
                    if(Feature.type == 'CDS' and "protein_id" in Feature.qualifiers and "translation" in Feature.qualifiers):
                        ContainProtAnnotation = True
                        try: # TODO Triple exception makes me sad
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
                        if len(ProtSeq) >= self.payload['ProteinLength_Cutoff']:
                            ProtRecord = SeqRecord(Seq(ProtSeq),
                                                #    id   = f"{Family}_{Genus}_{VirusName}_{GenBankID}|{ProtID}", # Print taxonomic information to mashup subjects file
                                                #    name =f"{Family}_{Genus}_{VirusName}_{GenBankID}|{ProtID}",
                                                    id=f"{GenBankID}|{ProtID}",
                                                    name=f"{GenBankID}|{ProtID}",
                                                    description=ProtName,
                                                    annotations={'taxonomy': [BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping]})
                            ProtList.append(ProtRecord)
                            ProtIDList.append(f"{GenBankID}|{ProtID}")

                '''If the genome isn't annotated with any ORFs, find some'''
                if not ContainProtAnnotation:
                    prots, prot_ids, raw_nas = find_orfs(GenBankID, GenBankRecord.seq, TranslTable, self.payload['ProteinLength_Cutoff'],
                                                    taxonomy_annots=[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping])
                    ProtList += prots
                    ProtIDList += prot_ids
                    raw_seqs.update(raw_nas)
            Virus_i += 1
        clean_stdout()

        if len(ProtList) == 0:
            raise_gravity_error(f"No protein sequences could be extracted from input sequences. Did you input untranslated RNA sequences?")

        with open(self.fnames['MashSubjectFile'], "w") as MashSubject_txt:
            SeqIO.write(ProtList, MashSubject_txt, "fasta") # RM < TODO Swap to pure py?

        '''Save ref seqs for comparative analysis'''
        ref_seqs = [[i, str(GenBankDict[i].seq)] for i in GenBankDict.keys()]
        TaxoLabelList = TaxoLabel_Constructor(SeqIDLists=self.genomes["SeqIDLists"],
                                              FamilyList=self.genomes["FamilyList"],
                                              GenusList=self.genomes["GenusList"],
                                              VirusNameList=self.genomes["VirusNameList"]
                                              )
        for ref_idx in range(len(ref_seqs)):
            hits = [s for s in TaxoLabelList if ref_seqs[ref_idx][0] in s]
            if len(hits) > 0 and len(hits) < 2:
                ref_seqs[ref_idx][0] = hits[0]
        with open(self.fnames['RefSeqFile'], "w") as f: [f.write(f">{i[0]}\n{i[1]}\n") for i in ref_seqs]

        return ProtList, np.array(ProtIDList)

    def mash_analysis(self, ProtList):
        '''6/10: Perform ALL-VERSUS-ALL Mash analysis'''
        progress_msg("Creating Mash sketches")
        SeenPair, SeenPair_i, MashMatrix, N_ProtSeqs = {}, 0, [], len(ProtList)
        shell(f"mash-Linux64-v2.3/mash sketch -p {self.payload['N_CPUs']} -a -i {self.fnames['MashSubjectFile']}") # TODO PARAMETERISE MASH CALL; add parallelism (-p)
        for ProtSeq_i in alive_it(range(N_ProtSeqs)):
            '''Mash query fasta file'''
            MashQuery = ProtList[ProtSeq_i]
            with open(self.fnames['MashQueryFile'], "w") as MashQuery_txt:
                MashQuery_txt.write(f">{MashQuery.name} {MashQuery.description}\n{str(MashQuery.seq)}")

            mash_fname = f'{"/".join(self.fnames["MashOutputFile"].split("/")[:-1])}/mashup_scores.tab'
            shell(f"mash-Linux64-v2.3/mash sketch -p {self.payload['N_CPUs']} -a -i {self.fnames['MashQueryFile']}")
            shell(f"mash-Linux64-v2.3/mash dist -p {self.payload['N_CPUs']} -i {self.fnames['MashSubjectFile']}.msh {self.fnames['MashQueryFile']}.msh > {mash_fname}")
            # TODO Pandas read might be slowing this down a lot? Can we pipe STDOUT to Python?
            mash_df = pd.read_csv(mash_fname, sep="\t", header=None, names=["orf", "query", "dist", "p", "hashes"])

            mash_df["query"] = ProtList[ProtSeq_i].id
            mash_df = mash_df[mash_df["orf"] != mash_df["query"]]
            mash_df["p"] = mash_df["p"].astype(float)
            mash_df["mash_sim"] = np.abs(1 - mash_df["dist"])

            mash_iter = mash_df.to_dict(orient="records")
            for i in mash_iter:
                pair = ", ".join(sorted([ProtList[ProtSeq_i].id, i["orf"]]))
                if  ((i["p"] < self.payload['Mash_p_val_cutoff']) # TODO Put reporting thresholds into mash call
                    and (i["mash_sim"] > self.payload['Mash_sim_score_cutoff'])):
                    '''Query must: not match subject, have identity > thresh, have query coverage > thresh and query coverage normalised to subject length > thresh'''
                    if pair in SeenPair:
                        '''If the pair has already been seen...'''
                        if i["mash_sim"] > MashMatrix[SeenPair[pair]][2]:
                            '''...and if the new mash_sim is higher'''
                            MashMatrix[SeenPair[pair]][2] = i["mash_sim"]
                    else:
                        SeenPair[pair] = SeenPair_i
                        MashMatrix.append(
                            [pair.split(", ")[0], pair.split(", ")[1], i["mash_sim"]])
                        SeenPair_i += 1

        MashMatrix = np.array(MashMatrix)

        if MashMatrix.shape[0] == 0:
            # TODO PARAMETERISE THIS!!
            try:
                blastp_analysis(ProtList, self.fnames, self.payload)
            except:
                raise_gravity_error(f"No Mash results were extracted from protein sequences. Check Mash is installed and that your parameters aren't too restrictive."
                                    f"Turning mash similarity below 0.4 can make false positives frequent: possibly your sequence is very unlike what's in your VMR?"
                                    f"If so, try enabling BLASTp instead, which will replace the Mash call.")
        np.savetxt(fname=self.fnames['MashSimFile'],
                   X=MashMatrix,
                   fmt='%s',
                   delimiter="\t",
                   header="SeqID_I\tSeqID_II\tBit score"
                   )

    def mcl_clustering(self, ProtIDList):
        '''7/10: Use Mcl to do clustering on Mash bit scores'''
        progress_msg("- Doing protein sequence clustering based on Mash bit scores, using the MCL algorithm")
        shell(f"mcl {self.fnames['MashSimFile']} --abc -o {self.fnames['MashProtClusterFile']} -I {self.payload['ProtClustering_MCLInflation']}",
                    "PPHMMDB Contruction: mcl clustering, call mcl")

        '''For each cluster found, pull out seen proteins by ID'''
        SeenProtIDList = []
        with open(self.fnames['MashProtClusterFile'], 'r') as MashProtCluster_txt:
            for Cluster in MashProtCluster_txt.readlines():
                SeenProtIDList.extend(Cluster.split(
                    "\r\n")[0].split("\n")[0].split("\t"))

        '''Append seen protein IDs to clusters'''
        with open(self.fnames['MashProtClusterFile'], 'a') as MashProtCluster_txt:
            MashProtCluster_txt.write(
                "\n".join(list(set(ProtIDList)-set(SeenProtIDList))))

    def make_alignments(self, ProtList, ProtIDList):
        '''8/10: Do protein alignments with Mafft, make cluster alignment annotations'''
        progress_msg("- Make protein alignments")
        _, Cluster_i, Cluster_MetaDataDict = LineCount(
            self.fnames['MashProtClusterFile'])+1, 0, {}
        with open(self.fnames['MashProtClusterFile'], 'r') as MashProtCluster_txt:
            for Cluster in alive_it(MashProtCluster_txt.readlines()):
                HitList, TaxoLists, DescList, Cluster = [], [], [], Cluster.split("\n")[
                    0].split("\t")
                for ProtID in Cluster:
                    HitList.append(
                        ProtList[np.where(ProtIDList == ProtID)[0][0]])
                    TaxoLists.append(HitList[-1].annotations['taxonomy'])
                    DescList.append(HitList[-1].description.replace(", ", " ").replace(",", " ").replace(": ", "_").replace(
                        ":", "_").replace("; ", " ").replace(";", " ").replace(" (", "/").replace("(", "/").replace(")", ""))

                '''Cluster file'''
                AlnClusterFile = f"{self.fnames['ClustersDir']}/Cluster_{Cluster_i}.fasta"
                with open(AlnClusterFile, "w") as UnAlnClusterTXT: ################# TODO this likely very slow.
                    p = SeqIO.write(HitList, UnAlnClusterTXT, "fasta")
                # with open(AlnClusterFile, "w") as UnAlnClusterTXT: # In progress......
                #     UnAlnClusterTXT.write(f">{AlnClusterFile.name} {AlnClusterFile.description}\n{str(AlnClusterFile.seq)}")

                temp_aln_fname = f"{self.fnames['ClustersDir']}/temp.fasta"
                '''Align cluster using Mafft'''
                shell(f"mafft --thread {self.payload['N_CPUs']} --localpair --maxiterate 1000  {AlnClusterFile} > {temp_aln_fname}",
                        "PPHMMDB Construction: make alignments, main Mafft call")
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

    def pphmm_and_merge_alignments(self, Cluster_MetaDataDict):
        '''9/10 (OPTIONAL): Make HHSuite PPHMMs and DB, merge & make protein alignmens'''
        if self.payload['N_AlignmentMerging'] > 0:
            progress_msg(
                f"- Merge protein alignments, {self.payload['N_AlignmentMerging']} rounds of merging")
        elif self.payload['N_AlignmentMerging'] < 0:
            progress_msg("- Merge protein alignments until exhausted")

        AlignmentMerging_i_round = 0
        wdir = f'{"/".join(self.fnames["HHsuite_PPHMMDB"].split("/")[:-2])}/HHSuiteDB/'
        fname = f"{wdir}/mycluster"
        shell(f"mkdir {wdir}")

        '''Make HHsuite PPHMMs from protein alignments'''
        for Cluster_i in alive_it(range(len(Cluster_MetaDataDict))):
            AlnClusterFile = f"{self.fnames['ClustersDir']}/Cluster_{Cluster_i}.fasta"
            HHsuite_PPHMMFile = f"{self.fnames['HHsuite_PPHMMDir']}/PPHMM_{Cluster_i}.hmm"
            shell(f"hhmake -i {AlnClusterFile} -o {HHsuite_PPHMMFile} -seq {len(Cluster_MetaDataDict[Cluster_i]['Cluster'])+1} -name Cluster_{Cluster_i} -id 100 -M 50 -v 0",
                                 "PPHMMDB Construction: merge alignments, hhmake")

        '''Make PPHMM DB'''
        AlignmentMerging_i_round = self.rebuild_hhsuite_db(AlignmentMerging_i_round, fname, first=True)

        '''Merge protein alignments'''
        while True:
            if self.payload['HMMER_PPHMMDb_ForEachRoundOfPPHMMMerging'] == True:
                progress_msg(
                    f"\t\t - Building HMMER PPHMM DB...")
                _ = self.Make_HMMER_PPHMM_DB(HMMER_PPHMMDb=f"{self.fnames['HMMER_PPHMMDbDir']}/HMMER_PPHMMDb_{AlignmentMerging_i_round}",
                                             Cluster_MetaDataDict=Cluster_MetaDataDict)

                shell(f"find {self.fnames['HMMER_PPHMMDir']} -type f -name '*.hmm' -delete",
                        "PPHMMDB Construction: merge alignments, hhsearch")

            '''Inter-PPHMM similarity scoring'''
            progress_msg(f"\t - Round {AlignmentMerging_i_round + 1} - Determine PPHMM-PPHMM similarity scores (ALL-VERSUS-ALL hhsearch)")
            hhsearchDir = f"{self.fnames['HHsuiteDir']}/hhsearch_{''.join(random.choice(string.ascii_uppercase + string.digits)for _ in range(10))}"
            os.makedirs(hhsearchDir)
            hhsearchOutFile = f"{hhsearchDir}/hhsearch.txt"

            '''Build hhm DBs'''
            progress_msg("\t\t - Building HMM Databases...")
            SeenPair, SeenPair_i, PPHMMSimScoreCondensedMat = {}, 0, []
            N_PPHMMs = LineCount(f"{fname}_hhm.ffindex")
            for PPHMM_i in alive_it(range(0, N_PPHMMs)):
                HHsuite_PPHMMFile = f'{self.fnames["HHsuite_PPHMMDir"]}/PPHMM_{PPHMM_i}.hmm'
                hit_match = re.compile(
                    r"[A-Z0-9]+\|[A-Z0-9]+.[0-9]{1}\s[a-zA-Z0-9_ ]{0,10}")
                # RM < TODO PARAMETERISE MAXMEM
                shell(f"hhsearch -maxmem 9.0 -i {HHsuite_PPHMMFile} -d {fname} -o {hhsearchOutFile} -e {self.payload['HHsuite_evalue_Cutoff']} -E {self.payload['HHsuite_evalue_Cutoff']} -cov {self.payload['HHsuite_SubjectCoverage_Cutoff']} -z 1 -b 1 -id 100 -global -v 0 -cpu {self.payload['N_CPUs']}",
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
                            if ( evalue <= self.payload['HHsuite_evalue_Cutoff'] and # RM < TODO Why would we need these when in HHsuite search terms?
                                 pvalue <= self.payload['HHsuite_pvalue_Cutoff'] and # RM < TODO ... so check output
                                 qcovs >= self.payload['HHsuite_QueryCoverage_Cutoff']):
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

            '''Cluster PPHMMs with Mcl'''
            progress_msg(
                "\t\t - Cluster PPHMMs based on hhsearch scores, using the MCL algorithm")
            PPHMMClustersFile = f"{hhsearchDir}/PPHMMClusters.txt"
            shell(f"mcl {PPHMMSimScoreCondensedMatFile} --abc -o {PPHMMClustersFile} -I {self.payload['PPHMMClustering_MCLInflation']}",
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

            if N_PPHMMs_AfterMerging == N_PPHMMs or AlignmentMerging_i_round == self.payload['N_AlignmentMerging']:
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
                                shell(f"mafft --thread {self.payload['N_CPUs']} --localpair --maxiterate 1000  {ClusterFile_i} > {mafft_temp_fname}",
                                        calling_fn="PPHMMDB Construction: merge alignments, mafft refine")
                                shell(f"rm {ClusterFile_i} && mv {mafft_temp_fname} {ClusterFile_i}",
                                        "PPHMMDB Construction: merge alignments, rm mafft intermediates")
                                break

                            PPHMM_i, PPHMM_j = sorted(
                                [int(m.group(1)), int(m.group(2))])
                            PPHMMTreeNewick = re.sub(
                                r"\((\d+),(\d+)\)", str(PPHMM_i), PPHMMTreeNewick, count=1)

                            ClusterFile_i = f"{self.fnames['ClustersDir']}/Cluster_{PPHMM_i}.fasta"
                            ClusterFile_j = f"{self.fnames['ClustersDir']}/Cluster_{PPHMM_j}.fasta"
                            HHsuite_PPHMMFile_j = f"{self.fnames['HHsuite_PPHMMDir']}/PPHMM_{PPHMM_j}.hmm"

                            # RM < TODO REPLACE WITH MAFFT
                            shell(f"muscle -profile -in1 {ClusterFile_i} -in2 {ClusterFile_j} -out {mafft_temp_fname} -gapopen -3.0 -gapextend -0.0",
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

                        HHsuite_PPHMMFile_i = f"{self.fnames['HHsuite_PPHMMDir']}/PPHMM_{PPHMM_i}.hmm"
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

                ClusterFile_i = f"{self.fnames['ClustersDir']}/Cluster_{PPHMM_i}.fasta"
                ClusterFile_j = f"{self.fnames['ClustersDir']}/Cluster_{AfterMergingPPHMM_i}.fasta"
                shell(f"mv {ClusterFile_i} {ClusterFile_j}",
                            "PPHMMDB Construction: merge alignments, move cluster files")

                HHsuite_PPHMMFile_i = f"{self.fnames['HHsuite_PPHMMDir']}/PPHMM_{PPHMM_i}.hmm"
                HHsuite_PPHMMFile_j = f"{self.fnames['HHsuite_PPHMMDir']}/PPHMM_{AfterMergingPPHMM_i}.hmm"
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
            AlignmentMerging_i_round = self.rebuild_hhsuite_db(AlignmentMerging_i_round, fname)
            shell(f"rm -rf {hhsearchDir}",
                        "PPHMMDB Contruction: rebuild hhsite db, remove hhsuite intermediates")

        progress_msg("\tAlignment merging is done.")
        '''Delete dir'''
        shell(f"rm -rf {self.fnames['HHsuiteDir']}",
              "PPHMMDB Construction: merge alignments, remove hhsuite dir")

    def rebuild_hhsuite_db(self, AlignmentMerging_i_round, fname, first=False):
        '''Rebuild the HHsuite PPHMM database'''
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

        if not first:
            AlignmentMerging_i_round += 1

        return AlignmentMerging_i_round


    def Make_HMMER_PPHMM_DB(self, PphmmDb, Cluster_MetaDataDict):
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
            AlnClusterFile = f"{self.fnames['ClustersDir']}/Cluster_{PPHMM_i}.fasta"
            HMMER_PPHMMFile = f"{self.fnames['HMMER_PPHMMDir']}/PPHMM_{PPHMM_i}.hmm"
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
        shell(f"find {self.fnames['HMMER_PPHMMDir']} -name '*.hmm' -exec cat {{}} \; > {PphmmDb}", "PPHMMDB Construction: make HMMER PPHMMDB, find db")
        shell(f"hmmpress -f {PphmmDb}", "PPHMMDB Construction: make HMMER PPHMMDB, HMMPress")

        '''Make a PPHMMDBSummary file'''
        ClusterIDList = [
            f"Cluster_{Cluster_i}" for Cluster_i in range(N_PPHMMs)]
        np.savetxt(fname=f'{PphmmDb}_Summary.txt',
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
        pickle.dump(pphhmmdb_construction_out, open(self.fnames["PphmmdbPickle"], "wb"))

    def main(self):
        '''Entrypoint to PPHMMDB construction functions'''
        section_header(
            "Build a database of virus protein profile hidden Markov models (PPHMMs)")

        '''3/10: Get GenBank files if not exists'''
        GenBankDict = self.get_genbank()

        '''4/10: Sequence extraction'''
        ProtList, ProtIDList = self.sequence_extraction(GenBankDict)

        '''6/10: Do Mash analysis, save output'''
        self.mash_analysis(ProtList)

        '''7/10: Cluster using Mcl'''
        self.mcl_clustering(ProtIDList)

        '''8/10: Make Alignments'''
        Cluster_MetaDataDict = self.make_alignments(
            ProtList, ProtIDList)

        '''9/10: Make PPHMMs, DB and Merge Alignments'''
        if self.payload['N_AlignmentMerging'] != 0:
            self.pphmm_and_merge_alignments(Cluster_MetaDataDict)

        '''10/10: Make HMMER DB, save to persistent storage'''
        self.Make_HMMER_PPHMM_DB(self.fnames['HMMER_PPHMMDb'], Cluster_MetaDataDict)
