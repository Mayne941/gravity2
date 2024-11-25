from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random, string, os
from alive_progress import alive_it
import time
from tqdm import tqdm
from multiprocessing import Pool

from app.utils.line_count import LineCount
from app.utils.shell_cmds import shell
from app.utils.stdout_utils import warning_msg, progress_msg
from app.utils.orf_identifier import find_orfs
from app.utils.error_handlers import raise_gravity_error, error_handler_hmmscan

def PPHMMSignatureTable_Constructor(
            genomes,
            payload,
            fnames,
            GenomeSeqFile,
            HMMER_PPHMMDB,
            Pl2=False,
        ):
    progress_msg("- Generating PPHMM signature table and PPHMM location table")
    '''All Hmmer dirs are temporary, so are generated dynamically then removed'''
    PPHMMDB_Summary = f"{fnames['HMMER_PPHMMDb']}_Summary.txt"
    HMMER_hmmscanDir = f"{fnames['HMMERDir']}/hmmscan_{''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))}"
    os.makedirs(HMMER_hmmscanDir)

    '''Load GenBank record'''
    GenBankDict = SeqIO.index(GenomeSeqFile, "fasta" if os.path.splitext(GenomeSeqFile)[1] in [".fas", ".fst", ".fasta"] else "gb")
    Records_dict = {}
    for i in GenBankDict.items():
        ### RM < TODO: EITHER HARMONISE FUNTION CALL WITH PPHMMDB OR LOAD GB FILE
        if "u" in str(i[1].seq).lower():
            i[1].seq = i[1].seq.back_transcribe()
        Records_dict[i[0].split(".")[0]] = i[1]
    Records_dict = {k.split(".")[0]: v for k, v in Records_dict.items()}

    N_PPHMMs = LineCount(PPHMMDB_Summary)-1
    PPHMMSignatureTable = np.empty((0, N_PPHMMs))
    PPHMMLocMiddleBestHitTable = np.empty((0, N_PPHMMs))
    NaiveLocationTable = np.empty((0, N_PPHMMs))

    clf = Pphmm_Sig_Gen(payload, Records_dict, HMMER_PPHMMDB, N_PPHMMs, HMMER_hmmscanDir)
    pool = Pool(payload["N_CPUs"])

    progress_msg(f"-  Spinning up {payload['N_CPUs']} workers to generate PPHMM signatures. This may take a while...")
    with pool as p, tqdm(total=genomes["SeqIDLists"].shape[0]) as pbar:
        res = [p.apply_async( # can specify specific inputs with args in a set!
            clf.generate_sigs_for_genome, args=(i,), callback=lambda _: pbar.update(1)) for i in genomes["SeqIDLists"]]
        results = [r.get() for r in res]

    for result in results:
        # Naive.. must be updated before PPhmmloc..
        NaiveLocationTable = np.vstack((PPHMMLocMiddleBestHitTable, result[2]))
        PPHMMLocMiddleBestHitTable = np.vstack((PPHMMLocMiddleBestHitTable, result[1]))
        PPHMMSignatureTable = np.vstack((PPHMMSignatureTable, result[3]))

    '''Delete temp HMMER dir'''
    shell(f"rm -rf {HMMER_hmmscanDir}")
    # RM < TODO SAVE TO PICKLE FOR HOT START
    return PPHMMSignatureTable, PPHMMLocMiddleBestHitTable, NaiveLocationTable

class Pphmm_Sig_Gen:
    def __init__(self, payload, Records_dict, HMMER_PPHMMDB, N_PPHMMs, HMMER_hmmscanDir) -> None:
        self.payload = payload
        self.Records_dict = Records_dict
        self.HMMER_PPHMMDB = HMMER_PPHMMDB
        self.HMMER_hmmscanDir = HMMER_hmmscanDir
        self.N_PPHMMs = N_PPHMMs

    def generate_sigs_for_genome(self, SeqIDList, TranslTable=1):
        PPHMMQueryFile = f'{self.HMMER_hmmscanDir}/QProtSeqs_{SeqIDList[0]}.fasta'
        PPHMMScanOutFile = f'{self.HMMER_hmmscanDir}/PPHMMScanOut_{SeqIDList[0]}.txt'
        GenBankSeqList, GenBankIDList = [], []
        for SeqID in SeqIDList:
            GenBankRecord = self.Records_dict[SeqID]
            GenBankSeqList.append(GenBankRecord.seq)
            GenBankIDList.append(GenBankRecord.id)

        '''Sort lists by sequence/segment lengths, then concat to single seq'''
        _, GenBankSeqList, GenBankIDList = list(zip(
            *sorted(zip(list(map(len, list(map(str, GenBankSeqList)))), GenBankSeqList, GenBankIDList), reverse=True))) # TODO WAS TRUE 17/02

        if len(GenBankSeqList) > 1:
            GenBankSeqList = sum(GenBankSeqList, Seq(""))
            GenBankIDList = "/".join([i.split(".")[0] for i in SeqIDList])
        else:
            GenBankSeqList = GenBankSeqList[0]
            GenBankIDList = GenBankIDList[0]

        '''Get each orf for a genome'''
        ProtList, ProtIDList = [], []
        prot, prot_id, _ = find_orfs(GenBankIDList, GenBankSeqList, TranslTable, self.payload['ProteinLength_Cutoff'], call_locs=True)
        ProtList += prot
        ProtIDList += prot_id

        if len(ProtList) < 1:
            raise_gravity_error(f"GRAViTy couldn't detect any reading frames in your input sequence(s) (Protein IDs {ProtIDList}; Accessions {SeqIDList}). Check that your sequences are labelled properly; some viroids can break this process if they have no detectable ORFs.")

        # BEFORE TODO
        with open(PPHMMQueryFile, "w") as f:
            '''Write each translated ORF to the PPHMM Query File'''
            for i in range(len(ProtIDList)):
                f.write(f">{ProtIDList[i]}\n{str(ProtList[i].seq)}\n")

        out = shell(f"hmmscan --cpu 1 -E {self.payload['HMMER_C_EValue_Cutoff']} --noali --nobias --domtblout {PPHMMScanOutFile} {self.HMMER_PPHMMDB} {PPHMMQueryFile}",
                    ret_output=True)
        error_handler_hmmscan(out, "hmmscan (PPHMM signature table constructor, PphmmSignatureTable_Constructor_DEPRECATED())")

        PPHMMIDList, PPHMMScoreList, FeatureFrameBestHitList, FeatureLocFromBestHitList, \
            FeatureLocToBestHitList, FeatureDescList = [], [], [], [], [], []
        with open(PPHMMScanOutFile, "r") as PPHMMScanOut_txt:
            for Line in PPHMMScanOut_txt:
                if Line[0] == "#":
                    '''If header'''
                    continue
                Line = Line.split()
                '''Concatenate the cluster description back'''
                try:
                    Line[22] = " ".join(Line[22:])
                except:
                    break # TODO TEST - sometimes the line is only half formed. Not sure why, possibly if no matches?
                Line = Line[:23]
                C_EValue = float(Line[11])
                HitScore = float(Line[7])
                OriAASeqlen = float(len(GenBankSeqList))/3

                if HitScore <= self.payload['HMMER_HitScore_Cutoff']:
                    '''Threshold at user-set values'''
                    continue

                '''Determine the frame and the location of the hit'''
                iden = int(Line[0].split('_')[-1])
                HitFrom = int(Line[3].split('|')[-1].replace("START",""))
                HitTo = HitFrom + int(Line[5])
                HitMid = float(HitFrom+HitTo)/2
                Frame = int(np.ceil(HitMid/OriAASeqlen)) if np.ceil(HitMid/OriAASeqlen) <= 3 else int(-(np.ceil(HitMid/OriAASeqlen)-3))
                LocFrom = int(HitFrom % OriAASeqlen)

                if LocFrom == 0:
                    '''if the hit occurs preciously from the end of the sequence'''
                    LocFrom = int(OriAASeqlen)
                LocTo = int(HitTo % OriAASeqlen)

                if LocTo == 0:
                    '''if the hit occurs preciously to the end of the sequence'''
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
                    Best_C_EValue = C_EValue
                    PPHMMIDList.append(iden)
                    PPHMMScoreList.append(HitScore)
                    FeatureDescList.append(Line[22].split('|')[0])
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
            NaiveLocationList = np.zeros(self.N_PPHMMs)
            NaiveLocationList[PPHMMIDList] = np.mean(np.array([FeatureLocFromBestHitList, FeatureLocToBestHitList]), axis=0)
            FeatureLocMiddleBestHitList = np.zeros(self.N_PPHMMs)
            FeatureLocMiddleBestHitList[PPHMMIDList] = np.mean(np.array([FeatureLocFromBestHitList, FeatureLocToBestHitList]), axis=0)*(
                np.array(FeatureFrameBestHitList)/abs(np.array(FeatureFrameBestHitList)))
            FeatureValueList = np.zeros(self.N_PPHMMs)
            FeatureValueList[PPHMMIDList] = PPHMMScoreList

        return (SeqIDList, FeatureLocMiddleBestHitList, NaiveLocationList, FeatureValueList)
        #            0              1                           2               3
        # NaiveLocationTable = np.vstack((PPHMMLocMiddleBestHitTable, result[2]))
        # PPHMMLocMiddleBestHitTable = np.vstack((PPHMMLocMiddleBestHitTable, result[1]))
        # PPHMMSignatureTable = np.vstack((PPHMMSignatureTable, result[3]))
