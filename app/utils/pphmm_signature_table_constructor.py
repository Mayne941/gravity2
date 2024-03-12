from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random, string, os
from alive_progress import alive_it

from .line_count import LineCount
from app.utils.shell_cmds import shell
from app.utils.stdout_utils import warning_msg, progress_msg
from app.utils.orf_identifier import find_orfs
from app.utils.error_handlers import raise_gravity_error

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
    if Pl2:
        PPHMMDB_Summary = f"{HMMER_PPHMMDB}_Summary.txt"
        HMMER_hmmscanDir = f"{fnames['HMMERDir_UcfVirus']}/hmmscan_{''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))}"
    else:
        PPHMMDB_Summary = f"{fnames['HMMER_PPHMMDb']}_Summary.txt"
        HMMER_hmmscanDir = f"{fnames['HMMERDir']}/hmmscan_{''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))}"
    os.makedirs(HMMER_hmmscanDir)
    PPHMMQueryFile = f'{HMMER_hmmscanDir}/QProtSeqs.fasta'
    PPHMMScanOutFile = f'{HMMER_hmmscanDir}/PPHMMScanOut.txt'

    '''Load GenBank record'''
    GenBankDict = SeqIO.index(GenomeSeqFile, "fasta" if os.path.splitext(GenomeSeqFile)[1] in [".fas", ".fst", ".fasta"] else "gb")
    Records_dict = {}
    for i in GenBankDict.items():
        ### RM < TODO: TEST BACK TRANSCRIBE IF RNA SEQS USED
        if "u" in str(i[1].seq).lower():
            i[1].seq = i[1].seq.back_transcribe()
        Records_dict[i[0].split(".")[0]] = i[1]
    Records_dict = {k.split(".")[0]: v for k, v in Records_dict.items()}

    N_PPHMMs = LineCount(PPHMMDB_Summary)-1
    PPHMMSignatureTable = np.empty((0, N_PPHMMs))
    PPHMMLocMiddleBestHitTable = np.empty((0, N_PPHMMs))

    for SeqIDList, TranslTable, BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping in alive_it(zip(genomes["SeqIDLists"], genomes["TranslTableList"], genomes["BaltimoreList"], genomes["OrderList"], genomes["FamilyList"], genomes["SubFamList"], genomes["GenusList"], genomes["VirusNameList"], genomes["TaxoGroupingList"])):
        GenBankSeqList, GenBankIDList = [], []
        for SeqID in SeqIDList:
            GenBankRecord = Records_dict[SeqID]
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
        # for GenBankSeq, GenBankID in zip(GenBankSeqList, GenBankIDList):
        prot, prot_id, _ = find_orfs(GenBankIDList, GenBankSeqList, TranslTable, payload['ProteinLength_Cutoff'], call_locs=True,
                                        taxonomy_annots=[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping])
        ProtList += prot
        ProtIDList += prot_id ## TODO de nest this

        if len(ProtList) < 1:
            raise_gravity_error(f"GRAViTy couldn't detect any reading frames in your input sequence(s) ({ProtIDList}). Check that your sequences are labelled properly.")

        with open(PPHMMQueryFile, "w") as f:
            '''Write each translated ORF to the PPHMM Query File'''
            for i in range(len(ProtIDList)):
                f.write(f">{ProtIDList[i]}\n{str(ProtList[i].seq)}\n")

        shell(f"hmmscan --cpu {payload['N_CPUs']} -E {payload['HMMER_C_EValue_Cutoff']} --noali --nobias --domtblout {PPHMMScanOutFile} {HMMER_PPHMMDB} {PPHMMQueryFile}")

        PPHMMIDList, PPHMMScoreList, FeatureFrameBestHitList, FeatureLocFromBestHitList, \
            FeatureLocToBestHitList, FeatureDescList = [], [], [], [], [], []
        with open(PPHMMScanOutFile, "r") as PPHMMScanOut_txt:
            for Line in PPHMMScanOut_txt:
                if Line[0] == "#":
                    '''If header'''
                    continue
                Line = Line.split()
                '''Concatenate the cluster description back'''
                Line[22] = " ".join(Line[22:])
                Line = Line[:23]
                C_EValue = float(Line[11])
                HitScore = float(Line[7])
                OriAASeqlen = float(len(GenBankSeqList))/3

                if HitScore <= payload['HMMER_HitScore_Cutoff']:
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
        FeatureLocMiddleBestHitList = np.zeros(N_PPHMMs)
        FeatureLocMiddleBestHitList[PPHMMIDList] = np.mean(np.array([FeatureLocFromBestHitList, FeatureLocToBestHitList]), axis=0)*(
            np.array(FeatureFrameBestHitList)/abs(np.array(FeatureFrameBestHitList)))
        PPHMMLocMiddleBestHitTable = np.vstack(
            (PPHMMLocMiddleBestHitTable, FeatureLocMiddleBestHitList))
        FeatureValueList = np.zeros(N_PPHMMs)
        FeatureValueList[PPHMMIDList] = PPHMMScoreList
        PPHMMSignatureTable = np.vstack(
            (PPHMMSignatureTable, FeatureValueList))

    '''Delete temp HMMER dir'''
    shell(f"rm -rf {HMMER_hmmscanDir}")

    return (PPHMMSignatureTable, PPHMMLocMiddleBestHitTable)
