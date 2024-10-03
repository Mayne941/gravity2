import numpy as np
import re
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import warnings
from Bio import BiopythonWarning # Silence biopython warning when len(seq) % 3 != 0
warnings.simplefilter('ignore', BiopythonWarning)

from app.utils.error_handlers import raise_gravity_warning

def get_orf_trasl_table():
    return {1: "---M------**--*----M---------------M----------------------------",
            2: "----------**--------------------MMMM----------**---M------------",
            3: "----------**----------------------MM----------------------------",
            4: "--MM------**-------M------------MMMM---------------M------------",
            5: "---M------**--------------------MMMM---------------M------------",
            6: "--------------*--------------------M----------------------------",
            7: "--MM------**-------M------------MMMM---------------M------------",
            8: "---M------**--*----M---------------M----------------------------",
            9: "----------**-----------------------M---------------M------------",
            10: "----------**-----------------------M----------------------------",
            11: "---M------**--*----M------------MMMM---------------M------------",
            12: "----------**--*----M---------------M----------------------------",
            13: "---M------**----------------------MM---------------M------------",
            14: "-----------*-----------------------M----------------------------",
            15: "----------*---*--------------------M----------------------------",
            16: "----------*---*--------------------M----------------------------",
            17: "---M------**--*----M---------------M----------------------------",
            18: "---M------**--*----M---------------M----------------------------",
            19: "---M------**--*----M---------------M----------------------------",
            20: "---M------**--*----M---------------M----------------------------",
            21: "----------**-----------------------M---------------M------------",
            22: "------*---*---*--------------------M----------------------------",
            23: "--*-------**--*-----------------M--M---------------M------------",
            24: "---M------**-------M---------------M---------------M------------",
            25: "---M------**-----------------------M---------------M------------",
            26: "----------**--*----M---------------M----------------------------",
            27: "--------------*--------------------M----------------------------",
            28: "----------**--*--------------------M----------------------------",
            29: "--------------*--------------------M----------------------------",
            30: "--------------*--------------------M----------------------------",
            31: "----------**-----------------------M----------------------------",
        }

def no_orf_match():
    return "---M------**--*----M---------------M----------------------------"


def find_orfs(seq_id, GenBankSeq, TranslTable, protein_length_cutoff, call_locs = False, taxonomy_annots=[
            "None", "None", "None", "None", "None", "None", "None"
        ]):
    orf_tranl_table = get_orf_trasl_table()
    orf_no_match = no_orf_match()
    ProtList, ProtIDList, raw_nas = [], [], {}
    raw_nas[seq_id] = []
    '''Identifying ORFs'''
    try:
        Starts = orf_tranl_table[TranslTable]
    except KeyError:
        '''No ORF found, use standard code'''
        Starts = orf_no_match

    '''Generate start and stop codon lists'''
    CodonList = [
        Base1+Base2+Base3 for Base1 in "TCAG" for Base2 in "TCAG" for Base3 in "TCAG"]
    StartCodonList, StopCodonList = [], []
    for i, j in enumerate(Starts):
        if j == "M":
            StartCodonList.append(CodonList[i])
        if j == "*":
            StopCodonList.append(CodonList[i])

    SeqLength, ORF_i = len(GenBankSeq), 0
    for _, nuc in [(+1, GenBankSeq), (-1, GenBankSeq.reverse_complement())]:
        '''Split into multiple of 3, get in-frame nucleotide seq, split sequence into codons'''
        raw_nas[seq_id].append(str(GenBankSeq.translate()))
        for frame in range(3):
            length = 3 * ((SeqLength-frame) // 3)
            nuc_inframe = nuc[frame:(frame+length)]
            nuc_codonList = [str(nuc_inframe[i:i+3])
                                for i in range(0, length, 3)]

            '''Find stop codons'''
            StopCodon_indices = [i for i, codon in enumerate(
                nuc_codonList) if codon in StopCodonList]
            Coding_Start_IndexList = np.array(
                [-1]+StopCodon_indices)+1
            Coding_End_IndexList = np.array(
                StopCodon_indices+[len(nuc_codonList)])

            ProtSeqList, loc_list = [], []
            for i, j in zip(Coding_Start_IndexList, Coding_End_IndexList):
                for k, codon in enumerate(nuc_codonList[i:j]):
                    if codon in StartCodonList:
                        ProtSeqList.append(
                            Seq("".join(nuc_codonList[i:j][k:])).translate(table=TranslTable))
                        loc_list.append(i)
                        break


            for idx, ProtSeq in enumerate(ProtSeqList):
                '''Exclude protein sequences with <'ProteinLength_Cutoff' aa'''
                if not call_locs:
                    '''If on first call, count Xs and report to user to alert to mistranslated regions'''
                    x_count = ProtSeq.count("X")
                    if x_count > 10:
                        raise_gravity_warning(f"(PPHMMDB Construction: orf_identifier()). Putative ORF {seq_id}_{idx} has {x_count} mistranslated residues.")

                if len(ProtSeq.replace("X","")) >= protein_length_cutoff:
                    if ProtSeq.count("X") > 100:
                        raise_gravity_warning(f"TRIMMING {seq_id}|ORF{ORF_i} OF {ProtSeq.count('X')} MISTRANSLATED RESIDUES. Recommend user checks input sequence.")
                        ProtSeq = ProtSeq.replace("X","")
                    n = 1000000 # RM < TODO WIP, subsetting of large orfs
                    if len(ProtSeq) > n:
                        ProtSeqs = [ProtSeq[i:i + n] for i in range(0, len(ProtSeq), n)]
                        for orf_idx, ProtSeq in enumerate(ProtSeqs):
                            if not len(ProtSeq) < protein_length_cutoff:
                                ProtRecord = SeqRecord(ProtSeq,
                                                        id=f"{seq_id}|ORF{ORF_i}.{orf_idx}",
                                                        name=f"{seq_id}|ORF{ORF_i}.{orf_idx}",
                                                        description="~",
                                                        annotations={'taxonomy': taxonomy_annots})
                                ProtList.append(ProtRecord)
                                if call_locs:
                                    ProtIDList.append(
                                        f"{seq_id}|ORF{ORF_i}.{orf_idx}|START{loc_list[idx]}")
                                else:
                                    ProtIDList.append(
                                        f"{seq_id}|ORF{ORF_i}.{orf_idx}")

                    else:
                        ProtRecord = SeqRecord(ProtSeq,
                                                id=f"{seq_id}|ORF{ORF_i}.0",
                                                name=f"{seq_id}|ORF{ORF_i}.0",
                                                description="~",
                                                annotations={'taxonomy': taxonomy_annots})
                        ProtList.append(ProtRecord)
                        if call_locs:
                            ProtIDList.append(
                                f"{seq_id}|ORF{ORF_i}.0|START{loc_list[idx]}")
                        else:
                            ProtIDList.append(
                                f"{seq_id}|ORF{ORF_i}.0")
                    ORF_i += 1

    return ProtList, ProtIDList, raw_nas
