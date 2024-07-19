from Bio import SeqIO
from alive_progress import alive_it
import numpy as np
import pandas as pd

from app.utils.stdout_utils import clean_stdout, progress_msg
from app.utils.shell_cmds import shell
from app.utils.error_handlers import raise_gravity_warning, raise_gravity_error

def blastp_analysis(ProtList, fnames, payload):
    '''6/10: Perform ALL-VERSUS-ALL BLASTp analysis'''
    raise_gravity_warning("No similar matches were found using Mash. Entering failover with a less sensitive search method (BLASTp)...")
    progress_msg("Performing ALL-VERSUS-ALL BLASTp analysis")
    with open(fnames['MashSubjectFile'], "w") as BLASTSubject_txt:
        SeqIO.write(ProtList, BLASTSubject_txt, "fasta")
    shell(f"makeblastdb -in {fnames['MashSubjectFile']} -dbtype prot", "PPHMMDB Construction: make BLASTp db")

    BitScoreMat, SeenPair, SeenPair_i, N_ProtSeqs = [
        ], {}, 0, len(ProtList)
    for ProtSeq_i in alive_it(range(N_ProtSeqs)):#############
    # for ProtSeq_i in range(N_ProtSeqs):
        '''BLAST query fasta file'''
        BLASTQuery = ProtList[ProtSeq_i]
        with open(fnames['MashQueryFile'], "w") as BLASTQuery_txt: _ = SeqIO.write(BLASTQuery, BLASTQuery_txt, "fasta")
        mash_fname = f'{"/".join(fnames["MashOutputFile"].split("/")[:-1])}/mashup_scores.tab'
        '''Perform BLASTp, load output to dataframe'''
        outfmat = '"6 qseqid sseqid pident qcovs qlen slen evalue bitscore"'
        shell(f'blastp -query {fnames["MashQueryFile"]} -db {fnames["MashSubjectFile"]} -out {mash_fname} -evalue 1E-6 -outfmt {outfmat} -num_alignments 1000000 -num_threads {payload["N_CPUs"]}',
                f"PPHMMDB Construction: BLASTp (protein ID = {ProtList[ProtSeq_i].id})")

        try:
            blast_df = pd.read_csv(mash_fname, sep="\t", names=["qseqid", "sseqid", "pident", "qcovs", "qlen", "slen", "evalue", "bitscore"])
        except Exception as e:
            raise raise_gravity_error(f"Could not open BLASTp output with exception: {e}")

        if blast_df.empty:
            '''If no hits in sequence, continue'''
            continue

        '''Load BLASTp results, collate bit scores for matches'''
        blast_iter = blast_df.to_dict(orient="records")
        for i in blast_iter:
            pair = ", ".join(sorted([i["qseqid"], i["sseqid"]]))

            # if i["slen"] / i["qlen"] < 0.6 and i["slen"] / i["qlen"] > 0.10: # TODO PUT IN SEPARATE LOGIC FOR V DIFF LEN ORFS
            if min(i["slen"],i["qlen"]) / max(i["slen"],i["qlen"]) < 0.6 and i["slen"] / i["qlen"] > 0.10: # TODO PUT IN SEPARATE LOGIC FOR V DIFF LEN ORFS

                if ((i["qseqid"] != i["sseqid"]) and
                    # ((i["qcovs"]*i["qlen"]/i["slen"]) >= 30)):
                    ((i["pident"]*i["qlen"]/i["slen"]) >= 30)): # WAS 60
                    # ((i["qcovs"]/i["pident"]*i["qlen"]/i["slen"]) >= 0.5)): # TODO THIS WORKED BUT WAS GREEDY

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
            else:
                if ((i["qseqid"] != i["sseqid"]) and
                    (i["pident"] >= 70) and  # was 50 / 70
                    (i["qcovs"] >= 90)       # was 75 / 90
                    and ((i["qcovs"]*i["qlen"]/i["slen"]) >= 95)): # was 75 / 95
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
    return np.array(BitScoreMat)
