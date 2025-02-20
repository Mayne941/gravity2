from Bio import SeqIO
from alive_progress import alive_it
import numpy as np
import pandas as pd

from app.utils.stdout_utils import progress_msg
from app.utils.shell_cmds import shell
from app.utils.error_handlers import raise_gravity_warning, raise_gravity_error, error_handler_blast

def eval_blast_query(i, SeenPair, SeenPair_i, BitScoreMat):
    pair = ", ".join(sorted([i["qseqid"], i["sseqid"]]))
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

    return BitScoreMat, SeenPair, SeenPair_i

def blastp_analysis(ProtList, fnames, payload, is_test=False):
    '''6/10: Perform ALL-VERSUS-ALL BLASTp analysis'''
    raise_gravity_warning("Performing all-vs-all BLASTp analysis. If you didn't select this as an option, it's because there were no Mash hits (you might need to refine your settings).")
    progress_msg("Performing ALL-VERSUS-ALL BLASTp analysis")
    with open(fnames['MashSubjectFile'], "w") as BLASTSubject_txt:
        SeqIO.write(ProtList, BLASTSubject_txt, "fasta")
    shell(f"makeblastdb -in {fnames['MashSubjectFile']} -dbtype prot", "PPHMMDB Construction: make BLASTp db")

    BitScoreMat, SeenPair, SeenPair_i, N_ProtSeqs = [
        ], {}, 0, len(ProtList)
    for ProtSeq_i in alive_it(range(N_ProtSeqs)):
        '''BLAST query fasta file'''
        BLASTQuery = ProtList[ProtSeq_i]
        with open(fnames['MashQueryFile'], "w") as BLASTQuery_txt: _ = SeqIO.write(BLASTQuery, BLASTQuery_txt, "fasta")
        mash_fname = f'{"/".join(fnames["MashOutputFile"].split("/")[:-1])}/mashup_scores.tab'
        '''Perform BLASTp, load output to dataframe'''
        outfmat = '"6 qseqid sseqid pident qcovs qlen slen evalue bitscore"'
        out = shell(f'blastp -query {fnames["MashQueryFile"]} -db {fnames["MashSubjectFile"]} -out {mash_fname} -evalue 1E-6 -outfmt {outfmat} -num_alignments 1000000 -num_threads {payload["N_CPUs"]}',
                ret_output=True)
        if not is_test:
            error_handler_blast(out, "BLASTp, PPHMMDB construction")

        try:
            blast_df = pd.read_csv(mash_fname, sep="\t", names=["qseqid", "sseqid", "pident", "qcovs", "qlen", "slen", "evalue", "bitscore"])
        except Exception as e:
            raise raise_gravity_error(f"Could not open BLASTp output with exception: {e}")

        if blast_df.empty:
            '''If no hits in sequence, continue'''
            continue

        '''Load BLASTp results, collate bit scores for matches'''
        blast_iter = blast_df.to_dict(orient="records")
        DISCREPANT_SIZE_RANGE = (0.10, 0.60)
        NORM_PIDENT_discrepant_sizes = 30 # WAS 30 - 60
        NORM_PIDENT_similar_sizes = 95    # was 75 - 95
        PIDENT = 70                       # was 50 - 70
        QCOV = 90                         # was 75 - 90
        for i in blast_iter:

            if i["slen"] / i["qlen"] > DISCREPANT_SIZE_RANGE[0] and (min(i["slen"],i["qlen"]) / max(i["slen"],i["qlen"])) <= DISCREPANT_SIZE_RANGE[1]:
                if ((i["qseqid"] != i["sseqid"]) and
                    ((i["pident"]*i["qlen"]/i["slen"]) >= NORM_PIDENT_discrepant_sizes)):

                    BitScoreMat, SeenPair, SeenPair_i = eval_blast_query(i, SeenPair, SeenPair_i, BitScoreMat)
            else:
                if ((i["qseqid"] != i["sseqid"]) and
                    (i["pident"] >= PIDENT) and
                    (i["qcovs"] >= QCOV) and
                    ((i["qcovs"]*i["qlen"]/i["slen"]) >= NORM_PIDENT_similar_sizes)):
                    '''Query must: not match subject, have identity > thresh, have query coverage > thresh and query coverage normalised to subject length > thresh'''

                    BitScoreMat, SeenPair, SeenPair_i = eval_blast_query(i, SeenPair, SeenPair_i, BitScoreMat)

    return np.array(BitScoreMat)
