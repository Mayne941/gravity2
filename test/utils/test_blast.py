import pytest
from unittest import mock
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from app.utils.blast import blastp_analysis
from app.utils.blast import eval_blast_query


@mock.patch("app.utils.blast.shell")
@mock.patch("app.utils.blast.SeqIO.write")
@mock.patch("app.utils.blast.pd.read_csv")
def test_blastp_analysis(mock_read_csv, mock_seqio_write, mock_shell):
    mock_prot_list = [SeqRecord(Seq("MKTIIALSYIFCLVFADYKDDDDK"), id="prot1"),
                      SeqRecord(Seq("MKTIIALSYIFCLVFADYKDDDDK"), id="prot2")]
    mock_blast_df = pd.DataFrame({
        "qseqid": ["prot1", "prot1"],
        "sseqid": ["prot2", "prot1"],
        "pident": [99.0, 99.0],
        "qcovs": [100, 100],
        "qlen": [24, 24],
        "slen": [24, 24],
        "evalue": [0.0, 0.0],
        "bitscore": [50.0, 50.0]
    })
    mock_read_csv.return_value = mock_blast_df
    mock_shell.return_value = ""

    fnames = {
        "MashSubjectFile": "test/data/blast_subject.fasta",
        "MashQueryFile": "test/data/blast_query.fasta",
        "MashOutputFile": "test/data/MashOutputFile.tab"
    }
    payload = {"N_CPUs": 1}

    result = blastp_analysis(mock_prot_list, fnames, payload, is_test=True)

    mock_seqio_write.assert_any_call(mock_prot_list, mock.ANY, "fasta")
    assert mock_shell.called
    assert isinstance(result, np.ndarray)
    assert result.shape == (1, 3)
    assert result[0][0] == "prot1"
    assert result[0][1] == "prot2"
    assert result[0][2] == "50.0"


def test_eval_blast_query():
    i = {"qseqid": "seq1", "sseqid": "seq2", "bitscore": 50}
    SeenPair = {}
    SeenPair_i = 0
    BitScoreMat = []

    BitScoreMat, SeenPair, SeenPair_i = eval_blast_query(i, SeenPair, SeenPair_i, BitScoreMat)

    assert SeenPair == {"seq1, seq2": 0}
    assert SeenPair_i == 1
    assert BitScoreMat == [["seq1", "seq2", 50]]

    i = {"qseqid": "seq1", "sseqid": "seq2", "bitscore": 60}
    BitScoreMat, SeenPair, SeenPair_i = eval_blast_query(i, SeenPair, SeenPair_i, BitScoreMat)

    assert SeenPair == {"seq1, seq2": 0}
    assert SeenPair_i == 1
    assert BitScoreMat == [["seq1", "seq2", 60]]

    i = {"qseqid": "seq1", "sseqid": "seq2", "bitscore": 40}
    BitScoreMat, SeenPair, SeenPair_i = eval_blast_query(i, SeenPair, SeenPair_i, BitScoreMat)

    assert SeenPair == {"seq1, seq2": 0}
    assert SeenPair_i == 1
    assert BitScoreMat == [["seq1", "seq2", 60]]

    i = {"qseqid": "seq3", "sseqid": "seq4", "bitscore": 70}
    BitScoreMat, SeenPair, SeenPair_i = eval_blast_query(i, SeenPair, SeenPair_i, BitScoreMat)

    assert SeenPair == {"seq1, seq2": 0, "seq3, seq4": 1}
    assert SeenPair_i == 2
    assert BitScoreMat == [["seq1", "seq2", 60], ["seq3", "seq4", 70]]
