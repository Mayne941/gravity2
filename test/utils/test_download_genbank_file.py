import pytest
import os
from unittest.mock import patch, mock_open, MagicMock
from app.utils.download_genbank_file import DownloadGenBankFile
from app.utils.error_handlers import raise_gravity_error

@patch("app.utils.download_genbank_file.Entrez.efetch")
@patch("builtins.open", new_callable=mock_open)
@patch("os.makedirs")
@patch("os.path.exists")
def test_DownloadGenBankFile(mock_exists, mock_makedirs, mock_open, mock_efetch):
    GenomeSeqFile = "/path/to/genome_seq_file.gb"
    SeqIDLists = [["ABC123", "DEF456"], ["GHI789"]]
    email = "test@example.com"

    mock_exists.return_value = False
    mock_handle = MagicMock()
    mock_handle.read.return_value = "GenBank data"
    mock_efetch.return_value = mock_handle

    DownloadGenBankFile(GenomeSeqFile, SeqIDLists, email)

    mock_exists.assert_called_once_with("/path/to")
    mock_makedirs.assert_called_once_with("/path/to")
    mock_efetch.assert_called_once_with(db="nucleotide", id="ABC123, DEF456, GHI789", rettype="gb", retmode="text")
    mock_open.assert_called_once_with(GenomeSeqFile, "w")
    mock_handle.read.assert_called_once()
    mock_handle.close.assert_called_once()
