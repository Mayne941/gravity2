import pytest
from unittest.mock import patch, MagicMock
from app.utils.error_handlers import (
    raise_gravity_error,
    raise_gravity_warning,
    error_handler_virus_classifier,
    main_error_msg,
    decode_stdout,
    error_handler_hmmbuild,
    error_handle_mafft,
    error_handler_mash_sketch,
    error_handler_mash_dist,
    error_handler_hmmscan,
    error_handler_blast,
    error_handler_mcl,
    error_handler_hhsuite
)

def test_raise_gravity_error():
    with pytest.raises(SystemExit):
        raise_gravity_error("Test error message")

def test_raise_gravity_warning(capfd):
    raise_gravity_warning("Test warning message")
    out, err = capfd.readouterr()
    assert "WARNING" in out

def test_error_handler_virus_classifier():
    payload = {
        'UseUcfVirusPPHMMs': True,
        'GenomeSeqFiles_RefVirus': None,
        'GenomeSeqFile_UcfVirus': None
    }
    with pytest.raises(ValueError):
        error_handler_virus_classifier(payload)

def test_main_error_msg():
    assert main_error_msg("test_command") == "Critical error in command line call to test_command."

def test_decode_stdout():
    assert decode_stdout(b"test output") == "test output"

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handler_hmmbuild(mock_decode_stdout):
    mock_decode_stdout.return_value = "Usage: hmmbuild"
    error_handler_hmmbuild(b"test output", "test message")

    mock_decode_stdout.return_value = "input parse error"
    with pytest.raises(SystemExit):
        error_handler_hmmbuild(b"test output", "test message")

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handle_mafft(mock_decode_stdout):
    mock_decode_stdout.return_value = "mafft --maxiterate 1000"
    error_handle_mafft(b"test output", "dependency test")

    mock_decode_stdout.return_value = "killed"
    with pytest.raises(SystemExit):
        error_handle_mafft(b"test output", "test message")

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handler_mash_sketch(mock_decode_stdout):
    mock_decode_stdout.return_value = "Sketching Writing"
    error_handler_mash_sketch(b"test output", "test message")

    mock_decode_stdout.return_value = ""
    with pytest.raises(SystemExit):
        error_handler_mash_sketch(b"test output", "test message")

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handler_mash_dist(mock_decode_stdout):
    mock_decode_stdout.return_value = "mash <command> [options]"
    error_handler_mash_dist(b"test output", "test message")

    mock_decode_stdout.return_value = "exception"
    with pytest.raises(SystemExit):
        error_handler_mash_dist(b"test output", "test message")

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handler_hmmscan(mock_decode_stdout):
    mock_decode_stdout.return_value = "[ok]"
    error_handler_hmmscan(b"test output", "test message")

    mock_decode_stdout.return_value = "exception"
    with pytest.raises(SystemExit):
        error_handler_hmmscan(b"test output", "test message")

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handler_blast(mock_decode_stdout):
    mock_decode_stdout.return_value = "BLAST query/options error:"
    error_handler_blast(b"test output", "test message")

    mock_decode_stdout.return_value = "not blank but not installed error"
    with pytest.raises(SystemExit):
        error_handler_blast(b"test output", "test message")

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handler_mcl(mock_decode_stdout):
    mock_decode_stdout.return_value = "[mcl] new tab created"
    error_handler_mcl(b"test output", "test message")

    mock_decode_stdout.return_value = ""
    with pytest.raises(SystemExit):
        error_handler_mcl(b"test output", "test message")

@patch("app.utils.error_handlers.decode_stdout")
def test_error_handler_hhsuite(mock_decode_stdout):
    mock_decode_stdout.return_value = "Usage: cstranslate"
    error_handler_hhsuite(b"test output", "test message")

    mock_decode_stdout.return_value = "not blank and not usage error"
    with pytest.raises(SystemExit):
        error_handler_hhsuite(b"test output", "test message")

if __name__ == "__main__":
    pytest.main()
