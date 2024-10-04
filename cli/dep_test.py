from app.utils.shell_cmds import shell
from app.utils.stdout_utils import progress_msg
from app.utils.error_handlers import error_handle_mafft, error_handler_hmmbuild, error_handler_mash_dist, error_handler_blast, error_handler_mcl

def test_imports():
    try:
        import numpy
        import ete3
        import Bio
        import scipy
        import matplotlib
        import sklearn
        import fastapi
        import uvicorn
        import pandas
        import bs4
        import requests
        import openpyxl
        import dendropy
        import alive_progress
        import pyarrow
        progress_msg("... Python library load successful")
    except:
        raise BaseException(f"Failed to import Python libraries. Ensure you've installed"
                            f"libraries with `pip install -r requirements.txt` and"
                            f"have activated your conda environment.")

class DepTest:
    def __init__(self) -> None:
        pass

    def test_mash(self):
        out = shell("mash", ret_output=True)
        error_handler_mash_dist(out, "Mash (dependency test)")
        progress_msg("... Mash call successful")

    def test_blast(self):
        out = shell("blastp", ret_output=True)
        error_handler_blast(out, "Blast (dependency test)")
        progress_msg("... BLAST call successful")

    def test_mcl(self):
        out = shell("mcl", ret_output=True)
        error_handler_mcl(out, "Mcl (dependency test)")
        progress_msg("... Mcl call successful")

    def test_mafft(self):
        out = shell("mafft --help", ret_output=True)
        error_handle_mafft(out, "Mafft (dependency test)")
        progress_msg("... Mafft call successful")

    def test_hmmer(self):
        out = shell("hmmpress", ret_output=True)
        error_handler_hmmbuild(out, "Hmmer (dependency test)")
        progress_msg("... Hmmer call successful")

    def test_hhsuite(self):
        out = shell("cstranslate", ret_output=True)
        error_handler_hmmbuild(out, "hhsuite (dependency test)")
        progress_msg("... hhsuite call successful")

    def main(self):
        progress_msg(f"Starting dependency tests...")
        test_imports()
        self.test_mash()
        self.test_mafft()
        self.test_blast()
        self.test_mcl()
        self.test_hmmer()
        self.test_hhsuite()
        progress_msg(f"Successfully passed all tests, you're set to start using GRAViTy-V2!")

if __name__ == "__main__":
    clf = DepTest()
    clf.main()
