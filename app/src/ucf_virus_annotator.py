from copy import copy
import os
import string
import random
import subprocess
import pickle
from scipy.sparse import coo_matrix

from app.utils.pphmm_signature_table_constructor import PPHMMSignatureTable_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_pickle
from app.utils.generate_fnames import generate_file_names
from app.utils.stdout_utils import progress_msg
from app.utils.shell_cmds import shell


class UcfVirusAnnotator:
    def __init__(self,
                 payload,
                 ExpDir,
                 ) -> None:
        '''fpaths'''
        self.payload = payload
        self.fnames = generate_file_names(payload, ExpDir, Pl2=True)
        self.genomes = retrieve_genome_vars(self.fnames['ReadGenomeDescTablePickle'])

    def annotate(self) -> None:
        '''3/3: Annotate unclassified viruses via PPHMM and GOM, generate loc/sig tables for unclassified, save to pickle'''
        progress_msg(f"Annotating unclassified viruses using the PPHMM and GOM databases of the reference viruses from Pipeline I")

        '''Retrieve variables related to the reference viruses'''
        pl1_ref_annotations = retrieve_pickle(self.fnames["Pl1RefAnnotatorPickle"])
        RefVirusGroup = self.fnames["ExpDir_Pl1"].split("/")[-1]
        if "GOMDB_coo" in pl1_ref_annotations.keys(): # RM < TODO Check if still needed
            pl1_ref_annotations["GOMDB"] = {GOMID: GOM_coo.toarray(
            ) for GOMID, GOM_coo in pl1_ref_annotations["GOMDB_coo"].items()}

        '''Generate PPHMMSignatureTable, PPHMMLocationTable, and GOMSignatureTable for unclassified viruses using the reference PPHMM and GOM database'''
        '''Make HMMER_hmmscanDir'''
        HMMER_hmmscanDir_RefVirus = f"{self.fnames['HMMERDir_RefVirus']}/hmmscan_"+''.join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        os.makedirs(HMMER_hmmscanDir_RefVirus)

        '''Generate PPHMMSignatureTable and PPHMMLocationTable'''
        PPHMMSignatureTable, PPHMMLocationTable = PPHMMSignatureTable_Constructor(
            SeqIDLists=self.genomes["SeqIDLists"],
            GenBankFile=self.payload['GenomeSeqFile_UcfVirus'],
            TranslTableList=self.genomes["TranslTableList"],
            HMMER_PPHMMDB=self.fnames['HMMER_PPHMMDB_RefVirus'],
            HMMER_hmmscanDir=HMMER_hmmscanDir_RefVirus,
            HMMER_N_CPUs=self.payload['N_CPUs'],
            HMMER_C_EValue_Cutoff=self.payload['HMMER_C_EValue_Cutoff'],
            HMMER_HitScore_Cutoff=self.payload['HMMER_HitScore_Cutoff']
        )

        '''Delete HMMER_hmmscanDir'''
        shell(f"rm -rf {HMMER_hmmscanDir_RefVirus}")

        GOMSignatureTable = GOMSignatureTable_Constructor(PPHMMLocationTable=PPHMMLocationTable,
                                                            GOMDB=pl1_ref_annotations["GOMDB"],
                                                            GOMIDList=pl1_ref_annotations["GOMIDList"])

        '''Construct out dict, dump to pickle'''
        # TODO Remove when checked that virus classifier still works with new use of "RefVirusGroup", i.e. no iteration
        # PPHMMSignatureTable_Dict, PPHMMLocationTable_Dict, GOMSignatureTable_Dict = {}, {}, {}
        # PPHMMSignatureTable_Dict[RefVirusGroup] = PPHMMSignatureTable
        # PPHMMLocationTable_Dict[RefVirusGroup] = PPHMMLocationTable
        # GOMSignatureTable_Dict[RefVirusGroup] = GOMSignatureTable
        # all_ucf_genomes["PPHMMSignatureTable_Dict_coo"] = {RefVirusGroup: coo_matrix(
        #     SigTable) for RefVirusGroup, SigTable in PPHMMSignatureTable}
        # all_ucf_genomes["PPHMMLocationTable_Dict_coo"] = {RefVirusGroup: coo_matrix(
        #     LocTable) for RefVirusGroup, LocTable in PPHMMLocationTable.items()}
        # all_ucf_genomes["GOMSignatureTable_Dict"] = GOMSignatureTable

        all_ucf_genomes = {}
        all_ucf_genomes["PPHMMSignatureTable_Dict_coo"] = {RefVirusGroup: coo_matrix(
            PPHMMSignatureTable)}
        all_ucf_genomes["PPHMMLocationTable_Dict_coo"] = {RefVirusGroup: coo_matrix(
            PPHMMLocationTable)}
        all_ucf_genomes["GOMSignatureTable_Dict"] = GOMSignatureTable

        pickle.dump(all_ucf_genomes, open(self.fnames['UcfAnnotatorPickle'], "wb"))

    def main(self) -> None:
        '''Generate PPHMM signature table, PPHMM location table, and GOM signature table for unclassified viruses, using PPHMM databases of reference viruses'''
        section_header(
            "Generate PPHMM signature table, PPHMM location table, and GOM signature table\nfor unclassified viruses, using PPHMM databases of reference viruses")
        '''Annotate unclassified viruses using PPHMM & GOM'''
        self.annotate()
