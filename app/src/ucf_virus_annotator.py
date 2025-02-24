import os
import string
import random
import pickle
from scipy.sparse import coo_matrix

from app.utils.pphmm_signature_table_constructor import PPHMMSignatureTable_Constructor
from app.utils.gom_signature_table_constructor import GOMSignatureTable_Constructor
from app.utils.console_messages import section_header
from app.utils.retrieve_pickle import retrieve_genome_vars, retrieve_pickle
from app.utils.generate_fnames import generate_file_names
from app.utils.stdout_utils import progress_msg


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
        '''GOM Coords might not be in PL1 data depending on settings used'''
        if "GOMDB_coo" in pl1_ref_annotations.keys():
            pl1_ref_annotations["GOMDB"] = {GOMID: GOM_coo.toarray(
            ) for GOMID, GOM_coo in pl1_ref_annotations["GOMDB_coo"].items()}

        '''Generate PPHMMSignatureTable, PPHMMLocationTable, and GOMSignatureTable for unclassified viruses using the reference PPHMM and GOM database'''
        '''Generate PPHMMSignatureTable and PPHMMLocationTable'''
        PPHMMSignatureTable, PPHMMLocationTable, NaivePPHMMLocationTable = PPHMMSignatureTable_Constructor(
                                                self.genomes,
                                                self.payload,
                                                self.fnames,
                                                GenomeSeqFile=self.payload['GenomeSeqFile_UcfVirus'],
                                                HMMER_PPHMMDB=self.fnames['HMMER_PPHMMDB_RefVirus'],
                                                Pl2=True
                                            )

        GOMSignatureTable = GOMSignatureTable_Constructor(PPHMMLocationTable=PPHMMLocationTable,
                                                            GOMDB=pl1_ref_annotations["GOMDB"],
                                                            GOMIDList=pl1_ref_annotations["GOMIDList"],
                                                            payload=self.payload)

        '''Construct out dict, dump to pickle'''
        all_ucf_genomes = {}
        all_ucf_genomes["NaivePPHMMLocationTable"] = NaivePPHMMLocationTable
        all_ucf_genomes["PPHMMSignatureTable_coo"] = coo_matrix(PPHMMSignatureTable).toarray()
        all_ucf_genomes["PPHMMLocationTable_coo"] = coo_matrix(PPHMMLocationTable).toarray()
        all_ucf_genomes["GOMSignatureTable_Dict"] = GOMSignatureTable
        pickle.dump(all_ucf_genomes, open(self.fnames['UcfAnnotatorPickle'], "wb"))

    def main(self) -> None:
        '''Generate PPHMM signature table, PPHMM location table, and GOM signature table for unclassified viruses, using PPHMM databases of reference viruses'''
        section_header(
            "Generate PPHMM signature table, PPHMM location table, and GOM signature table\nfor unclassified viruses, using PPHMM databases of reference viruses")
        '''Annotate unclassified viruses using PPHMM & GOM'''
        self.annotate()
