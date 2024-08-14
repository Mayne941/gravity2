from app.src.read_genome_desc_table import ReadGenomeDescTable
from app.src.pphmmdb_construction import PPHMMDBConstruction
from app.src.ucf_virus_annotator import UcfVirusAnnotator
from app.src.virus_classification import VirusClassificationAndEvaluation

from app.utils.str_to_bool import str2bool
from app.utils.timer import timing
from app.utils.generate_logs import Log_Generator_Pl2

import optparse
import os
import time
import json
import shutil

class Pipeline_II:
    def __init__(self, payload):
        self.options = payload
        '''Create directories'''
        if not os.path.exists(self.options["ShelveDir_UcfVirus"]):
            os.makedirs(self.options["ShelveDir_UcfVirus"])
        '''Logs'''
        self.log_gen = Log_Generator_Pl2(
            self.options, self.options["ShelveDir_UcfVirus"])
        self.actual_start = time.time()
        self.logs = self.log_gen.entrypoint()
        shutil.copyfile(self.options['GenomeDescTableFile_UcfVirus'], f"{self.options['ShelveDir_UcfVirus']}/PL2_vmr.csv")
        '''Seed payload with dummies so as to not hit PL1-specific options'''
        self.options["TaxoGroupingFile"] = None
        '''Save logs'''
        with open(f"{self.options['ShelveDir_UcfVirus']}/run_parameters.json", "w") as f:
            f.write(json.dumps(payload))

    @timing
    def read_genome_desc_table(self, refresh_genbank):
        '''I: Fire Read Genom Desc Table'''
        [print(log_text) for log_text in self.logs[0]]
        rgdt = ReadGenomeDescTable(
            payload=self.options,
            GenomeDescTableFile=self.options['GenomeDescTableFile_UcfVirus'],
            GenomeSeqFile=self.options['GenomeSeqFile_UcfVirus'],
            ExpDir=self.options['ShelveDir_UcfVirus'],
            RefreshGenbank=refresh_genbank,
            is_secondpass=True
        )
        rgdt.entrypoint()
        self.pphmmdb_construction()

    @timing
    def pphmmdb_construction(self):
        '''II: Fire PPHMDB Constructor'''
        if str2bool(self.options['UseUcfVirusPPHMMs']) == True:
            [print(log_text) for log_text in self.logs[1]]
            pc = PPHMMDBConstruction(
                payload=self.options,
                GenomeSeqFile=self.options['GenomeSeqFile_UcfVirus'],
                ExpDir=self.options['ShelveDir_UcfVirus'],
            )
            pc.main()
        self.ucf_virus_annotator()

    @timing
    def ucf_virus_annotator(self):
        '''III: Fire UCF Virus Annotator'''
        [print(log_text) for log_text in self.logs[2]]
        ucf = UcfVirusAnnotator(
            payload=self.options,
            ExpDir=self.options['ShelveDir_UcfVirus'], # TODO Is it necessary to pass this as it's in the payload??
        )
        ucf.main()
        self.virus_classification()

    @timing
    def virus_classification(self):
        '''IV: Fire Classifiers'''
        [print(log_text) for log_text in self.logs[3]]
        vce = VirusClassificationAndEvaluation(
            payload=self.options,
            ExpDir=self.options['ShelveDir_UcfVirus'],
        )
        vce.main()

        print(f"Time to complete: {time.time() - self.actual_start}")


if __name__ == '__main__':
    pl = Pipeline_II()
    pl.main()
