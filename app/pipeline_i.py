
import os
import time
import json
import shutil

from app.src.read_genome_desc_table import ReadGenomeDescTable
from app.src.pphmmdb_construction import PPHMMDBConstruction
from app.src.ref_virus_annotator import RefVirusAnnotator
from app.src.pl1_graphs import GRAViTyDendrogramAndHeatmapConstruction
from app.src.mutual_information_calculator import MutualInformationCalculator
from app.utils.generate_logs import Log_Generator_Pl1
from app.utils.timer import timing
from app.utils.error_handlers import raise_gravity_error


class Pipeline_I:
    def __init__(self, payload):
        self.options = payload
        '''Create directories'''
        if not os.path.exists(self.options['ExpDir']):
            os.makedirs(self.options['ExpDir'])

        '''Logs'''
        self.log_gen = Log_Generator_Pl1(
            self.options, self.options['ExpDir'])
        self.logs = self.log_gen.entrypoint()
        self.actual_start = time.time()
        with open(f"{self.options['ExpDir']}/run_parameters.json", "w") as f:
            f.write(json.dumps(payload))
        try:
            shutil.copyfile(self.options['GenomeDescTableFile'], f"{self.options['ExpDir']}/PL1_vmr.csv")
        except: raise_gravity_error(f"I couldn't find your Pipeline I VMR-like document.\n"
                                    f"This usually happens when you've copied settings from a previous run. Try deleting your .gb fine and start again.")

    @timing
    def read_genome_desc_table(self, refresh_genbank):
        '''I: Fire Read Genom Desc Table'''
        [print(log_text) for log_text in self.logs[0]]
        rgdt = ReadGenomeDescTable(
            payload=self.options,
            GenomeDescTableFile=self.options['GenomeDescTableFile'],
            GenomeSeqFile=self.options['GenomeSeqFile'],
            ExpDir=self.options['ExpDir'],
            RefreshGenbank=refresh_genbank
        )
        rgdt.entrypoint()
        self.pphmmdb_construction()

    @timing
    def pphmmdb_construction(self):
        '''II: Fire PPHMDB Constructor'''
        [print(log_text) for log_text in self.logs[1]]
        pph = PPHMMDBConstruction(
            payload=self.options,
            GenomeSeqFile=self.options['GenomeSeqFile'],
            ExpDir=self.options['ExpDir'],
        )
        pph.main()
        self.ref_virus_annotator()

    @timing
    def ref_virus_annotator(self):
        '''III: Fire Reference Virus Annotator'''
        [print(log_text) for log_text in self.logs[2]]
        rva = RefVirusAnnotator(
            payload=self.options,
            GenomeSeqFile=self.options['GenomeSeqFile'],
            ExpDir=self.options['ExpDir'],
        )
        rva.main()
        self.make_graphs()

    @timing
    def make_graphs(self):
        '''IV: Fire Heatmap and Dendrogram Constructors'''
        [print(log_text) for log_text in self.logs[3]]
        ghm = GRAViTyDendrogramAndHeatmapConstruction(
            payload=self.options,
            ExpDir=self.options['ExpDir'],
        )
        ghm.main()
        self.mutual_info_calculator()

    @timing
    def mutual_info_calculator(self):
        '''V: Fire Mutal Info Calculator'''
        [print(log_text) for log_text in self.logs[4]]
        mic = MutualInformationCalculator(
            payload=self.options,
            ExpDir=self.options['ExpDir'],
        )
        mic.main()

        print(f"Time to complete: {time.time() - self.actual_start}")


if __name__ == '__main__':
    pl = Pipeline_I()
    pl.main()
