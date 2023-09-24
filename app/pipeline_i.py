
import optparse
import os
import time
import json
import shutil

from app.src.read_genome_desc_table import ReadGenomeDescTable
from app.src.pphmmdb_construction import PPHMMDBConstruction
from app.src.ref_virus_annotator import RefVirusAnnotator
from app.src.graphing_tools import GRAViTyDendrogramAndHeatmapConstruction
from app.src.mutual_information_calculator import MutualInformationCalculator
from app.utils.str_to_bool import str2bool
from app.utils.generate_logs import Log_Generator_Pl1
from app.utils.timer import timing


class Pipeline_I:
    def __init__(self, payload):
        self.options = payload
        '''Create directories'''
        if not os.path.exists(self.options['ShelveDir']):
            os.makedirs(self.options['ShelveDir'])

        '''Logs'''
        self.log_gen = Log_Generator_Pl1(
            self.options, self.options['ShelveDir'])
        self.logs = self.log_gen.entrypoint()
        self.actual_start = time.time()
        with open(f"{self.options['ShelveDir']}/run_parameters.json", "w") as f:
            f.write(json.dumps(payload))
        shutil.copyfile(self.options['GenomeDescTableFile'], f"{self.options['ShelveDir']}/PL1_vmr.csv")

        '''Catch bad database flags'''
        if (self.options['Database'] != None and self.options['Database_Header'] == None):
            raise optparse.OptionValueError(
                f"You have specified DATABASE as {self.options['Database']}, 'Database_Header' cannot be 'None'")
        if (self.options['Database'] == None and self.options['Database_Header'] != None):
            Proceed = input(
                f"You have specified 'Database_Header' as {self.options['Database_Header']}, but 'Database' is 'None'. GRAViTy will analyse all genomes. Do you want to proceed? [Y/n]: ")
            if Proceed != "Y":
                raise SystemExit("GRAViTy terminated.")

    @timing
    def read_genome_desc_table(self, refresh_genbank):
        '''I: Fire Read Genom Desc Table'''
        [print(log_text) for log_text in self.logs[0]]
        rgdt = ReadGenomeDescTable(
            GenomeDescTableFile=self.options['GenomeDescTableFile'],
            GenomeSeqFile=self.options['GenomeSeqFile'],
            ShelveDir=self.options['ShelveDir'],
            Database=self.options['Database'],
            Database_Header=self.options['Database_Header'],
            TaxoGrouping_Header=self.options['TaxoGrouping_Header'],
            TaxoGroupingFile=self.options['TaxoGroupingFile'],
            RefreshGenbank=refresh_genbank
        )
        rgdt.entrypoint()
        self.pphmmdb_construction()

    @timing
    def pphmmdb_construction(self):
        '''II: Fire PPHMDB Constructor'''
        [print(log_text) for log_text in self.logs[1]]
        pph = PPHMMDBConstruction(
            genbank_email=self.options['genbank_email'],
            GenomeSeqFile=self.options['GenomeSeqFile'],
            ShelveDir=self.options['ShelveDir'],
            ProteinLength_Cutoff=self.options['ProteinLength_Cutoff'],
            IncludeIncompleteGenomes=str2bool(
                self.options['IncludeProteinsFromIncompleteGenomes']),
            BLASTp_evalue_Cutoff=self.options['BLASTp_evalue_Cutoff'],
            BLASTp_PercentageIden_Cutoff=self.options['BLASTp_PercentageIden_Cutoff'],
            BLASTp_QueryCoverage_Cutoff=self.options['BLASTp_QueryCoverage_Cutoff'],
            BLASTp_SubjectCoverage_Cutoff=self.options['BLASTp_SubjectCoverage_Cutoff'],
            BLASTp_num_alignments=self.options['BLASTp_num_alignments'],
            BLASTp_N_CPUs=self.options['N_CPUs'],
            MUSCLE_GapOpenCost=self.options['MUSCLE_GapOpenCost'],
            MUSCLE_GapExtendCost=self.options['MUSCLE_GapExtendCost'],
            ProtClustering_MCLInflation=self.options['ProtClustering_MCLInflation'],
            N_AlignmentMerging=self.options['N_AlignmentMerging'],
            HHsuite_evalue_Cutoff=self.options['HHsuite_evalue_Cutoff'],
            HHsuite_pvalue_Cutoff=self.options['HHsuite_pvalue_Cutoff'],
            HHsuite_N_CPUs=self.options['N_CPUs'],
            HHsuite_QueryCoverage_Cutoff=self.options['HHsuite_QueryCoverage_Cutoff'],
            HHsuite_SubjectCoverage_Cutoff=self.options['HHsuite_SubjectCoverage_Cutoff'],
            PPHMMClustering_MCLInflation=self.options['PPHMMClustering_MCLInflation_ForAlnMerging'],
            HMMER_PPHMMDb_ForEachRoundOfPPHMMMerging=str2bool(
                self.options['HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging']),
        )
        pph.main()
        self.ref_virus_annotator()

    @timing
    def ref_virus_annotator(self):
        '''III: Fire Reference Virus Annotator'''
        [print(log_text) for log_text in self.logs[2]]
        rva = RefVirusAnnotator(
            GenomeSeqFile=self.options['GenomeSeqFile'],
            ShelveDir=self.options['ShelveDir'],
            ProteinLength_Cutoff=self.options['ProteinLength_Cutoff'],
            IncludeIncompleteGenomes=str2bool(
                self.options['AnnotateIncompleteGenomes']),
            HMMER_N_CPUs=self.options['N_CPUs'],
            HMMER_C_EValue_Cutoff=self.options['HMMER_C_EValue_Cutoff'],
            HMMER_HitScore_Cutoff=self.options['HMMER_HitScore_Cutoff'],
            RemoveSingletonPPHMMs=str2bool(
                self.options['RemoveSingletonPPHMMs']),
            N_VirusesOfTheClassToIgnore=self.options['N_VirusesOfTheClassToIgnore'],
            PPHMMSorting=str2bool(self.options['PPHMMSorting']),
            HHsuite_evalue_Cutoff=self.options['HHsuite_evalue_Cutoff'],
            HHsuite_pvalue_Cutoff=self.options['HHsuite_pvalue_Cutoff'],
            HHsuite_N_CPUs=self.options['N_CPUs'],
            HHsuite_QueryCoverage_Cutoff=self.options['HHsuite_QueryCoverage_Cutoff'],
            HHsuite_SubjectCoverage_Cutoff=self.options['HHsuite_SubjectCoverage_Cutoff'],
            PPHMMClustering_MCLInflation=self.options['PPHMMClustering_MCLInflation_ForPPHMMSorting'],
        )
        rva.main()
        self.make_graphs()

    @timing
    def make_graphs(self):
        '''IV: Fire Heatmap and Dendrogram Constructors'''
        [print(log_text) for log_text in self.logs[3]]
        ghm = GRAViTyDendrogramAndHeatmapConstruction(
            ShelveDir=self.options['ShelveDir'],
            IncludeIncompleteGenomes=str2bool(
                self.options['AnnotateIncompleteGenomes']),
            SimilarityMeasurementScheme=self.options['SimilarityMeasurementScheme'],
            p=self.options['p'],
            Dendrogram=str2bool(self.options['Dendrogram']),
            Dendrogram_LinkageMethod=self.options['Dendrogram_LinkageMethod'],
            Bootstrap=str2bool(self.options['Bootstrap']),
            N_Bootstrap=self.options['N_Bootstrap'],
            Bootstrap_method=self.options['Bootstrap_method'],
            Bootstrap_N_CPUs=self.options['N_CPUs'],
            Heatmap=str2bool(self.options['Heatmap']),
            Heatmap_VirusOrderScheme=self.options['Heatmap_VirusOrderScheme'],
            Heatmap_WithDendrogram=str2bool(
                self.options['Heatmap_WithDendrogram']),
            Heatmap_DendrogramFile=self.options['Heatmap_DendrogramFile'],
            Heatmap_DendrogramSupport_Cutoff=self.options['Heatmap_DendrogramSupport_Cutoff'],
            VirusGrouping=str2bool(self.options['VirusGrouping']),
        )
        ghm.main()
        self.mutual_info_calculator()

    @timing
    def mutual_info_calculator(self):
        '''V: Fire Mutal Info Calculator'''
        [print(log_text) for log_text in self.logs[4]]
        mic = MutualInformationCalculator(
            ShelveDir=self.options['ShelveDir'],
            IncludeIncompleteGenomes=str2bool(
                self.options['AnnotateIncompleteGenomes']),
            VirusGroupingFile=self.options['VirusGroupingFile'],
            N_Sampling=self.options['N_Sampling'],
            SamplingStrategy=self.options['SamplingStrategy'],
            SampleSizePerGroup=self.options['SampleSizePerGroup'],
        )
        mic.main()

        print(f"Time to complete: {time.time() - self.actual_start}")


if __name__ == '__main__':
    pl = Pipeline_I()
    pl.main()
