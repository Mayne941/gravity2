import pytest
from app.utils.end_to_end_entrypoint import split_payloads_for_e2e

def test_split_payloads_for_e2e():
    payload = {
        "ExperimentName": "test_experiment",
        "GenomeDescTableFile_FirstPass": "first_pass_desc_table.txt",
        "GenomeDescTableFile_SecondPass": "second_pass_desc_table.txt"
    }

    payload_fp_pl1, payload_fp_pl2, payload_sp_pl1, payload_sp_pl2 = split_payloads_for_e2e(payload)

    firstpass_shelve_dir = "./output/test_experiment_firstpass_pipeline_1"
    secondpass_shelve_dir = "./output/test_experiment_secondpass_pipeline_1"
    firstpass_gen_seq_f = "./output/test_experiment_firstpass_pipeline_1.gb"
    secondpass_gen_seq_f = "./output/test_experiment_secondpass_pipeline_1.gb"

    assert payload_fp_pl1["GenomeDescTableFile"] == "first_pass_desc_table.txt"
    assert payload_fp_pl1["ExpDir"] == firstpass_shelve_dir
    assert payload_fp_pl1["GenomeSeqFile"] == firstpass_gen_seq_f

    assert payload_sp_pl1["GenomeDescTableFile"] == "second_pass_desc_table.txt"
    assert payload_sp_pl1["ExpDir"] == secondpass_shelve_dir
    assert payload_sp_pl1["GenomeSeqFile"] == secondpass_gen_seq_f

    assert payload_fp_pl2["ShelveDir_UcfVirus"] == "output/test_experiment_firstpass_pipeline_2"
    assert payload_fp_pl2["ExpDir_Pl1"] == firstpass_shelve_dir
    assert payload_fp_pl2["GenomeSeqFiles_RefVirus"] == firstpass_gen_seq_f
    assert payload_fp_pl2["GenomeSeqFile_UcfVirus"] == "output/test_experiment_unclassified_viruses.gb"

    assert payload_sp_pl2["ShelveDir_UcfVirus"] == "output/test_experiment_secondpass_pipeline_2"
    assert payload_sp_pl2["ExpDir_Pl1"] == secondpass_shelve_dir
    assert payload_sp_pl2["GenomeSeqFiles_RefVirus"] == secondpass_gen_seq_f
    assert payload_sp_pl2["GenomeSeqFile_UcfVirus"] == "output/test_experiment_unclassified_viruses.gb"
