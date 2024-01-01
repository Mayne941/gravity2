def split_payloads_for_e2e(payload):
    '''Construct input parameter payloads for each pipeline per pass, for end to end entrypoint'''
    payload_fp_pl1, payload_fp_pl2, payload_sp_pl1, payload_sp_pl2 = payload.copy(), payload.copy(), payload.copy(), payload.copy()

    firstpass_shelve_dir  = f"./output/{payload['ExperimentName']}_firstpass_pipeline_1"
    secondpass_shelve_dir = f"./output/{payload['ExperimentName']}_secondpass_pipeline_1"
    firstpass_gen_seq_f   = f"./output/{payload['ExperimentName']}_firstpass_pipeline_1.gb"
    secondpass_gen_seq_f  = f"./output/{payload['ExperimentName']}_secondpass_pipeline_1.gb"

    '''FP, PL1'''
    payload_fp_pl1["GenomeDescTableFile"] = payload["GenomeDescTableFile_FirstPass"]
    payload_fp_pl1["ExpDir"]              = firstpass_shelve_dir
    payload_fp_pl1["TaxoGroupingFile"]    = payload["TaxoGroupingFile_FirstPass"]
    payload_fp_pl1["GenomeSeqFile"]       = firstpass_gen_seq_f

    '''SP, PL1'''
    payload_sp_pl1["GenomeDescTableFile"] = payload["GenomeDescTableFile_SecondPass"]
    payload_sp_pl1["ExpDir"]              = secondpass_shelve_dir
    payload_sp_pl1["TaxoGroupingFile"]    = payload["TaxoGroupingFile_SecondPass"]
    payload_sp_pl1["GenomeSeqFile"]       = secondpass_gen_seq_f

    '''FP, PL2'''
    payload_fp_pl2["ShelveDir_UcfVirus"]      = f"output/{payload['ExperimentName']}_firstpass_pipeline_2"
    payload_fp_pl2["ExpDir_Pl1"]     = firstpass_shelve_dir
    payload_fp_pl2["GenomeSeqFiles_RefVirus"] = firstpass_gen_seq_f
    payload_fp_pl2["GenomeSeqFile_UcfVirus"]  = f"output/{payload['ExperimentName']}_unclassified_viruses.gb"

    '''SP, PL2'''
    payload_sp_pl2["ShelveDir_UcfVirus"]      = f"output/{payload['ExperimentName']}_secondpass_pipeline_2"
    payload_sp_pl2["ExpDir_Pl1"]     = secondpass_shelve_dir
    payload_sp_pl2["GenomeSeqFiles_RefVirus"] = secondpass_gen_seq_f
    payload_sp_pl2["GenomeSeqFile_UcfVirus"]  = f"output/{payload['ExperimentName']}_unclassified_viruses.gb"

    return payload_fp_pl1, payload_fp_pl2, payload_sp_pl1, payload_sp_pl2
