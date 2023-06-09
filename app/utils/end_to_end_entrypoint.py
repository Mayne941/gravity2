def split_payloads_for_e2e(payload):
    '''Construct input parameter payloads for each pipeline per pass, for end to end entrypoint'''
    payload_fp_pl1, payload_fp_pl2, payload_sp_pl1, payload_sp_pl2 = payload.copy(), payload.copy(), payload.copy(), payload.copy()

    '''FP, PL1'''
    payload_fp_pl1["GenomeDescTableFile"] = payload["GenomeDescTableFile_FirstPass"]
    payload_fp_pl1["ShelveDir"]           = payload["ShelveDir_FirstPass"]
    payload_fp_pl1["TaxoGroupingFile"]    = payload["TaxoGroupingFile_FirstPass"]
    payload_fp_pl1["GenomeSeqFile"]       = payload["GenomeSeqFile_FirstPass"]

    '''SP, PL1'''
    payload_sp_pl1["GenomeDescTableFile"] = payload["GenomeDescTableFile_SecondPass"]
    payload_sp_pl1["ShelveDir"]           = payload["ShelveDir_SecondPass"]
    payload_sp_pl1["TaxoGroupingFile"]    = payload["TaxoGroupingFile_SecondPass"]
    payload_sp_pl1["GenomeSeqFile"]       = payload["GenomeSeqFile_SecondPass"]

    '''FP, PL2'''
    payload_fp_pl2["ShelveDir_UcfVirus"]      = payload["ShelveDir_UcfVirus_FirstPass"]
    payload_fp_pl2["ShelveDirs_RefVirus"]     = payload["ShelveDir_FirstPass"]
    payload_fp_pl2["GenomeSeqFiles_RefVirus"] = payload["GenomeSeqFile_FirstPass"]

    '''SP, PL2'''
    payload_sp_pl2["ShelveDir_UcfVirus"]      = payload["ShelveDir_UcfVirus_SecondPass"]
    payload_sp_pl2["ShelveDirs_RefVirus"]     = payload["ShelveDir_SecondPass"]
    payload_sp_pl2["GenomeSeqFiles_RefVirus"] = payload["GenomeSeqFile_SecondPass"]

    return payload_fp_pl1, payload_fp_pl2, payload_sp_pl1, payload_sp_pl2
