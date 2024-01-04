from termcolor import colored

def raise_gravity_error(msg):
    raise SystemExit(f"{colored('FATAL ERROR', 'red')}: {msg}")

def error_handler_virus_classifier(payload):
    if payload['UseUcfVirusPPHMMs'] == True:
        if payload['GenomeSeqFiles_RefVirus'] == None or payload['GenomeSeqFile_UcfVirus'] == None:
            raise ValueError(
                "'UseUcfVirusPPHMMs' == True. 'GenomeSeqFiles_RefVirus' and 'GenomeSeqFile_UcfVirus' cannot be None")
