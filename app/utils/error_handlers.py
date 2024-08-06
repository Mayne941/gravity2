from termcolor import colored

def raise_gravity_error(msg):
    raise SystemExit(f"{colored('FATAL ERROR', 'red')}: {msg}")

def raise_gravity_warning(msg):
    print(f"{colored('WARNING', 'yellow')}: {msg}")

def error_handler_virus_classifier(payload):
    if payload['UseUcfVirusPPHMMs'] == True:
        if payload['GenomeSeqFiles_RefVirus'] == None or payload['GenomeSeqFile_UcfVirus'] == None:
            raise ValueError(
                "'UseUcfVirusPPHMMs' == True. 'GenomeSeqFiles_RefVirus' and 'GenomeSeqFile_UcfVirus' cannot be None")

def main_error_msg(out):
    return f"Critical error in command line call to {out}."

def decode_stdout(out):
    return out.decode("utf-8")

def error_handler_hmmbuild(out, msg):
    out = decode_stdout(out)
    if "input parse error" in out:
        raise_gravity_error(f"{main_error_msg(msg)}."
                            f"Ensure HMMER 3 is installed and active in your environment."
                            f"Check that input fasta files (Mash/Clusters/Cluster_x.fasta) are being generated and are not empty.")

def error_handle_mafft(out, msg):
    out = decode_stdout(out)
    if "killed" in out or not "done" in out:
        raise_gravity_error(f"{main_error_msg(msg)}."
                            f"Ensure mafft is installed and active in your environment."
                            f"Mafft is NOT memory safe: ensure you're not overflowing your memory (if you are, try reducing dataset size).")

def error_handler_mash_sketch(out, msg):
    out = decode_stdout(out)
    if not "Sketching" in out and not "Writing" in out:
        raise_gravity_error(f"{main_error_msg(msg)}."
                            f"Ensure Mash is installed and active in your environment.")

def error_handler_mash_dist(out, msg):
    out = decode_stdout(out)
    if out != "":
        if "mash <command> [options]" in out:
            '''For dependency test - empty call contents'''
            return
        else:
            raise_gravity_error(f"{main_error_msg(msg)}."
                                f"Ensure Mash is installed and active in your environment.")

def error_handler_hmmscan(out, msg):
    out = decode_stdout(out)
    if "exception" in out or not "[ok]" in out:
        if "Target sequence length > 100K, over comparison pipeline limit." in out:
            raise_gravity_error(f"{main_error_msg(msg)}."
                                f"You've tried to pass in a proteins sequence of length > 100k AAs, which the software can't handle."
                                f"There are a few viruses with very long ORFs (consider subsetting genomes), otherwise your input sequences might be incorrect.")
        else:
            raise_gravity_error(f"{main_error_msg(msg)}."
                                f"Ensure HMMER 3 is installed and active in your environment."
                                f"Check that input .hmm files are being generatedby hmmpress, in PPHMMDB construction step.")

def error_handler_blast(out, msg):
    out = decode_stdout(out)
    if not out == "":
        if "BLAST query/options error:" in out:
            '''For dependency test - empty call contents'''
            return
        else:
            raise_gravity_error(f"{main_error_msg(msg)}."
                                f"Ensure BLAST is installed and active in your environment.")

def error_handler_mcl(out, msg):
    out = decode_stdout(out)
    if not "[mcl] new tab created" in out:
        if "[mcl] usage:" in out:
            '''For dependency test - empty call contents'''
            return
        else:
            raise_gravity_error(f"{main_error_msg(msg)}."
                                f"Ensure Mcl is installed and active in your environment.")
