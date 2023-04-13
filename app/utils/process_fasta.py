from Bio import SeqIO
import pandas as pd

def get_vmr_cols():
    return ["Unnamed: 0", "Sort", "Realm", "Subrealm", "Kingdom", "Subkingdom",
            "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", "Family",
            "Subfamily", "Genus", "Subgenus", "Species", "Exemplar or additional isolate",
            "Virus name(s)", "Virus name abbreviation(s)", "Virus isolate designation",
            "Virus GENBANK accession", "Genome coverage", "Genome composition", "Host source",
            "Baltimore Group", "Genetic code table", "Taxonomic grouping"]

def fasta_to_genbank(payload):
    '''Convert FASTA-style base sequences into genbank files, generate VMR-like object'''
    '''Parse input, prepare data structures'''
    input_handle = open(
        f"{payload['save_path']}/{payload['fasta_fname']}", "r")
    # output_handle = open(
    #     f"{payload['save_path']}/{payload['genbank_fname']}", "w")
    sequences = list(SeqIO.parse(input_handle, "fasta"))
    seq_names, seq_codes = [], []

    for i in range(len(sequences)):
        sequences[i].annotations['molecule_type'] = 'DNA'
        seq_codes.append(sequences[i].id)
        seq_names.append(f"Query_{i+1}_{sequences[i].id}")
        sequences[i].id = f"Query_{i+1}"

    # RM < DISABLED AS OUTPUT FILE SEEMS INCOMPATIBLE WITH GRAVITY
    #count = SeqIO.write(sequences, output_handle, "genbank")

    '''Build VMR-like csv'''
    df = pd.DataFrame(columns=get_vmr_cols())
    df["Virus GENBANK accession"] = seq_codes
    df["Virus name(s)"] = seq_names
    df["Genetic code table"] = 1

    df.to_csv(f"{payload['save_path']}/{payload['vmr_fname']}")

    # output_handle.close()
    input_handle.close()
    return f"Successfully converted {df.shape[0]} records"

def combine_segments(payload):
    input_handle = open(
        f"{payload['save_path']}/{payload['fasta_fname']}", "r")
    output_handle = open(
        f"{payload['save_path']}/{payload['genbank_fname']}", "w")
    sequences = list(SeqIO.parse(input_handle, "fasta"))
    seq_names, seq_codes = [], []

    for i in range(len(sequences)):
        sequences[i].annotations['molecule_type'] = 'DNA'
        seq_codes.append(sequences[i].id)
        seq_names.append(f"Query_{i+1}_{sequences[i].id}")
        sequences[i].id = f"Query_{i+1}"

    count = SeqIO.write(sequences, output_handle, "genbank")

    '''Build VMR-like csv'''
    df = pd.DataFrame(columns=get_vmr_cols())
    df["Virus GENBANK accession"] = [", ".join(seq_codes).replace(",", ";")]
    df["Virus name(s)"] = f"Query_{payload['combined_seq_name']}"
    df["Genetic code table"] = 1
    df.to_csv(f"{payload['save_path']}/{payload['vmr_fname']}")

    output_handle.close()
    input_handle.close()
    return "Successfully combined {} records".format(count)
