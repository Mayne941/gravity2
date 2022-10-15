from Bio import SeqIO
import pandas as pd


def fasta_to_genbank(payload):
    '''Convert FASTA-style base sequences into genbank files, generate VMR-like object'''
    '''Parse input, prepare data structures'''
    input_handle = open(
        f"{payload['save_path']}/{payload['fasta_fname']}", "r")
    output_handle = open(
        f"{payload['save_path']}/{payload['genbank_fname']}", "w")
    sequences = list(SeqIO.parse(input_handle, "fasta"))
    seq_names, seq_codes = [], []

    for i in range(len(sequences)):
        sequences[i].annotations['molecule_type'] = 'DNA'
        seq_names.append(sequences[i].id)
        seq_codes.append(f"GRAVITY_{i+1}")
        sequences[i].id = f"GRAVITY_{i+1}"

    count = SeqIO.write(sequences, output_handle, "genbank")

    '''Build VMR-like csv'''
    cols = ["Unnamed: 0", "Sort", "Realm", "Subrealm", "Kingdom", "Subkingdom",
            "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", "Family",
            "Subfamily", "Genus", "Subgenus", "Species", "Exemplar or additional isolate",
            "Virus name(s)", "Virus name abbreviation(s)", "Virus isolate designation",
            "Virus GENBANK accession", "Genome coverage", "Genome composition", "Host source",
            "Baltimore Group", "Genetic code table", "Taxonomic grouping"]
    df = pd.DataFrame(columns=cols)
    df["Virus GENBANK accession"] = seq_codes
    df["Virus name(s)"] = seq_names
    df["Genetic code table"] = "Complete coding genome"
    df.to_csv(f"{payload['save_path']}/{payload['vmr_fname']}")

    output_handle.close()
    input_handle.close()
    print("Successfully converted {} records".format(count))
