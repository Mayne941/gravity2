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
        f"{payload['fasta_fname']}", "r")
    output_handle = open(
        f"{payload['genbank_fname']}", "w")
    raw_seqs = SeqIO.parse(input_handle, "fasta")
    sequences = list(raw_seqs)
    seq_names, seq_codes, output_seqs = [], [], []

    for i in range(len(sequences)):
        '''Annotate each sequence and create VMR entries'''
        sequences[i].annotations['molecule_type'] = 'DNA'
        if not sequences[i].id in seq_codes:
            '''Don't process duplicate Acc IDs - omit all but first'''
            output_seqs.append(sequences[i])
            seq_codes.append(sequences[i].id)
            seq_names.append(f"Query_{i+1}_{sequences[i].id}")
            sequences[i].id = sequences[i].id
        else:
            print(f"WARNING: I detected an entry in your FASTA file that has a duplicate accession ID: {sequences[id].id}\n GRAViTy has only kept the first entry!")

    '''Write genbank'''
    _ = SeqIO.write(output_seqs, output_handle, "genbank")
    output_handle.close()
    input_handle.close()

    '''Build VMR-like csv'''
    df = pd.DataFrame(columns=get_vmr_cols())
    df["Virus GENBANK accession"] = seq_codes
    df["Virus name(s)"] = seq_names
    df["Genetic code table"] = 1

    df.to_csv(f"{payload['vmr_fname']}")

    return f"Successfully converted {df.shape[0]} records"

def combine_segments(payload):
    input_handle = open(
        f"{payload['fasta_fname']}", "r")
    output_handle = open(
        f"{payload['genbank_fname']}", "w")
    sequences = list(SeqIO.parse(input_handle, "fasta"))
    seq_names, seq_codes = [], []


    '''Here is the spreadsheet I have been using for the taxonomy proposal comparisons.
       These have all been run and it includes several segmented viruses that were joined
       by a script I wrote (eg. 2022.018M.A.Orthobunyavirus_29nsp_abolish4sp). If you could test of these, we could compare directly with the outputs run previously. '''
    '''Note that the program would be best if it could process several segmented virus sequences at the same time (as in the list),
      rather than individually. As mentioned  there needs to be some kind of notation so that accession numbers for segments are separated
      from each using a different delimiter than between those of different viruses. I suppose the other issue is working out how to generate
      similar sets of concatenated sequences when they are not on GenBank. There would have to be some way to indicate which sequences in a
      combined FASTA file should be joined up and which should be in separate sets.'''
    # RM <TODO WIP
    # len_and_ids = sorted((len(rec), rec.id)
    # for rec in SeqIO.parse("mycobt.fasta","fasta"))
    # ids = reversed([id for (length, id) in len_and_ids])
    # del len_and_ids
    # record_index = SeqIO.index("mycobt.fasta", "fasta")
    # records = (record_index[id] for id in ids)
    # SeqIO.write(records, "sorted.fasta", "fasta")

    for i in range(len(sequences)):
        sequences[i].annotations['molecule_type'] = 'DNA'
        seq_codes.append(sequences[i].id)
        seq_names.append(f"Query_{i+1}_{sequences[i].id}")
        sequences[i].id = sequences[i].id

    count = SeqIO.write(sequences, output_handle, "genbank")

    '''Build VMR-like csv'''
    df = pd.DataFrame(columns=get_vmr_cols())
    df["Virus GENBANK accession"] = [", ".join(seq_codes).replace(",", ";")]
    df["Virus name(s)"] = f"Query_{payload['combined_seq_name']}"
    df["Genetic code table"] = 1
    df.to_csv(f"{payload['vmr_fname']}")

    output_handle.close()
    input_handle.close()
    return f"Successfully combined {count} records"
