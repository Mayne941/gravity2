import pickle

from app.utils.error_handlers import raise_gravity_error

def retrieve_genome_vars(read_desc_table_p) -> dict:
    '''Read pickle file to dict containing VMR data'''
    f = pickle.load(open(read_desc_table_p, "rb"))
    if len(f["VirusNameList"]) == 0:
        raise_gravity_error("The list of processed genomes I retrieved is empty!\nThis usually results from errors reading the VMR or you're working with incomplete genomes with AnnotateIncompleteGenomes=false (set to true).")
    return f

def retrieve_pickle(file) -> dict:
    '''Generic retrieve with no error handler'''
    return pickle.load(open(file, "rb"))
