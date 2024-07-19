import os
from fastapi import FastAPI, BackgroundTasks
from fastapi.encoders import jsonable_encoder

from app.pipeline_i import Pipeline_I
from app.pipeline_ii import Pipeline_II
from app.utils.load_premade_pipeline import main as load_premade
from app.utils.banner import print_banner
from app.utils.scrape_vmr import scrape, first_pass_taxon_filter, second_pass
from app.utils.process_fasta import fasta_to_genbank, combine_segments
from app.utils.stdout_utils import progress_msg
from app.utils.end_to_end_entrypoint import split_payloads_for_e2e
from app.utils.api_classes import (Pipeline_i_data, Pipeline_ii_data, Endp_data_scrape_data,
                                   Endp_data_first_pass_taxon_filter, Endp_data_second_pass_filter,
                                   Endp_data_fasta_to_gb, CombineGenomeSegs, E2e_data, Premade_data)

print_banner()
progress_msg("Access the GRAViTy-V2 UI in your browser here (or control-click link):\n\thttp://127.0.0.1:8000/docs")
description = """
    Genome Relationships Applied to Virus Taxonomy is a software framework for identifying and classifying viruses, based on analysis of entire genomes. Developed by Simmonds Group, Peter Medawar Building for Pathogen Research, Nuffield Department of Medicine, University of Oxford.
"""

tags_metadata = [
    {
        "name": "Premade pipelines",
        "description": "Premade all-in-one pipelines for a variety of use cases.",
    },
    {
        "name": "Custom pipelines",
        "description": "Users can fully specify their parameters and run pipeline from any stage.",
    },
    {
        "name": "Pipeline II",
        "description": "PL2 description & user info.",
    },
    {
        "name": "VMR Utilities",
        "description": "Refresh Virus Metadata Resource (VMR) from ICTV and construct datasets."
    },
    {
        "name": "Automated Workflows",
        "description": "Pre-made workflows to automate the transition between a first pass PL2 run to a second pass PL1-2 run."
    },
]

app = FastAPI(
    title="GRAViTy V2",
    version="1.0",
    description=description,
    contact={
        "name": "Rich Mayne (lead developer)",
        "url": "https://www.medawar.ox.ac.uk/research/research-groups/simmonds-group",
        "email": "richard.mayne@ndm.ox.ac.uk"
    },
    license_info={
        "name": "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "url": "https://creativecommons.org/licenses/by-nc/4.0/",
    },
    openapi_tags=tags_metadata
)

'''Premade pipelines'''

@app.post("/similar_viruses/", tags=["Premade pipelines"])
async def similar_viruses(payload: Premade_data):
    payload = process_json(jsonable_encoder(payload) | load_premade("dev/premade_pl_files/similar_viruses"))
    run_pipeline_i_full(payload, refresh_genbank=False)

@app.post("/divergent_viruses/", tags=["Premade pipelines"])
async def similar_viruses(payload: Premade_data):
    payload = process_json(jsonable_encoder(payload) | load_premade("dev/premade_pl_files/divergent_viruses"))
    run_pipeline_i_full(payload, refresh_genbank=False)

@app.post("/long_single_orf_viruses/", tags=["Premade pipelines"])
async def similar_viruses(payload: Premade_data):
    payload = process_json(jsonable_encoder(payload) | load_premade("dev/premade_pl_files/similar_viruses"))
    run_pipeline_i_full(payload, refresh_genbank=False)

'''Dev Endpoints/fns'''

@app.get("/", tags=["Dev Utilities"])
async def read_root():
    return {"response": "API is healthy. Append the current URL to include '/docs/' at the end to visit the GUI."}

def process_json(payload):
    payload = jsonable_encoder(payload)

    try:
        if type(payload["NThreads"]) == str:
            if payload["NThreads"] == "auto":
                payload["N_CPUs"] = os.cpu_count() - 1
            elif payload["NThreads"] == "hpc":
                payload["N_CPUs"] = 1
            else:
                raise Exception("Didn't recognise NThreads argument, defaulting to n-1")
        elif type(payload["NThreads"]) == int:
            if payload["NThreads"] > os.cpu_count():
                raise Exception("You specified to use more CPUs than your computer has. Defauling to n-1")
            payload["N_CPUs"] = payload["NThreads"]
        else:
            raise Exception("Didn't recognise NThreads argument, defaulting to n-1")
    except:
        payload["N_CPUs"] = os.cpu_count() - 1

    return payload

'''Unified Endpoints'''

@app.post("/end_to_end/", tags=["End to end"])
async def end_to_end(payload: E2e_data, background_tasks: BackgroundTasks):
    payload = process_json(payload)
    background_tasks.add_task(run_end_to_end, payload)
    return "Task fired successfully, running in background"

@app.post("/end_to_end_programmatic/", tags=["Dev Utilities"])
async def end_to_end(payload: E2e_data):
    payload = process_json(payload)
    run_end_to_end(payload)
    return True

@app.post("/pipeline_i_programmatic/", tags=["Dev Utilities"])
async def end_to_end_pr(payload: Pipeline_i_data):
    import time
    st = time.time()
    payload = process_json(payload)
    run_pipeline_i_full(payload)
    t = time.time() - st

    import numpy as np
    with open(f'{payload["ExpDir"]}/output/ref_seqs.fasta', "r") as f: seqs = f.readlines()
    seqlen = np.mean([len(i) for i in seqs if not i[0] == ">"])
    with open(f'{payload["ExpDir"]}/BENCHMARKING.txt', "w") as f:
        f.write(f"Time for end-to-end: {t}")
        f.write(f"Mean seq len: {seqlen}")
    print(f"Time for end-to-end: {time.time() - st}")

    return True


def run_end_to_end(payload):
    import time
    st = time.time()
    payload_fp_pl1, payload_fp_pl2, payload_sp_pl1, payload_sp_pl2 = split_payloads_for_e2e(payload)
    if not payload["SkipFirstPass"]:
        run_pipeline_i_full(payload_fp_pl1)
        # run_pipeline_ii_full(payload_fp_pl2)
    run_pipeline_i_full(payload_sp_pl1, refresh_genbank=True) ############
    # run_pipeline_ii_full(payload_sp_pl2, refresh_genbank=True)

    ## RM < TODO BENCHMARKING FOR PAPER
    import numpy as np
    with open(f'{payload_sp_pl2["ExpDir_Pl1"]}/output/ref_seqs.fasta', "r") as f: refs = f.readlines()
    with open(f'{payload_sp_pl2["ExpDir_Pl1"].replace("pipeline_1","pipeline_2")}/output/ref_seqs.fasta', "r") as f: ucfs = f.readlines()
    seqs = refs + ucfs
    seqlen = np.mean([len(i) for i in seqs if not i[0] == ">"])

    with open(f'{payload_sp_pl2["ExpDir_Pl1"].replace("pipeline_1", "pipeline_2")}/BENCHMARKING.txt', "w") as f:
        f.write(f"Time for end-to-end: {time.time() - st}")
        f.write(f"Mean seq len: {seqlen}")
    print(f"Time for end-to-end: {time.time() - st}")
    ##


'''PL1 Entrypoints'''

@app.post("/pipeline_i_full/", tags=["Custom pipelines"])
async def pipeline_i_full(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = process_json(payload)
    background_tasks.add_task(run_pipeline_i_full, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_full(payload, refresh_genbank=False):
    pl = Pipeline_I(payload)
    pl.read_genome_desc_table(refresh_genbank)


@app.post("/pipeline_i_from_pphmmdb_construction/", tags=["Custom pipelines"])
async def pipeline_i__from_pphmmdb_construction(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = process_json(payload)
    background_tasks.add_task(
        run_pipeline_i_from_pphmmdb_construction, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_pphmmdb_construction(payload):
    pl = Pipeline_I(payload)
    pl.pphmmdb_construction()


@app.post("/pipeline_i_from_ref_virus_annotator/", tags=["Custom pipelines"])
async def pipeline_i__from_ref_virus_annotator(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = process_json(payload)
    background_tasks.add_task(run_pipeline_i_from_ref_virus_annotator, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_ref_virus_annotator(payload):
    pl = Pipeline_I(payload)
    pl.ref_virus_annotator()


@app.post("/pipeline_i_from_graph_generator/", tags=["Custom pipelines"])
async def pipeline_i__from_graph_generator(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = process_json(payload)
    background_tasks.add_task(run_pipeline_i_from_graph_generator, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_graph_generator(payload):
    pl = Pipeline_I(payload)
    pl.make_graphs()


@app.post("/pipeline_i_from_mutual_info_calculator/", tags=["Custom pipelines"])
async def pipeline_i_from_mutual_info_calculator(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = process_json(payload)
    background_tasks.add_task(
        run_pipeline_i_from_mutual_info_calculator, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_mutual_info_calculator(payload):
    pl = Pipeline_I(payload)
    pl.mutual_info_calculator()


'''PL2 entrypoints'''


# @app.post("/pipeline_ii_full/", tags=["Pipeline II"])
# async def pipeline_ii_full(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
#     payload = process_json(payload)
#     background_tasks.add_task(run_pipeline_ii_full, payload)
#     return "Task fired successfully, running in background"


# def run_pipeline_ii_full(payload, refresh_genbank=False):
#     pl = Pipeline_II(payload)
#     pl.read_genome_desc_table(refresh_genbank)


# @app.post("/pipeline_ii_from_pphmmdb_construction/", tags=["Pipeline II"])
# async def pipeline_ii_from_pphmmdb_construction(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
#     payload = process_json(payload)
#     background_tasks.add_task(
#         run_pipeline_ii_from_pphmmdb_construction, payload)
#     return "Task fired successfully, running in background"


# def run_pipeline_ii_from_pphmmdb_construction(payload):
#     pl = Pipeline_II(payload)
#     pl.pphmmdb_construction()


# @app.post("/pipeline_ii_from_ucf_virus_annotator/", tags=["Pipeline II"])
# async def pipeline_ii_from_ucf_virus_annotator(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
#     payload = process_json(payload)
#     background_tasks.add_task(
#         run_pipeline_ii_from_ucf_virus_annotator, payload)
#     return "Task fired successfully, running in background"


# def run_pipeline_ii_from_ucf_virus_annotator(payload):
#     pl = Pipeline_II(payload)
#     pl.ucf_virus_annotator()


# @app.post("/pipeline_ii_from_virus_classification/", tags=["Pipeline II"])
# async def pipeline_ii_from_virus_classification(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
#     payload = process_json(payload)
#     background_tasks.add_task(
#         run_pipeline_ii_from_virus_classification, payload)
#     return "Task fired successfully, running in background"


# def run_pipeline_ii_from_virus_classification(payload):
#     pl = Pipeline_II(payload)
#     pl.virus_classification()


'''Utility entrypoints'''


@app.post("/scrape_vmr/", tags=["VMR Utilities"])
async def run_vmr_scrape(trigger: Endp_data_scrape_data):
    payload = process_json(trigger)
    return scrape(payload)

@app.post("/construct_first_pass_vmr_taxon_filter/", tags=["VMR Utilities"])
async def vmr_first_pass_taxon_filter(trigger: Endp_data_first_pass_taxon_filter):
    payload = process_json(trigger)
    return first_pass_taxon_filter(payload)

@app.post("/construct_second_pass_vmr/", tags=["VMR Utilities"])
async def vmr_second_pass(trigger: Endp_data_second_pass_filter):
    payload = process_json(trigger)
    return second_pass(payload)

@app.post("/convert_fasta_to_genbank_and_vmr/", tags=["VMR Utilities"])
async def convert_fasta_to_genbank(trigger: Endp_data_fasta_to_gb):
    payload = process_json(trigger)
    return fasta_to_genbank(payload)

@app.post("/combine_genome_segments/", tags=["VMR Utilities"])
async def convert_genome_segments(trigger: CombineGenomeSegs):
    payload = process_json(trigger)
    return combine_segments(payload)

'''Deprecated Endpoints'''
# Baltimore filter
# @app.post("/construct_first_pass_vmr_baltimore_filter/", tags=["VMR Utilities"])
# async def vmr_first_pass_baltimore_filter(trigger: FirstPassBaltimoreFilter):
#     payload = process_json(trigger)
#     return first_pass_baltimore_filter(payload)

# Generic filter
# @app.post("/filter_vmr/", tags=["VMR Utilities"])
# async def filter_that_vmr(trigger: VmrFilter):
#     payload = process_json(trigger)
#     return vmr_filter(payload)
