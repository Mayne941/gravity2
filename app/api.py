from app.pipeline_i import Pipeline_I
from app.pipeline_ii import Pipeline_II
from app.utils.api_classes import Pipeline_i_data, Pipeline_ii_data, ScrapeData, FirstPass, SecondPass
from app.utils.scrape_vmr import scrape, first_pass, second_pass

from fastapi import FastAPI, BackgroundTasks
from fastapi.encoders import jsonable_encoder

description = """
    App in development
"""

tags_metadata = [
    {
        "name": "Pipeline I",
        "description": "PL1 description & user info",
    },
    {
        "name": "Pipeline II",
        "description": "PL2 description & user info.",
    },
    {
        "name": "VMR Utilities",
        "description": "Refresh Virus Metadata Resource (VMR) from ICTV and construct first/second pass sets."
    },
    {
        "name": "Dev Utilities",
        "description": "Endpoints for use by developers."
    }
]

app = FastAPI(
    title="GRAViTy V2",
    version="1.0",
    description=description,
    contact={
        "name": "Mayne Bio Analytics Ltd",
        "url": "http://www.mayneba.com",
        "email": "director@mayneba.com"
    },
    # RM < Check this license is appropriate
    license_info={
        "name": "Apache 2.0",
        "url": "https://www.apache.org/licenses/LICENSE-2.0.html",
    },
    openapi_tags=tags_metadata
)


'''Dev Endpoints'''


@app.get("/", tags=["Dev Utilities"])
async def read_root():
    return {"response": "API is healthy. Append the current URL to include '/docs/' at the end to visit the GUI."}


'''PL1 Entrypoints'''


@app.post("/pipeline_i_full/", tags=["Pipeline I"])
async def pipeline_i_full(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(run_pipeline_i_full, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_full(payload):
    pl = Pipeline_I(payload)
    pl.read_genome_desc_table()


@app.post("/pipeline_i_from_pphmmdb_construction/", tags=["Pipeline I"])
async def pipeline_i__from_pphmmdb_construction(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(
        run_pipeline_i_from_pphmmdb_construction, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_pphmmdb_construction(payload):
    pl = Pipeline_I(payload)
    pl.pphmmdb_construction()


@app.post("/pipeline_i_from_ref_virus_annotator/", tags=["Pipeline I"])
async def pipeline_i__from_ref_virus_annotator(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(run_pipeline_i_from_ref_virus_annotator, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_ref_virus_annotator(payload):
    pl = Pipeline_I(payload)
    pl.ref_virus_annotator()


@app.post("/pipeline_i_from_graph_generator/", tags=["Pipeline I"])
async def pipeline_i__from_graph_generator(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(run_pipeline_i_from_graph_generator, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_graph_generator(payload):
    pl = Pipeline_I(payload)
    pl.make_graphs()


@app.post("/pipeline_i_from_mutual_info_calculator/", tags=["Pipeline I"])
async def pipeline_i_from_mutual_info_calculator(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(
        run_pipeline_i_from_mutual_info_calculator, payload)
    return "Task fired successfully, running in background"


def run_pipeline_i_from_mutual_info_calculator(payload):
    pl = Pipeline_I(payload)
    pl.mutual_info_calculator()


'''PL2 entrypoints'''


@app.post("/pipeline_ii_full/", tags=["Pipeline II"])
async def pipeline_ii_full(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(run_pipeline_ii_full, payload)
    return "Task fired successfully, running in background"


def run_pipeline_ii_full(payload):
    pl = Pipeline_II(payload)
    pl.read_genome_desc_table()


@app.post("/pipeline_ii_from_pphmmdb_construction/", tags=["Pipeline II"])
async def pipeline_ii_from_pphmmdb_construction(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(
        run_pipeline_ii_from_pphmmdb_construction, payload)
    return "Task fired successfully, running in background"


def run_pipeline_ii_from_pphmmdb_construction(payload):
    pl = Pipeline_II(payload)
    pl.pphmmdb_construction()


@app.post("/pipeline_ii_from_ucf_virus_annotator/", tags=["Pipeline II"])
async def pipeline_ii_from_ucf_virus_annotator(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(
        run_pipeline_ii_from_ucf_virus_annotator, payload)
    return "Task fired successfully, running in background"


def run_pipeline_ii_from_ucf_virus_annotator(payload):
    pl = Pipeline_II(payload)
    pl.ucf_virus_annotator()


@app.post("/pipeline_ii_from_virus_classification/", tags=["Pipeline II"])
async def pipeline_ii_from_virus_classification(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(
        run_pipeline_ii_from_virus_classification, payload)
    return "Task fired successfully, running in background"


def run_pipeline_ii_from_virus_classification(payload):
    pl = Pipeline_II(payload)
    pl.virus_classification()


'''Utility entrypoints'''


@app.post("/scrape_vmr/", tags=["VMR Utilities"])
async def run_vmr_scrape(trigger: ScrapeData):
    payload = jsonable_encoder(trigger)
    status = scrape(payload)
    return status


@app.post("/construct_first_pass_vmr/", tags=["VMR Utilities"])
async def vmr_first_pass(trigger: FirstPass):
    payload = jsonable_encoder(trigger)
    status = first_pass(payload)
    return status


@app.post("/construct_second_pass_vmr/", tags=["VMR Utilities"])
async def vmr_second_pass(trigger: SecondPass):
    payload = jsonable_encoder(trigger)
    status = second_pass(payload)
    return status
