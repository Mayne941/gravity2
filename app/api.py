from app.pipeline_i import Pipeline_I
from app.pipeline_ii import Pipeline_II
from app.utils.api_classes import Pipeline_i_data, Pipeline_ii_data

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
]

app = FastAPI(
    title = "GRAViTy V2",
    version = "1.0",
    description = description,
    contact = {
        "name": "Mayne Bio Analytics Ltd",
        "url": "http://www.mayneba.com",
        "email": "director@mayneba.com"
    },
    # RM < Check this license is appropriate
    license_info={
        "name": "Apache 2.0",
        "url": "https://www.apache.org/licenses/LICENSE-2.0.html",
    },
    openapi_tags = tags_metadata
)


@app.get("/")
async def read_root():
    return {"response": "healthy"}

@app.post("/pipeline_i/", tags=["Pipeline I"])
async def pipeline_i(payload: Pipeline_i_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(run_pipeline_i, payload)
    return "Task fired successfully, running in background"

def run_pipeline_i(payload):
    pl = Pipeline_I(payload)
    pl.main()

@app.post("/pipeline_ii/", tags=["Pipeline II"])
async def pipeline_ii(payload: Pipeline_ii_data, background_tasks: BackgroundTasks):
    payload = jsonable_encoder(payload)
    background_tasks.add_task(run_pipeline_ii, payload)
    return "Task fired successfully, running in background"

def run_pipeline_ii(payload):
    pl = Pipeline_II(payload)
    pl.main()
