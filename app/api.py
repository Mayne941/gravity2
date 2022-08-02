from app.pipeline_i import Pipeline_I
from app.pipeline_ii import Pipeline_II
from app.utils.api_classes import Pipeline_i_data, Pipeline_ii_data

from fastapi import FastAPI

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
        "name": "Richard Mayne",
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
def read_root():
    return {"response": "healthy"}

@app.post("/pipeline_i/", tags=["Pipeline I"])
def pipeline_i(payload: Pipeline_i_data):
    #pl = Pipeline_I(payload)
    #pl.main(payload)
    # Background task
    return []

@app.post("/pipeline_ii/", tags=["Pipeline II"])
def pipeline_ii(payload: Pipeline_ii_data):
    #pl = Pipeline_I(payload)
    #pl.main(payload)
    # Background task
    return "Running in background"