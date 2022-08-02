from app.pipeline_i import Pipeline_I
from app.pipeline_ii import Pipeline_II
from app.utils.api_classes import Pipeline_i_data, Pipeline_ii_data

from fastapi import FastAPI

app = FastAPI()


@app.get("/")
def read_root():
    return {"response": "healthy"}

@app.post("/pipeline_i/")
def pipeline_i(payload: Pipeline_i_data):
    #pl = Pipeline_I(payload)
    #pl.main(payload)
    # Background task
    return []

@app.post("/pipeline_ii/")
def pipeline_ii(payload: Pipeline_ii_data):
    #pl = Pipeline_I(payload)
    #pl.main(payload)
    # Background task
    return "Running in background"