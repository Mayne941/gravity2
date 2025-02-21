import pytest
from app.utils.load_premade_pipeline import main as load_premade

def test_load_premade():
    pipelines = [
        "dev/premade_pl_files/similar_viruses",
        "dev/premade_pl_files/divergent_viruses",
        "dev/premade_pl_files/long_single_orf_viruses"
    ]

    for pipeline in pipelines:
        result = load_premade(pipeline)
        assert isinstance(result, dict)
