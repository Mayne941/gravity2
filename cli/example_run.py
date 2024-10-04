import json
from app.api import process_json, run_pipeline_i_full

if __name__ == "__main__":
    with open(f"data/eval/example_runfile.json") as f:
        payload = process_json(json.load(f))
        run_pipeline_i_full(payload)
    print(f"GRAViTy-V2 test complete. Check ./output/GRAViTyV2_example_run/ for output! Full output descriptions may be found on the GitHub readme.")
