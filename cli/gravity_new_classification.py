import json
import argparse
from app.api import process_json, run_pipeline_i_full

def parse():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog='Welcome to GRAViTy-V2 CLI')
    parser.add_argument('--runfile', required=True,
                        default="", help="File path for your GRAViTy-V2 run parameters file (json format).")
    return parser

def fire(parser):
    with open(parser.runfile) as f:
        payload = process_json(json.load(f))
        run_pipeline_i_full(payload)
    print(f"GRAViTy-V2 run complete.")


if __name__ == "__main__":
    fire(parse().parse_args())
