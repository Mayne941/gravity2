import json
import requests
import os
import time

def batch(indir, host):
    errors = []
    for fname in os.listdir(indir):
        st = time.time()
        with open(f"{indir}/{fname}") as f:
            payload = json.load(f)

            response = requests.post(f"{host}end_to_end_programmatic/", json=payload)
            if response:
                continue
            else:
                errors.append(fname)
        print(f"TTC {fname}: {st - time.time()}")

    print(f"FINISHED PROCESSING {len(os.listdir(indir))} SAMPLES\n{len(errors)} ERRORS:\n{errors}")

if __name__ == "__main__":
    run_files_dir = "data/2023_tp_runfiles"
    host = "http://127.0.0.1:8001/"
    batch(run_files_dir, host)
