import json

def main(pl_type):
    with open(f"{pl_type}.json", "r") as f:
        payload = json.loads(f.read())
    return payload
