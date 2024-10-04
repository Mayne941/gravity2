import pandas as pd

def join_input(payload):
    '''Takes two input vmr-like documents, joins and saves as csv'''
    df = pd.concat([pd.read_csv(payload["vmr_1"]), pd.read_csv(payload["vmr_2"])])
    crap_pandas = [i for i in df.columns if "Unnamed" in i]
    for i in crap_pandas:
        df = df.drop(columns=i)
    df.to_csv(payload["vmr_joined"], index=False)
    finished = f"CSV join completed successfully. Output saved to {payload['vmr_joined']}"
    print(finished)
    return finished
