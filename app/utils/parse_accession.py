import re

def get_accession_regex():
    return re.compile(r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{2}_[0-9]{6}|SRR[0-9]{7,8}|[a-zA-Z]{4,6}[0-9]{6,9}")
