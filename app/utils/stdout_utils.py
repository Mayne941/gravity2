import sys 

def progress_bar(out):
    sys.stdout.write(out)
    sys.stdout.flush()

def clean_stdout():
    sys.stdout.write("\033[K")
    sys.stdout.flush()

def error_handler(out, err, name):
    if err != b"":
        raise SystemExit(f"Something went wrong with {name}: {out, err}")
