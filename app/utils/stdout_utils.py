import sys


def progress_bar(out):
    '''Print progress bar to terminal.'''
    sys.stdout.write(out)
    sys.stdout.flush()


def clean_stdout():
    '''Flush stdout.'''
    sys.stdout.write("\033[K")
    sys.stdout.flush()


def error_handler(out, err, name):
    '''Kill program if error found; used in combo with POpen commands.'''
    if "Segmentation fault" in str(out):
        raise SystemExit(f"ERROR: Muscle ran out of memory :( Details follow:\n{str(out)[-100:]}")
    if err != b"":
        raise SystemExit(f"Something went wrong with {name}:\n {out, err}")
