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
    if err != b"":
        raise SystemExit(f"Something went wrong with {name}: {out, err}")
