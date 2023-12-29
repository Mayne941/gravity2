import subprocess as sp
from termcolor import colored

def shell(args, calling_fn="Misc shell function", ret_output=False):
    '''Call Bash shell with input string as argument'''
    whitelist = ["mafft"]
    _ = sp.Popen(args, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = _.communicate()
    if not any(x.lower() in whitelist for x in whitelist):
        sp_error_handler(
            out, err, f"{colored('ERROR', 'red')}: Something went wrong with CLI call: {calling_fn}, dumping STDOUT/STDERR to shell.")
    if ret_output:
        return out

def sp_error_handler(out, err, name):
    '''Kill program if error found; used in combo with POpen commands.'''
    if err != b"":
        raise SystemExit(f"***\n{name}:\nOut: {out}\nErr: {err}\n***")
