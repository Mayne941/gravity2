from termcolor import colored

def raise_gravity_error(msg):
    raise SystemExit(f"{colored('FATAL ERROR', 'red')}: {msg}")
