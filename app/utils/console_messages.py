from termcolor import colored

def section_header(msg):
    '''Print a section header to console and output log'''
    print(colored("#"*90, 'green'))
    print(colored(msg, 'green'))
    print(colored("#"*90, 'green'))
