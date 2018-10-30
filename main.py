import os, sys
sys.stdout.write('\x1b[s..loading..')
sys.stdout.flush()
import logging, argparse, readline

from src.fitsclass import FITS
from src.alsclass import ALS_DATA
from src.sourceclass import Source
from src.common_tasks import *

logging.basicConfig( level='DEBUG', format="%(message)s")

def exit():
    quit('Bye!')

def manual():
    print('Help Page:')
    print('help: this page')
    print('pre_adjust: adjustments to fits file before data reduction')
    print('exit: leaves StarBug')

def pre_adjust():
    pre_adjustments( fitsfromtxt(raw_input('File.txt >> ')))

commands = {'pre_adjust': pre_adjust,
            'help':manual,
            'exit':exit}

def readin(string=''):
    try: return(raw_input(string))
    except: return(input(string))

def complete(text, state):
    for cmd in commands.keys():
        if cmd.startswith(text):
            if not state:
                return cmd
            else:
                state -= 1

def mainloop():
    print('\x1b[u\x1b[1;36mHello! Welcome to \x1b[1;32mStarBug\x1b[0m\n')
    command=manual
    while command != exit:
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete)
        cmd_in = readin('> ')#raw_input('> ')
        if cmd_in in commands:
            commands[cmd_in]()
        else: print('Command not recognised: type "help" for manual')




mainloop()

