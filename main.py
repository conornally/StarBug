import os, sys, logging, argparse

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

def pre_adjust():
    pre_adjustments( fitsfromtxt(raw_input('File.txt >> ')))

commands = {'pre_adjust': pre_adjust,
            'help':manual,
            'exit':exit}



def mainloop():
    print('\x1b[1;36mHello! Welcome to \x1b[1;32mStarBug\x1b[0m')
    command=manual
    while command != exit:
        cmd_in = raw_input('> ')
        if cmd_in in commands:
            commands[cmd_in]()
        else: print('Command not recognised: type "help" for manual')




mainloop()

