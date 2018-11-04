import os, sys
sys.stdout.write('\x1b[s..loading..')
sys.stdout.flush()
import logging, argparse, readline

from src.fitsclass import FITS
from src.alsclass import ALS_DATA
from src.sourceclass import Source
from src.common_tasks import *

logging.basicConfig( level='DEBUG', format="%(message)s")
class StarBug:
    def __init__(self):
        self.fitslist={'Flat':[], 'Dark':[]}
        self.alslist=[]
        self.sourcelist=[]
        self.commands ={'pre_adjust': self.pre_adjust,
                        'build_flat': self.build_flat,
                        'build_dark': self.build_dark,
                        'load': self.file_loadin,
                        'show': self.display_loaded,
                        'save': self.save,
                        'clean': self.clean,
                        'terminal': self.terminal,
                        'help': self.manual,
                        'exit': self.exit}
        self.mainloop()


    ##########################
    # Data Reduction Regimes #
    ##########################
    def pre_adjust(self):
        pre_adjustments( self.fitslist[ self.readin('Loaded Group Name >> ')])#fitsfromtxt(self.readin('File.txt >> ')))
    def build_flat(self):
        flat = flatFrame_build( self.fitslist[ self.readin('Loaded Group Name >> ')], self.fitslist['Dark'][0])   #fitsfromtxt(self.readin('Flatsfiles.txt >> ')), FITS( self.readin('Dark Frame.fits >> '))) )
        self.fitslist['Flat'] = [flat]
        self.display_loaded()
    def build_dark(self):
        dark= darkFrame_build( self.fitslist[ self.readin('Loaded Group Name >> ')])  # fitsfromtxt( self.readin('DarksFiles.txt >> '))))
        self.fitslist['Dark']= [dark]
        self.display_loaded()

    def darkframe_subtract(self):
        darkframe_subtract( self.fitslist[ self.readin('Loaded Group Name >> ')]
        


    #############################
    # General Menu / Navigation #
    #############################
    def manual(self):
        print('\nHelp Page:\n----------')
        print('\x1b[1;37mpre_adjust\x1b[0m: adjustments to fits file before data reduction')
        print('\x1b[1;37mbuild_flat\x1b[0m: Build flat field frame from list of raw flats and a dark frame')
        print('\x1b[1;37mbuild_dark\x1b[0m: Build dark frame from list of dark fields')
        print('\x1b[1;37mload\x1b[0m: Give list, or pathfile of fits fits to be loaded into program')
        print('\x1b[1;37mshow\x1b[0m: Display the currently loaded files')
        print('\x1b[1;37msave\x1b[0m: Saves all files from selected load group')
        print('\x1b[1;37mclean\x1b[0m: Deletes out/ directory')
        print('\x1b[1;37mterminal\x1b[0m: Enter terminal mode')
        print('\x1b[1;37mhelp\x1b[0m: this page')
        print('\x1b[1;37mexit\x1b[0m: leaves StarBug\n')

    def file_loadin(self):
        inString = self.readin('Load Files [Single OR Comma Separated OR From File or Paths]\n>> ') 
        try: inlist = [FITS(filename) for filename in inString.split(',')]
        except: inlist = fitsfromtxt( inString )
        self.fitslist[ self.readin('Name for group of files >> ') ] = inlist


    def display_loaded(self):
        print('\nLoaded Files')
        for item in self.fitslist:
            print(item +str(':') + str(self.fitslist[item] ) )

    def save(self):
        try: os.system('mkdir out')
        except: pass
        group = self.readin('Group Name >> ')
        for f in self.fitslist[group]:
            f.export('out/'+f.filename)

    def clean(self):
        if self.readin('Delete out/ folder [y/n]') =='y':
            os.system('rm -r out/')

    def terminal(self):
        print('Entering Terminal Mode. To exit type "EXIT"')
        while True:
            cmd=self.readin('$~ ')
            if cmd=='EXIT': break
            else: os.system(cmd)
        print('Terminal Mode Exitted')



    def exit(self):
        quit('Bye!')

    def readin(self, string=''):
        try: return(raw_input(string))
        except: return(input(string))

    def complete(self, text, state):
        for cmd in self.commands.keys():
            if cmd.startswith(text):
                if not state:
                    return cmd
                else:
                    state -= 1

    def mainloop(self):
        print('\x1b[u\x1b[1;36mHello! Welcome to \x1b[1;32mStarBug\x1b[0m\n')
        command=self.manual
        while command != exit:
            readline.parse_and_bind("tab: complete")
            readline.set_completer(self.complete)
            cmd_in = self.readin('> ')#raw_input('> ')
            if cmd_in in self.commands:
                self.commands[cmd_in]()
            #elif cmd_in == '\r': print('r')
            #elif cmd_in == '\n': print('n')
            #else: print('Command not recognised: type "help" for manual')



StarBug()

