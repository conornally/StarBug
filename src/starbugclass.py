import os, sys
sys.stdout.write('\x1b[s..loading..')
sys.stdout.flush()
import logging, argparse, readline, glob


from fitsclass import FITS
from catclass import CATALOG
from sourceclass import Source
from common_tasks import *
import parse_config as parse_config

logging.basicConfig( level='DEBUG', format="%(message)s")
class StarBug:
    def __init__(self):
        self.fitslist={'Flat':[], 'Dark':[]}
        self.cataloglist={}
        self.commands ={'pre_adjust': self.pre_adjust,
                        'build_flat': self.build_flat,
                        'build_dark': self.build_dark,
                        'subtract_dark': self.darkframe_subtract,
                        'flatfielding': self.flatfield_divide,
                        'align': self.calculate_offset,
                        'stack': self.stack,
                        # analysis
                        'stats': self.stats,
                        'find': self.find,
                        # file manipulations
                        'update_header': self.update_header,
                        'dtype': self.convert_dtype,
                        'scale': self.scale,
                        'add': self.add,
                        'cut_below':self.cut_below,
                        # io
                        'load': self.file_loadin,
                        'loadcat': self.loadcat,
                        'show': self.display_loaded,
                        'header': self.display_header,
                        'd': self.delete,
                        'move': self.move,
                        'save': self.save,
                        'delete': self.delete_group,
                        'clean': self.clean,
                        #debug file i/o
                        'exportoffset':self.exportoffset,
                        'display':self.display,
                        # utils
                        'terminal': self.terminal,
                        'options': self.display_options,
                        'help': self.manual,
                        'exit': self.exit}
        self.reset_completer()
        if os.path.exists('config'): self.options = parse_config.load()
        else: logging.warning("\x1b[u\x1b[1;31mNo config file found in directory\x1b[0m")
        self.mainloop()


    ##########################
    # Data Reduction Regimes #
    ##########################
    def pre_adjust(self):
        pre_adjustments( self.get_group())#fitsfromtxt(self.readin('File.txt >> ')))
    def build_flat(self):
        flat = flatFrame_build( self.get_group(), self.fitslist['Dark'][0])   #fitsfromtxt(self.readin('Flatsfiles.txt >> ')), FITS( self.readin('Dark Frame.fits >> '))) )
        self.fitslist['Flat'] = [flat]
        self.display_loaded()
    def build_dark(self):
        dark= darkFrame_build( self.get_group())  # fitsfromtxt( self.readin('DarksFiles.txt >> '))))
        self.fitslist['Dark']= [dark]
        self.display_loaded()

    def darkframe_subtract(self):
        if self.fitslist['Dark']:
            darkFrame_subtract( self.get_group(), self.fitslist['Dark'])
        else: print('No Dark frame loaded')

    def flatfield_divide(self):
        if self.fitslist['Flat']:
            flatField_divide( self.get_group(), self.fitslist['Flat'])
        else: print('No flat field loaded')

    def calculate_offset(self):
        align( self.get_group() )

    def stack(self):
        fitslist = self.get_group()
        offsetfile = self.readin("offset file (auto)>> ")
        if offsetfile: readinOFFSET(fitslist, offsetfile)
        if len(fitslist) >=2:
            fitslist[0].stack( fitslist[1:], crop=True)
            #fitslist[0].export('tmp.fits', overwrite=True)
        

    
    ##################
    # Image Analysis #
    ##################

    def stats(self):
        basic_stats( self.get_group(), float(self.readin('Sigma Clip >> ')), int(self.readin('Iterations >> ')))

    def find(self):
        for f in self.get_group():
            if hasattr(self, 'options'):
                f.load_options(self.options)
            f.find()



    ######################
    # File Manipulations #
    ######################

    def convert_dtype(self):
        """Converts data type of fits pixel arrays
        """
        dtype_convert( self.get_group(), self.readin('New Data Type >> '))

    def update_header(self):
        """
        """
        group = self.get_group()
        key = self.readin('Header Key >> ')
        if key in group[0].header:
            print('%s: %s'%(key, group[0].header[key]))
            value = self.readin('Value >> ')
            for f in group:
                f.header[key] = value
                logging.debug('%s'%(f))
        else: print('No key %s'%key)

    def scale(self):
        a = float(self.readin("Scale Factor>> "))
        for f in self.get_group("group to scale> "):
            f.scale(a)

    def add(self):
        group = self.get_group()
        a = float(self.readin("+ "))
        for f in group:
            np.add(f.data, a, out=f.data)
            logging.debug(f)

    def cut_below(self):
        group = self.get_group()
        v = self.readin("absolute value >> ")
        for f in group:
            f.clip_below(v)
            logging.info(f)




    #############################
    # General Menu / Navigation #
    #############################
    def manual(self):
                    
        print('\nHelp Page:\n----------')
        print('\x1b[1;37mTEMPLATE\x1b[0m: Here is some words')

        print('\x1b[4;37m\nImage Reduction\x1b[0m')
        print('\x1b[1;37mpre_adjust\x1b[0m: Adjustments to fits file before data reduction')
        print('\x1b[1;37mbuild_flat\x1b[0m: Build flat field frame from list of raw flats and a dark frame')
        print('\x1b[1;37mbuild_dark\x1b[0m: Build dark frame from list of dark fields')
        print('\x1b[1;37msubtract_dark\x1b[0m: Subtracts Dark frame (group) from fits group')
        print('\x1b[1;37mflatfield\x1b[0m: Divides fits group by loaded Flat group')
        print('\x1b[1;37malign\x1b[0m: Get alignment offsets from reference image')
        print('\x1b[1;37mexportoffset\x1b[0m: Export offsets to file')
        print('\x1b[1;37mstack\x1b[0m: Stacks images based on alignment offset')

        print('\x1b[4;37m\nFILE/PIXEL Manipulations\x1b[0m')
        print('\x1b[1;37mdtype\x1b[0m: Change data type of pixel arrays')
        print('\x1b[1;37mheader\x1b[0m: Display header files')
        print('\x1b[1;37mupdate_header\x1b[0m: Change value in header files')
        print('\x1b[1;37madd\x1b[0m: Add value to pixels')
        print('\x1b[1;37mscale\x1b[0m: Scale pixels to power of scaling factor')

        print('\x1b[4;37m\nANALYSIS\x1b[0m')
        print('\x1b[1;37mstats\x1b[0m: Get sigma clipped stats of arrays')

        print('\x1b[4;37m\nBASIC I/O\x1b[0m')
        print('\x1b[1;37mload\x1b[0m: Give list, or pathfile of fits fits to be loaded into program')
        print('\x1b[1;37mloadcat\x1b[0m: Load catalog, and give catalog style and associated fitsfile')
        print('\x1b[1;37mshow\x1b[0m: Display the currently loaded files')
        print('\x1b[1;37mdelte\x1b[0m: Delete group from loaded files')
        print('\x1b[1;37msave\x1b[0m: Saves all files from selected load group')
        print('\x1b[1;37mmove\x1b[0m: Move fits between loaded groups')
        print('\x1b[1;37md\x1b[0m: Delete fits from loaded groups')
        print('\x1b[1;37mclean\x1b[0m: Deletes out/ directory')
        print('\x1b[1;37mterminal\x1b[0m: Enter terminal mode')
        print('\x1b[1;37mdisplay\x1b[0m: Open fits displaying program [specified in config file]')
        print('\x1b[1;37mhelp\x1b[0m: Prints this page')
        print('\x1b[1;37mexit\x1b[0m: Leaves StarBug\n')

    def file_loadin(self):
        self.complete_style = "PATH"
        inString = self.readin('Load Files [Single \x1b[1;37mOR\x1b[0m Space Separated \x1b[1;37mOR\x1b[0m *.fits globbed files \x1b[1;37mOR\x1b[0m From File or Paths]\n>> ') 
        try:
            inlist=[]
            for FILES in inString.split(' '):
                for filename in sorted(glob.glob(FILES)):
                    inlist.append(FITS(filename))
        #try: inlist = [FITS(filename) for filename in sorted(glob.glob(inString.split(' ')))] #single or space separated
        except: 
            inlist = fitsfromtxt( inString ) #path file

        self.fitslist[ self.readin('Name for group of files >> ') ] = inlist
        self.reset_completer()

    def loadcat(self):
        self.complete_style = "PATH"
        #self.extension = '.cat'
        infile = self.readin('Load in Catalog >> ')
        catalog_style = self.readin('Catalog style [ sextractor , custom ] >> ')
        if infile: 
            name = self.readin('Catalog Name >> ')
            #self.extension = '.fits'
            fitsfilename = self.readin('Assosiated fitsfile (if applicable) >> ')
            self.cataloglist[name] = CATALOG(fitsfile = fitsfilename, catalog_style=catalog_style, catalog_filename=infile)
        self.reset_completer()



    def display_loaded(self):
        print('\nLoaded FITS Files')
        for item in self.fitslist:
            print(item +str(':') + str(self.fitslist[item] ) )
        print('\nLoaded Catalogs')
        for item in self.cataloglist:
            print("%s: %s"%(item, self.cataloglist[item]))

    def display_header(self):
        for f in self.get_group():
            print(f)
            print(repr(f.header))


    def get_group(self, string='Name of loaded group >> '):
        self.complete_list = self.fitslist.keys()
        instring = self.readin(string)
        return_list = []
        if instring in self.fitslist.keys():
            return_list= self.fitslist[instring]
        else: 
            print("Group %s doesn't exist"%instring)
            if self.readin('Create %s y/n: '%instring)=='y': 
                self.fitslist[instring]=[]
                return_list= self.fitslist[instring]
        self.reset_completer()
        return return_list

    def delete_group(self):
        #Deletes one of the loaded groups
        group = self.readin('Delete Group Name >> ')
        if group in self.fitslist.keys():
            del self.fitslist[group]
            logging.info('%s Deleted'%group)
        else: logging.info('No group named %s'%group)


    def yank(self):
        pass

    def move(self):
        group = self.get_group()
        selected = self.single_select(group)
        new_group = self.get_group('Move to >> ')
        print(new_group)
        i = int(self.readin('Insert to index >> '))
        if i <= len(new_group) and i >=0:
            new_group.insert(i, selected)
            group.pop(group.index(selected))
        else: print('Cannae dae that')


    def delete(self):
        group = self.get_group() 
        selected = self.single_select(group)
        if self.readin('Delete %s y/n: '%selected) == 'y': group.pop( group.index(selected))
        print(group)

    def single_select(self, group):
        for i,f in enumerate(group): print('%d: %s'%(i,f))
        i = int(self.readin('Select Index >> '))
        try: return(group[i])
        except: return None




    def save(self, filename=''):
        if not os.path.isdir('out'): os.system('mkdir out')
        overwrite = False
        for f in self.get_group():
            exists = os.path.exists('out/'+f.name)
            if exists and not overwrite:
                o = self.readin('Overwrite %s [Y/n/(a)]: '%f)
                if o=='y' or o=='Y': f.export('out/'+f.name, overwrite=True)
                elif o=='a' or o=='A':
                    overwrite = True
                    f.export('out/'+f.name, overwrite=True)
            elif not exists: f.export('out/'+f.name, overwrite=True)

    def saveas(self):
        filename = self.readin('Filename / filename base >> ')
        
        



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

    def display_options(self):
        """
        Prints options from config file
        If no config file was found, then iit gives option to find one
        Options can be altered here to
        """
        continue_flag = True
        if hasattr(self, 'options'): parse_config.display(self.options)
        else:
            path = self.readin('No config file found, enter path >> ')
            if os.path.exists(path): 
                self.options = parse_config.load(path)
                parse_config.display(self.options)
            else: 
                print('Cant open file')
                continue_flag = False
        while continue_flag:
            in_string = self.readin('Option = value >> ')
            if not in_string: continue_flag = False
            else: 
                key = parse_config.parse_key(in_string)
                value = parse_config.parse_value(in_string)
                if key in self.options.keys(): self.options[key] = value
        if hasattr(self, 'options'):
            print('')
            parse_config.display(self.options)
            for flist in self.fitslist.values():
                for f in flist: 
                    if f: f.load_options(self.options)


    ###############
    #DEBUG FILE IO#
    ###############

    def exportoffset(self):
        """FUNC: exports offsets in fitslist
        """
        exportOFFSET(self.get_group())

    def display(self):
        """FUNC: displays loaded group
        """
        display = self.options['DISPLAY']
        os.system("%s -console &"%(display))


    def exit(self):
        quit('Bye!')

    def readin(self, string=''):
        try: return(raw_input(string))
        except: return(input(string))

    def complete(self, text, state):
        if self.complete_style == "PATH":
            self.complete_list = glob.glob(text+"*")
        for cmd in self.complete_list:
            if cmd.startswith(text):
                if not state: return cmd
                else: state -= 1

    def reset_completer(self):
        self.complete_list = self.commands.keys()
        self.complete_style = "STANDARD"
        self.extension = '.fits'

    def _completeGROUP(self, text, state):
        for group in self.fitslist.keys():
            if group.startswith(text):
                if not state: return group
                else: state =-1

    def mainloop(self):
        print('\x1b[u\x1b[1;36mHello! Welcome to \x1b[1;32mStarBug\x1b[0m\n')
        command=self.manual
        readline.set_completer_delims('')
        #readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(self.complete)
        while command != exit:
            cmd_in = self.readin('> ')
            if cmd_in in self.commands:
                self.commands[cmd_in]()



StarBug()

