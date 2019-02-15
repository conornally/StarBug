import os, numpy as np, logging, time
np.warnings.filterwarnings('ignore')
from scipy import optimize
import matplotlib.pyplot as plt

from astropy.coordinates import match_coordinates_sky, Angle, SkyCoord
import astropy.units as u
from astropy import wcs

try:from sourceclass import Source
except:from src.sourceclass import Source

try: from fitsclass import FITS
except: from src.fitsclass import FITS
try:from parse_config import *
except:from src.parse_config import *

class CATALOG(object):
    def __init__(self, fitsfile=None, catalog_style='custom', catalog_filename='', catalog=None, configfile='config'):
        """
        INPUT: fitsfile: filename of fitsfile OR FITSCLASS object
               catalog_style: 'sextractor' 
                              'DAOPHOT find'
                              'DAOPHOT allstar'
                              'DAOStarFinder'
                              'photutils DAOPHOT'
                              'custom'
               catalog_filename: filename of catalog
               catalog: catalog data, if loaded from any photutils regimes?
        """
        if type(fitsfile)==None or fitsfile=='': logging.debug("No fitsfile loaded")
        elif type(fitsfile)==str: self.fits = FITS(fitsfile)
        elif type(fitsfile)==FITS: self.fits = fitsfile

        self.configfile = configfile
        if hasattr(self, 'fits'): self.loadconfig() #FOR NOW!! need to allow not including fits properlly

        if(catalog_style=="starbug"):
            self._loadStarBugData(catalog_filename)
        
        elif(catalog_style=='pleiades'):
            self.edgecase_constructPleiades()

        else:
            self.raw_data = self.construct_raw_data(catalog_style, catalog_filename, catalog)
            self.sourcelist = np.empty( (len(self.raw_data)), dtype=object)
            self.build_sourcelist()

        name = catalog_filename
        while '/' in name: name = name[name.index('/')+1:]
        self.name = name

    def edgecase_constructPleiades(self, catalog_filename=''):
        """
        THis is temporary
        Constructs pleides mainsequence data catalog, with colours not mags
        should i have colours in source?
        """
        data = np.genfromtxt("catalogs/pleiades_sloane")[:10]
        print(data.shape)
        self.sourcelist = np.empty(len(data), dtype=object)
        for i, line in enumerate(data):
            self.sourcelist[i] = Source(mag=line[2:], bands=2)
            self.sourcelist[i].set_colours(line[0], line[1])
        print(self.sourcelist)

    def construct_raw_data(self, catalog_style='custom', catalog_filename='', catalog=None):
        """INPUT: same as __init__
            FUNC: creates raw_data numpy array in correct format
        """
        skip_header=0
        skip_footer=0
        delimiter=None
        comments='#'

            
        if(catalog_style=='sextractor'):
            ID = 0
            flux = 1
            fluxerr=2
            mag = 3
            magerr=4
            ra = 5
            dec =6
            x = 7
            y = 8
            skip_header = 9
            #skip_footer = 0
        
        else: #custom
            ID =    get_value('IDcol',      configfile=self.configfile, dtype=int)
            flux =  get_value('Fluxcol',    configfile=self.configfile, dtype=int)
            fluxerr=get_value('FluxErrcol', configfile=self.configfile, dtype=int)
            mag =   get_value('Magcol',     configfile=self.configfile, dtype=int)
            magerr =get_value('MagErrcol',  configfile=self.configfile, dtype=int)
            ra =    get_value('RAcol',      configfile=self.configfile, dtype=int)
            dec =   get_value('DECcol',     configfile=self.configfile, dtype=int)
            x =     get_value('Xcol',       configfile=self.configfile, dtype=int)
            y =     get_value('Ycol',       configfile=self.configfile, dtype=int)
            skip_header = get_value('skip_header', configfile=self.configfile,dtype=int)
            skip_footer = get_value('skip_footer', configfile=self.configfile,dtype=int)
            #comments    = get_value('Comments', configfile=self.configfile)
            #delimiter =   get_value('delimiter', configfile=self.configfile)

        load_in_data = np.genfromtxt( catalog_filename, skip_header=skip_header, skip_footer=skip_footer, comments=comments, delimiter=delimiter)
        raw_data = np.zeros((np.shape(load_in_data)[0],9))

        raw_data[:,0] = load_in_data[:,ID]
        raw_data[:,1] = load_in_data[:,x]
        raw_data[:,2] = load_in_data[:,y]
        raw_data[:,3] = load_in_data[:,ra]
        raw_data[:,4] = load_in_data[:,dec]
        raw_data[:,5] = load_in_data[:,flux]
        raw_data[:,6] = load_in_data[:,fluxerr]
        raw_data[:,7] = load_in_data[:,mag]
        raw_data[:,8] = load_in_data[:,magerr]
        return raw_data

    def _loadStarBugData(self, filename):
        """TEMPORARY FIX
            loads in starbug *.sb files, but only with export(style=mean tiles) for now
        """
        with open(filename,'r') as catalog:
            name= catalog.readline().split(' ')
            length= int(catalog.readline().split(' ')[-1])
            numEpochs = int(catalog.readline().split(' ')[-1])
            numBands = int(catalog.readline().split(' ')[-1])
            self.sourcelist = np.empty((length), dtype=object)
            flux = [ 4*i + 3 for i in range(numBands)]
            fluxerr = [ 4*i + 4 for i in range(numBands)]
            mag = [ 4*i + 5 for i in range(numBands)]
            magerr = [ 4*i + 6 for i in range(numBands)]
        data = np.genfromtxt(filename, skip_header=5)
        for i, line in enumerate(data):
            self.sourcelist[i] = Source(ID=i, ra=line[1], dec=line[2],
                                        flux=line[flux], fluxerr=line[fluxerr],
                                        mag=line[mag], magerr=line[magerr],
                                        epochs=numEpochs, bands=numBands)

            

    def setup(self):
        """
        __init__ second round, without filenames
        """
        self.cat = np.zeros((len(self.sourcelist), 2))
        for i,s in enumerate(self.sourcelist):
            s._voidCalcTileMeans()
            self.cat[i]=[s.RA, s.DEC]
        self.cat = SkyCoord( ra=Angle( self.cat[:,0], unit=u.deg), dec=Angle( self.cat[:,1], unit=u.deg) )

    def fromfile(self, filename, skipheader=3):
        try: return(np.genfromtxt(filename, skip_header=skipheader))
        except: return(np.array([]))

    def fromfits(self, filename):
        self.fits = FITS( filename )    

    def loadconfig(self):
        if hasattr(self, 'fits'):
            if 'CHNLNUM' in self.fits.header: i=self.fits.header['CHNLNUM']
            elif 'FILTER' in self.fits.header: i=self.fits.header['FILTER']
            else:i = str(input("Wavelength Band [int]: "))
        else:    i = str(input("Wavelength Band [int]: "))
        if not i: i='0'
        i = i.replace(' ','')

        config = {'BAND':i}
        config['ZP_Flux'] = float(get_value("Band%s Zero Point"%i, configfile=self.configfile))
        config['ZP_ERR'] = float(get_value("Band%s Zero Point Error"%i, configfile=self.configfile))
        config['Colour Correction'] = float(get_value("Band%s Colour Correction"%i, configfile=self.configfile))
        config['Max Mag'] = float(get_value("Band%s Max Mag"%i, configfile=self.configfile))
        config['Min Mag'] = float(get_value("Band%s Min Mag"%i, configfile=self.configfile))
        config['Mag Error'] = float(get_value("Band%s Mag Error"%i, configfile=self.configfile))
        config['Max Sharp'] = float(get_value("Band%s Max Sharp"%i, configfile=self.configfile))
        config['Min Sharp'] = float(get_value("Band%s Min Sharp"%i, configfile=self.configfile))
        config['Max Round'] = float(get_value("Band%s Max Round"%i, configfile=self.configfile))
        config['Min Round'] = float(get_value("Band%s Min Round"%i, configfile=self.configfile))

        self.config = config

    def xy2radec(self):
        WCS = wcs.WCS(self.fits.header, relax=False)
        self.raw_data[:,1:3] = WCS.all_pix2world(self.raw_data[:,1:3],0)

    def build_sourcelist(self):
        #self.xy2radec()
        for i, l in enumerate(self.raw_data):
            self.sourcelist[i] = Source( *l )


    def combine(self, CAT):
        """
        INPUT: list or single instance of ALS_DATA()
        FUNC: Appends the sourcelist onto bottom of self,sourcelist
        """
        if type(CAT) == CATALOG: CAT = [CAT]
        for cat in CAT: 
            self.sourcelist = np.append(self.sourcelist, cat.sourcelist)
            self.size = len(self.sourcelist)
            logging.info("COMBINING: {} += {} --> Sources: {}".format(self, cat, self.size))

    def calibrate(self):
        """
        TEMP
        """
        logging.info('This is not implemented')


    def xcalibrate(self, unitscale='None'):
        """
        >Converts DN/s to mJy //needs to be generalised
        >Applies colour correction
        >Applies zero point correction
        // currently only setup for spitzer
        """
        logging.info('CALIBRATING: %s'%self)

        """
        DN2MJy = self.fits.header['FLUXCONV'] / self.fits.header['EXPTIME']
        MJy2mJy = (1e9/4.254517e10)*(self.fits.header['PXSCAL1']**2.)

        factor=1. 
        if unitscale =='Dns-1_mJy':
            factor = DN2MJy * MJy2mJy 
        """

        ZP_Flux = self.config['ZP_Flux']
        ZP_ERR = self.config['ZP_ERR']
        colour_correction = self.config['Colour Correction']
        logging.debug("ZP: %f, dZP: %f, CC: %f"%(ZP_Flux, ZP_ERR, colour_correction))

        for s in self.sourcelist:
            s.flux = factor* 10.**((s.mag-25.)/-2.5) 
            s.fluxerr= (s.magerr * s.flux) / (1.0857*colour_correction)
            s.flux/=colour_correction
            s.mag = -2.5*np.log10( s.flux / ZP_Flux)
            s.magerr= 1.0857 * (s.fluxerr / s.flux)


    def data_crop(self):
        #need to get config stats at the very start...
        logging.info("CROPPING: %s"%self) 
        old = self.size
        for s in self.sourcelist:
            s._voidCalcTileMeans()
            if s.MAGERR > self.config['Mag Error']: s.quality = False
            elif s.MAG > self.config['Max Mag']: s.quality = False
            elif s.MAG < self.config['Min Mag']: s.quality = False
            elif s.SHARP > self.config['Max Sharp']: s.quality = False
            elif s.SHARP < self.config['Min Sharp']: s.quality = False
            #elif s.ROUND > self.config['Max Round']: s.quality = False
            elif s.ROUND < self.config['Min Round']: s.quality = False
          
        self.sourcelist = [s for s in self.sourcelist if s.quality]
        self.size = len(self.sourcelist)
        logging.info("Source: %d --> %d"%(old, self.size))


    def NullMatch(self):
        """FUNC: this is a temporary fix
                 if ou do epoch/tile match on catalogs with different numbers of bands, it will cash as the arrays are different sizes
                 this adds a completely NULL band onto the sourcelist
        """
        for i, source in enumerate(self.sourcelist):
            source.append_Band(bad=True)
            logging.debug(source)

    def BandMatch(self, CAT):
        logging.info("\x1b[1;33mMATCHING\x1b[0m Bands: %s <-- %s"%(self, CAT))
        for s in self.sourcelist: s._voidCalcTileMeans()
        self.match(CAT, regime='band')
        
    def EpochMatch(self, CAT):
        logging.info("\x1b[1;33mMATCHING\x1b[0m Epochs: %s <-- %s"%(self, CAT))
        for s in self.sourcelist: s._voidCalcTileMeans()
        self.match(CAT, regime='epoch')

    def TileMatch(self, CAT):
        logging.info("\x1b[1;33mMATCHING\x1b[0m Tiles: %s <-- %s"%(self, CAT))
        for s in self.sourcelist: s._voidCalcTileMeans()
        self.match(CAT, regime='tile')

    def match(self, CAT, regime='band'):
        """
        INPUT: >List or single instance of ALS_DATA(), 
               >regime is 'band' or 'epoch' matching 'tile'
               >> for now, i think all epoch matching should be done first
        FUNC: adds additional data onto the sources in l=source list
        (should i retake averages inbetween matches?)
        this will be a loooot slower though, but maybe more representative
        """
        if type(CAT) == CATALOG: CAT = [CAT]
        
        for cat in CAT:
            logging.info('...%s'%cat)
            radec = np.array( [ [s.RA, s.DEC ] for s in self.sourcelist])
            radec = SkyCoord( ra=Angle( radec[:,0], unit=u.deg), dec=Angle( radec[:,1], unit=u.deg) )
            for s in cat.sourcelist: s._voidCalcTileMeans()
            catradec = np.array( [ [s.RA, s.DEC] for s in cat.sourcelist] )
            catradec = SkyCoord( ra=Angle( catradec[:,0], unit=u.deg), dec=Angle( catradec[:,1], unit=u.deg) )
            
            idx, d2d, d3d = match_coordinates_sky(radec, catradec)
            bad = np.where( d2d.arcsec > 1.0)[0]
            for i, IDX in enumerate(idx):
                if i not in bad:
                    ## need to put more conditions here
                    if d2d[i] == min(d2d[np.where( idx==IDX)]):
                        if regime=='epoch': self.sourcelist[i].append_Epoch(cat.sourcelist[IDX])
                        elif regime=='band': self.sourcelist[i].append_Band(cat.sourcelist[IDX])
                        elif regime=='tile': self.sourcelist[i].append_Tile(cat.sourcelist[IDX])
                        self.sourcelist[i].quality= True
                if((i in bad) and True):
                    if regime=='epoch': self.sourcelist[i].append_Epoch(bad=True)
                    elif regime=='band': self.sourcelist[i].append_Band(bad=True)
                    elif regime=='tile':
                        self.sourcelist[i].append_Tile(bad=True)
                        #cat.sourcelist[IDX].append_Tile(bad=True)
                        #self.sourcelist = np.append(self.sourcelist, cat.sourcelist[IDX] )
            if regime=='tile':
                for i in range(len(cat.sourcelist)):
                    if i not in idx:
                        self.sourcelist = np.append(self.sourcelist, cat.sourcelist[i])
                #If i dont crop out the bad sources i will NEED to append 9999 values to sourcelist objects..
            #self.sourcelist = [s for s in self.sourcelist if s.quality]
            for s in self.sourcelist: s._voidCalcTileMeans()
            logging.info("Sources: %d"%len(self.sourcelist))
            logging.info("Not Matched: %d"%len(bad))

    def single_match(self, cat):
        return match_coordinates_sky(self.cat, cat)

    def display(self):
        print("Catalog: %s sources: %i"%(self, len(self.sourcelist)))
        for i in range(len(self.sourcelist)):
            if i<len(self.sourcelist):
                print(self.sourcelist[i])

    def exportRegionfile(self):
        #Func: creates ds9 region file in WCS degrees
        if not os.path.isdir('out'):
            os.system('mkdir out/')
        with open("out/%s.reg"%self.name[:-4],'w') as outreg:
            logging.info("\x1b[1;33mWRITING\x1b[0m: to out/%s.reg"%self.name[:-4])
            for s in self.sourcelist:
                outreg.write("circle(%fd, %fd, 10)\n"%(s.RA, s.DEC))

    def exportcat(self):
        if not os.path.isdir('out'):
            os.system('mkdir out/')
        outfile = "out/%s.sb"%self.name[:-4]
        with open(outfile, 'w') as outcat:
            outcat.write("#Catalog name: %s\n"%self.name)
            outcat.write("#Data Length: %s\n"%len(self.sourcelist))
            outcat.write("#Number Epochs/Tiles: %d\n"%self.sourcelist[0].size[0])
            outcat.write("#Number Bands: %d\n"%self.sourcelist[0].size[1])
            outcat.write("ID RA DEC ")
            for b in range(self.sourcelist[0].size[1]):
                outcat.write("FLUX:%d "%b )
                outcat.write("FLUXERR:%d "%b )
                outcat.write("MAG:%d "%b )
                outcat.write("MAGERR:%d "%b )
            for i,s in enumerate(self.sourcelist,1):
                outcat.write("\n%d %s"%(i, s.createExportString(style='mean tiles')))
        """
        with open("out/%s.sb"%self.name[:-4], 'w') as outcat:
            logging.info("\x1b[1;33mWRITING\x1b[0m: to out/%s.sb"%self.name[:-4])
            outcat.write("ID RA DEC")
            for e in range(1, self.sourcelist[0].numEpochs+1):
                for i in range(1, self.sourcelist[0].numBands+1):
                    #outcat.write("FLUX%d FLUXERR%d "%(i,i))
                    outcat.write(" flux[e%d:b%d] fluxerr[e%d:b%d]"%(e,i,e,i))
            for e in range(1, self.sourcelist[0].numEpochs+1):
                for i in range(1, self.sourcelist[0].numBands+1):
                    #outcat.write("MAG%d MAGERR%d "%(i,i))
                    outcat.write(" mag[e%d:b%d] magerr[e%d:b%d]"%(e,i,e,i))
            outcat.write('\n')
            for i, s in enumerate(self.sourcelist,1):
                outcat.write("%d %s\n"%(i,s.createExportString(full=True)))
        """

    def dustCorrection(self):
        U = 0.3543 #um
        G = 0.4770 #um
        R = 0.6231 #um

        Rv= 3.1

        Cu = self.a(U) + self.b(U)/Rv
        Cg = self.a(G) + self.b(G)/Rv
        Cr = self.a(R) + self.b(R)/Rv

        print(Cu, Cg, Cr)

        logging.info("dust correction")
        mainsequence = np.genfromtxt("catalogs/pleiades_sloane")
        mask = ( mainsequence[:,1]<0.7)
        mainsequence = mainsequence[mask]
        coeffs = np.polyfit(mainsequence[:,1], mainsequence[:,0],6)
        x=np.arange(-0.3,0.75,0.1)
        #y=np.polyval(coeffs, x)
        
        Av = np.arange(0,5,0.01)
        Chivals = np.zeros(Av.shape)
        for i, av in enumerate(Av):
            chi=0
            for s in self.sourcelist:
                #s.construct_colours([2,1],[1,0])
                #if(np.sum(s.resolved) == (s.size[0]*s.size[1])):
                s.colours[0]+=((Cu-Cg)*av)
                s.colours[1]+=((Cg-Cr)*av)
                if(np.isfinite(sum(s.colours))):
                    chi+= ((s.colours[0] - np.polyval(coeffs, s.colours[1]))/(1))**2.0
            Chivals[i] = chi
        print(Chivals)
        print(Av[np.argmin(Chivals)])
        plt.plot(Av, np.log(Chivals))
        plt.show()

                #plt.scatter(s.colours[1], -s.colours[0],c='k', s=0.1)
        
        #plt.show()


    def a(self, x):
        x-=1.82
        coeffs = [0.32999, -0.77530, 0.01979, 0.72085, -0.02427, -0.50447, 0.17699, 1.0]
        return np.polyval(coeffs, x)

    def b(self, x):
        x-=1.82
        coeffs = [-2.09002, 5.30260, -0.62251, -5.38434, 1.07233, 2.28305, 1.41338, 0]
        return np.polyval(coeffs, x)







    def cut_lowconfidence(self):
        self.sourcelist = [s for s in self.sourcelist if s.quality]
        

    def __repr__(self):
        return("\x1b[1;32m%s\x1b[0m"%self.name)
    



if __name__=='__main__':
    #cat = CATALOG(fitsfile="../test/ngc884_g_radec.fits", configfile='../config', catalog_style='sextractor', catalog_filename='../test/ngc884_g.cat')
    cat = CATALOG(catalog_style='starbug', catalog_filename='test/ngc884.sb')
    testdata = np.genfromtxt("test/dereddening_teststars.txt", skip_header=1)
    print(testdata)
    cat.sourcelist = cat.sourcelist[:len(testdata)]
    for i in range(len(cat.sourcelist)):
        cat.sourcelist[i].set_colours(testdata[i,2], testdata[i,0])
        cat.sourcelist[i].resolved = np.ones(np.shape(cat.sourcelist[i].resolved))
        print(cat.sourcelist[i].colours)

    cat.dustCorrection()
