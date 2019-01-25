import numpy as np, logging, time
np.warnings.filterwarnings('ignore')
from scipy import optimize

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

        self.raw_data = self.construct_raw_data(catalog_style, catalog_filename, catalog)
        self.sourcelist = np.empty( (len(self.raw_data)), dtype=object)
        self.build_sourcelist()

        name = catalog_filename
        while '/' in name: name = name[name.index('/')+1:]
        self.name = name

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
            ID = int(get_value('IDcol', configfile=self.configfile))
            flux = int(get_value('Fluxcol', configfile=self.configfile))
            fluxerr = int(get_value('FluxErrcol', configfile=self.configfile))
            mag = int(get_value('Magcol', configfile=self.configfile))
            magerr = int(get_value('MagErrcol', configfile=self.configfile))
            ra = int(get_value('RAcol', configfile=self.configfile))
            dec = int(get_value('DECcol', configfile=self.configfile))
            x = int(get_value('Xcol', configfile=self.configfile))
            y = int(get_value('Ycol', configfile=self.configfile))
            skip_header = int(get_value('skip_header', configfile=self.configfile))
            skip_footer = int(get_value('skip_footer', configfile=self.configfile))
            comments = int(get_value('Comments', configfile=self.configfile))
            delimiter = int(get_value('delimiter', configfile=self.configfile))

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

    def setup(self):
        """
        __init__ second round, without filenames
        """
        self.cat = np.zeros((len(self.sourcelist), 2))
        for i,s in enumerate(self.sourcelist):
            s.set_means()
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
            #self.sourcelist[i] = Source( l[1], l[2], l[3], l[4], l[8], l[7])
            self.sourcelist[i] = Source( *l )


    def testing(self):
        start = time.time()
        for s in self.sourcelist:
            s.ra*=2
        print(time.time() - start)
        start = time.time()
        self.raw_data[:,1]*=2
        print(time.time() - start)


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

    def calibrate(self, units='Dns-1_mJy'):
        """
        >Converts DN/s to mJy //needs to be generalised
        >Applies colour correction
        >Applies zero point correction
        // currently only setup for spitzer
        """
        logging.info('CALIBRATING: %s'%self)

        DN2MJy = self.fits.header['FLUXCONV'] / self.fits.header['EXPTIME']
        MJy2mJy = (1e9/4.254517e10)*(self.fits.header['PXSCAL1']**2.)
        factor=1. 
        if units =='Dns-1_mJy':
            factor = DN2MJy * MJy2mJy 

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
            s.set_means()   #this sloooows things down a lot 
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



    def BandMatch(self, ALS):
        logging.info("\x1b[1;33mMATCHING\x1b[0m Bands: %s <-- %s"%(self, ALS))
        for s in self.sourcelist: s.set_means()
        self.match(ALS, regime='band')
        
    def EpochMatch(self, ALS):
        logging.info("\x1b[1;33mMATCHING\x1b[0m Epochs: %s <-- %s"%(self, ALS))
        for s in self.sourcelist: s.set_means()
        self.match(ALS, regime='epoch')

    def match(self, ALS, regime='epoch'):
        """
        INPUT: >List or single instance of ALS_DATA(), 
               >regime is 'band' or 'epoch' matching
               >> for now, i think all epoch matching should be done first
        FUNC: adds additional data onto the sources in l=source list
        (should i retake averages inbetween matches?)
        this will be a loooot slower though, but maybe more representative
        """
        if type(ALS) == ALS_DATA: ALS = [ALS]
        
        for cat in ALS:
            logging.info('...%s'%cat)
            radec = np.array( [ [s.RA, s.DEC ] for s in self.sourcelist])
            radec = SkyCoord( ra=Angle( radec[:,0], unit=u.deg), dec=Angle( radec[:,1], unit=u.deg) )
            for s in cat.sourcelist: s.set_means()
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
                        self.sourcelist[i].quality= True
                #If i dont crop out the bad sources i will NEED to append 9999 values to sourcelist objects..
            self.sourcelist = [s for s in self.sourcelist if s.quality]
            for s in self.sourcelist: s.set_means()
            logging.info("Sources: %d"%len(self.sourcelist))

    def single_match(self, cat):
        return match_coordinates_sky(self.cat, cat)

    def func(self, offsets, ALS=None):
        offset = SkyCoord( ra=offsets[0]*np.ones((ALS.size))*u.deg, dec=offsets[1]*np.ones((ALS.size))*u.deg)
        #tmpcat = SkyCoord( ra=Angle( [ra]*ALS.size, unit=u.deg), dec=Angle( [dec]*ALS.size, unit=u.deg))
        #tmpcat.ra += ALS.cat.ra
        
        tmpcat = SkyCoord( ra=ALS.cat.ra+offset.ra, dec=ALS.cat.dec+offset.dec)
        idx, d2d, d3d = self.single_match(tmpcat)
        bad = np.where(d2d.arcsec > 1)[0]
        badlen = len(bad)
        for i, IDX in enumerate(idx):
            if i not in bad:
                if d2d[i] != min(d2d[np.where( idx==IDX)]):
                    badlen+=1
        print(offsets, self.size, badlen)
        return(badlen)


    def align(self, ALS):
        """
        will need to pull out all the star data into an array, otherwise scipy will take hours
        match function will need to be v.fast
        """
        #out =optimize.minimize( self.func, [0.0,0.0], args=ALS) 
        ra = self.fits.header['CRVAL1']
        dec = self.fits.header['CRVAL2']
        print(ra,dec)
        perc=0.0001
        out = optimize.brute( self.func, [[ra-perc*ra, ra+perc*ra], [dec-perc*dec, dec+perc*dec]], args=[ALS])
        print(out)
        # print(self.size - out['fun'])





    def __repr__(self):
        return("\x1b[1;32m%s\x1b[0m"%self.name)
    



if __name__=='__main__':
    cat = CATALOG(fitsfile="../test/ngc884_g_radec.fits", configfile='../config', catalog_style='sextractor', catalog_filename='../test/ngc884_g.cat')
    print(cat)