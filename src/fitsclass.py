import sys, logging
import numpy as np
import astropy.io.fits as fits
import astropy.stats as stats
import matplotlib.pyplot as plt

logging.basicConfig(level='DEBUG')#, format="\x1b[1;%dm" % (32) + '%(message)s' + "\x1b[0m")

class FITS(object):
    def __init__(self, filename):
        img = fits.open(filename)
        # i think i ned to get support for multiple frames
        self.header = img[0].header
        self.data = img[0].data
        self.size=np.shape(self.data)
        img.close()

        self.filename = filename
        self.name = filename
        while('/' in self.name):
            self.name = self.name[self.name.index('/')+1:]
        logging.info('\x1b[1;33mINPUT\x1b[0m: %s'%self)

        #initialisaing yet unknown stats
        self.mean=0
        self.median=0
        self.std=0


    def add(self, fitsobj):
        """
        INPUT: another FITS instance
        FUNC: Adds the two datas together in self
        >> this will not change any values in fitclass instance
        """
        if np.shape(self.data) == np.shape(fitsobj.data):
            self.data = np.add(self.data, fitsobj.data)
            return 1
        else:
            logging.debug('%s data has different size to %s data'%(fitsobj, self))
            return 0

    def subtract(self, fitsobj):
        """
        INPUT: another fits instance
        FUNC: subtarcts fitsobj from self
        """
        if np.shape(self.data) == np.shape(fitsobj.data):
            self.data = np.subtract(self.data, fitsobj.data)
            return 1
        else: 
            logging.debug('%s data has different size to %s data'%(fitsobj, self))
            return 0

    def multiply(self, factor):
        """
        FUNC: multiplies self.data by fitsobj.data (element-wise)
        """
        if type(factor) == FITS:
            self.data = np.multiply(self.data, fitsobj.data, dtype='float64')
        else: self.data = np.multiply(self.data, factor, dtype='float64')

    def divide(self, factor):
        """
        FUNC: divide self.data by fitsobj.data (element-wise)
        """
        if type(factor) == FITS:
            self.data = np.divide(self.data, factor.data, dtype='float64')
        else: self.data = np.divide(self.data, factor, dtype='float64')

    def add_mean(self, fitsobj):
        """
        INPUT: list or single instance of FITS classes
        FUNC: Adds fitsdata together, then divides result by two
        // there is no safety net on this, if one fitsobj doesnt get added, then there is no warning
        """
        if type(fitsobj) == FITS: fitsobj = [fitsobj]
        count=1.0 #this will breake if fitsobj == None
        for f in fitsobj:
            if self.add(f): 
                count +=1.0
                logging.debug('%s contributing to add_mean on %s'%(f,self))

        self.divide(count)

    def add_median(self, fitsobj):
        """
        INPUT: list or single instance of FITS class
        FUNC: changes self.data pixels to be the median for each pixel over the frames
        """
        if type(fitsobj) == np.ndarray:
            fitsobj = np.insert(fitsobj, 0, self)
        if type(fitsobj) == FITS: fitsobj = np.array([self, fitsobj])
        if type(fitsobj) == list: 
            fitsobj.insert(0, self)
            fitsobj = np.array(fitsobj)
        for i, f in enumerate(fitsobj):#this is currently bad
            if np.shape(f.data) != np.shape(self.data):
                fitsobj = np.delete(fitsobj, i)
        fitsdata= np.array( [ f.data for f in fitsobj] )
        '''for i in range(len(self.data)):
            for j in range(len(self.data[0])):
                #self.data[i,j] = np.median(fitsdata[:,i,j])
                self.data[i,j] = np.median(fitsdata[:,i,j])
            print(i)
        print(self.data) 
        '''
        logging.debug("Taking median pixel values across: %s"%fitsobj)
        self.data = np.median(fitsdata[:,],0)
        

    def normalise(self, scale='max'):
        """
        INPUT: scale = scaling factor type
              >max = (default)
        FUNC: Normalises self.data to be between 0 (-1) and 1
        """
        logging.info("\x1b[1;33mNORMALISING\x1b[0m: %s"%self)
        if scale=='max': factor = np.nanmax(self.data)
        elif scale=='median': factor = np.nanmedian(self.data)
        self.data = np.divide(self.data, factor, dtype=np.float64)

    def sigma_clip(self, sigma):
        """
        INPUT: sigma value, above which will be cut
        This will be useful for flat fielding, as sometimes there are sources in the image, so i can cut them out here
        """
        pass

    def basic_stats(self, sigma, iters=3):
        logging.debug("%s Taking basic stats with sigma %f"%(self, sigma))
        self.mean, self.median, self.std = stats.sigma_clipped_stats( self.data, sigma=sigma, iters=iters )
        logging.debug("%s: \n   Mean: %f\n Median: %f\n    Std: %f\n    Min: %f\n    Max: %f"%(self, self.mean, self.median, self.std, np.min(self.data), np.max(self.data)))

    
    def unit_scale(self):
        """
        THIS IS NOT GENERAL YET
        It scales the units of the fits file from MJysr-1 to DNs-1
        """
        factor = float(self.header['EXPTIME']) / float(self.header['FLUXCONV'])
        logging.info("\x1b[1;33mSCALING\x1b[0m: %s *= %f"%(self, factor))
        self.header['BUNIT'] = 'DN/s'
        self.multiply(factor)

    def crop(self, auto=True):
        """
        INPUT: automatic cropping, will cut as close to FOV as possible
               IF not auto, then it will prompt inputs for top left and bottom right pixels
        FUNC: Crops the fits data value to ~match FOV
              Used if there is unwanted inf vaues surrounding FOV
        """
        if auto:
            mask = np.isfinite(self.data)

            left=0
            while 1 not in mask[:,left]: left+=1
            right=self.size[1]-1
            while 1 not in mask[:,right]: right-=1
            top=0
            while 1 not in mask[top]: top+=1
            bottom=self.size[0]-1
            while 1 not in mask[bottom]: bottom-=1

        else:
            left = int(input('left: '))
            right = int(input('right: '))
            top = int(input('top: '))
            bottom = int(input('bottom: '))
        
        sys.stdout.write('\x1b[1;31mCRPIX:\x1b[0m [{}, {}]'.format(self.header['CRPIX1'], self.header['CRPIX2'])) 
        self.header['CRPIX1']-=float(left)
        self.header['CRPIX2']-=float(top)
        sys.stdout.write(" --> [{}, {}]\n".format(self.header['CRPIX1'], self.header['CRPIX2']))
        self.data = self.data[top:bottom, left:right]
        self.size = np.shape(self.data)

    def nan_to_zero(self):
        """ FUNC: Converts nan values in data to int 0 """
        self.data[np.where( np.isnan(self.data))] =0

    def zero_to_nan(self):
        """ FUNC: Converts 0 in data to np.nan""" 
        self.data[np.where( self.data==0) ] = np.nan


    def export(self, filename='', overwrite=False):
        if filename=='': 
           filename = self.filename
        logging.info('\x1b[1;33mEXPORTING\x1b[0m: %s --> %s'%(self, filename))
        fits.PrimaryHDU(data=self.data, header=self.header).writeto(filename, overwrite=overwrite)
        
    def __repr__(self):
        return("\x1b[1;32m%s\x1b[0m"%self.name)


if __name__=='__main__':
    f1 = FITS("../test/raw1.fits")
    #f2 = FITS("../test/frame2.fits")
    #f3 = FITS("../test/frame3.fits")
    #fake=FITS("../test/frame1.fits")
    #fake.data = np.ones((2,2))
    #fake.name='FAKE'
    #f1.add_median([f2, f3, fake])
    f1.unit_scale()
    f1.zero_to_nan()
    f1.normalise('median')
    #f1.multiply(f2)
    #f1.divide(f2)
    #f1.basic_stats(300,1)
    #f1.crop()
    # f1.scale()
    #f1.export('tmp.fits', True) 
