import os, sys, logging
import numpy as np
import astropy.io.fits as fits
import astropy.stats as stats
import photutils
import parse_config

logging.basicConfig(level='DEBUG')#, format="\x1b[1;%dm" % (32) + '%(message)s' + "\x1b[0m")

class FITS(object):
    def __init__(self, filename):
        img = fits.open(filename)
        # i think i ned to get support for multiple frames
        self.header = img[0].header
        self.data = img[0].data
        self.size=np.shape(self.data)
        self.dtype = self.data.dtype
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

        self.options={'fwhm':0,'threshold':0,'sigma':0,'roundness':[0,0],'sharpness':[0,0]}



    #################################
    # Source Detection and Analysis #
    #################################

    def find(self):
        """InPUT:   
            FUNC:   Daofind routine to do initial pass on source detection
        """
        if self.options['threshold'] ==0: logging.warning('Threshold at default: 0, Change this for source detection')
        else:
            print(self.options)
            daofind = photutils.DAOStarFinder(  fwhm=self.options['fwhm'],
                                                threshold=self.options['threshold'],
                                                sharplo=self.options['sharpness'][0],
                                                sharphi=self.options['sharpness'][1],
                                                roundlo=self.options['roundness'][0],
                                                roundhi=self.options['roundness'][1])
            self.sourcelist = daofind(self.data)
            #this needs to be fixed!
            print(self.sourcelist)
            self.sourcelist.remove_rows( np.where( self.sourcelist['flux'] < 50))
            print(self.sourcelist)
            with open('tmp.reg','w') as reg:
                for line in self.sourcelist:
                    #reg.write('circle({}, {}, {})\n'.format(line['xcentroid'], line['ycentroid'], self.options['fwhm']))
                    reg.write('{} {}\n'.format(line['xcentroid'], line['ycentroid']))


    def get_offset(self, fitsobj):
        """INPUT:   instance of FITS
            FNUC:   gets pixel offset between fits objects
        """
        self.get_fft()
        fitsobj.get_fft()
        convolution = np.multiply(self.fft, np.conj(fitsobj.fft))
        inverse = np.fft.ifft2(convolution)
        print('argmax',np.argmax(inverse))
        print('shape',np.shape(inverse))
        print('argmax/shape0',np.argmax(inverse)/np.shape(self.data)[1])
        print('argmax0',np.argmax(inverse, 0))


    def get_fft(self):
        self.fft = np.fft.fft2(self.data)


    #############################
    # Pixel Array Manipulations #
    #############################

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
            self.data = np.multiply(self.data, fitsobj.data, dtype=self.data.dtype)
        else: self.data = np.multiply(self.data, factor, dtype=self.data.dtype)

    def divide(self, factor):
        """
        FUNC: divide self.data by fitsobj.data (element-wise)
        
        """
        logging.info('IMPORTANT, THIS MIGHT HAVE CAUSED ROUNDING!!')
        if type(factor) == FITS:
            self.data = np.divide(self.data, factor.data, dtype=float)
        else: self.data = np.divide(self.data, factor, dtype=float)

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
        self.data = np.divide(self.data, factor, dtype=self.data.dtype)

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

    def convert_dtype(self, dtype='float32'):
        """INPUT:   dtype, data type to convert self.data to
            FUNC:   changes data type of pixel array
        """
        if dtype in ['float16','float32','float64']:
            self.data = self.data.astype(dtype)
            logging.debug('--> (%s) %s'%(dtype, self))
            
        else: logging.info('Unknown data type')



    #######################
    # File I/O and config #
    #######################
    
    def load_options(self):
        """FUNC:loads in config file and stores options
        """
        
        config = parse_config.load()
        self.options['fwhm'] = parse_config.get_value('FWHM', config, dtype=float)
        self.options['threshold'] = parse_config.get_value('Threshold', config, dtype=float)
        self.options['sigma'] = parse_config.get_value('Sigma', config, dtype=float)
        self.options['sharpness'] = [parse_config.get_value('SharpLow', config, dtype=float), parse_config.get_value('SharpHigh', config, dtype=float)]
        self.options['roundness'] = [parse_config.get_value('RoundLow', config, dtype=float), parse_config.get_value('RoundHigh', config, dtype=float)]


    def display(self):
        try:
            os.system('ds9 %s -regions tmp.reg -zoom to fit'%self.filename)
        except: print('nope')

    def export(self, filename='', overwrite=False):
        if filename=='': 
           filename = self.filename
        logging.info('\x1b[1;33mEXPORTING\x1b[0m: %s --> %s'%(self, filename))
        fits.PrimaryHDU(data=self.data, header=self.header).writeto(filename, overwrite=overwrite)
        
    def __repr__(self):
        return("\x1b[1;32m%s\x1b[0m"%self.name)


if __name__=='__main__':
    f1 = FITS("../test/ngc869.fits")
    f2 = FITS("../test/ngc2.fits")
    #f1.load_options()
    #f1.options['fwhm']=20
    #f1.options['threshold']=100
    #f1.options['sharpness']=[0.3,0.7]
    #f1.convert_dtype('float32')
    #f1.find()
    #f1.display()
    #f1.get_offset(f2)
    #f1.add(f2)
    f1.display()
    f2.display()
    #f1.convert_dtype('float64')
    #f1.convert_dtype('float16')
    #f2 = FITS("../test/frame2.fits")
    #f3 = FITS("../test/frame3.fits")
    #fake=FITS("../test/frame1.fits")
    #fake.data = np.ones((2,2))
    #fake.name='FAKE'
    #f1.add_median([f2, f3, fake])
    #f1.unit_scale()
    #f1.zero_to_nan()
    #f1.normalise('median')
    #f1.multiply(f2)
    #f1.divide(f2)
    #f1.basic_stats(300,1)
    #f1.crop()
    # f1.scale()
    #f1.export('tmp.fits', True) 
