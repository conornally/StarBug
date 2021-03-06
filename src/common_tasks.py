import os
try:from fitsclass import FITS
except:from src.fitsclass import FITS
try:from fileio import *
except:from src.fileio import *
import numpy as np
import logging
from astropy.io import fits

###################
# Image Reduction #
###################

def pre_adjustments( fitslist ):
    """
    FUNC: conducts cropping, scaling, naming of fits files beforethey are put though daophot or photutils
    """
    if not os.path.isdir('out'): os.system('mkdir out')    
    if type(fitslist) == FITS: fitslist = [fitslist]
    for fits in fitslist:
        fits.unit_scale()
        fits.crop()


def darkFrame_build(darks):
    """
    INPUT: this is a list or single instance of dark frames
    FUNC: this return return the add_mean of the frames
    """
    if type(darks)==FITS: return darks
    else:
        darks[0].add_mean(darks[1:])
        return darks[0]

def flatFrame_build(flats, dark=False):
    """
    INPUT: this is a list or single instance of flats
    FUNC: returns a noralised addition of all the frames
    """
    #for now im going to do it in the case that there are multiple flat frames, i will generalise to one frame after
    medians = []
    for frame in flats:
        if dark: frame.subtract(dark)
        frame.median = np.nanmedian(frame.data)
        medians.append(frame.median)
    true_median = np.median(medians)
    for frame in flats:
        frame.multiply(true_median/frame.median)
    flats[0].add_median(flats[1:])
    flats[0].normalise('median')

    return(flats[0])

def darkFrame_subtract(fitslist, dark):
    """
    INPUT: Fits list and Fits instance as dark frame
    FUNC: Subtracts the fits pixel array from each of the fits files
    """
    if type(dark)==list: dark = dark[0]
    for fits in fitslist:
        logging.debug(fits)
        fits.subtract(dark)

def flatField_divide(fitslist, flat):
    """
    INPUTS: list of Fits instances and flat field fits instance
    FUNC: divides each pixel array in fitslist by the flat field
    """
    if type(flat)==list: flat=flat[0]
    for fits in fitslist:
        logging.debug(fits)
        fits.divide(flat)

def align(fitslist):
    """INPUT: fits instances
        FUNC: calculates offset between fitsimages
    """
        
    fitslist[0].get_fft()
    for i,fits in enumerate(fitslist[1:], start=1):
        fits.get_fft()
        fits.offset = fitslist[i-1].get_offset(fits) + fitslist[i-1].offset
        logging.info('%s offset %s'%(fits, fits.offset))
        if i >=2: fitslist[i-2].del_fft()
        
def exportOFFSET(fitslist, filename='auto'):
    """INPUT: list of single instance of FIST
              filename='auto' if derived from first FITS in list
              filename='path/to/file.offset' otherwise
        FUNC: saves list of offsets to filename
    """
    if type(fitslist)==FITS: fitslist=[fitslist]
    if filename=='auto': filename=fitslist[0].filename[:-5]
    filename+= str(".offset")
    with open(filename,'w')as offsetfile:
        for fits in fitslist:
            offsetfile.write("%s %d %d\n"%(fits, fits.offset[0], fits.offset[1]))

def readinOFFSET(fitslist, filename):
    """INPUT: offset filename
    """
    if type(fitslist)==FITS: fitslist=[fitslist]
    offset = np.genfromtxt(filename, usecols=(1,2), delimiter=' ')
    for i, fits in enumerate(fitslist):
        fits.offset = offset[i] 
        logging.info("%s <- [%d %d]"%(fits, offset[i,0], offset[i,1]))


     



#######################
# Array Manipulations #
#######################



############
# Analysis #
############

def basic_stats( fitslist, sigma=5, iters=3):
    """INPUT: list or single fits instant
                sigma clipping value
                clipping iterations
        Func:   Clips sigma value and returns mean, median stdev of fits image
    """
    if type(fitslist)==FITS: fitslist = [fitslist]
    for fits in fitslist:
        fits.basic_stats(sigma, iters)


###########
# File io #
###########

def dtype_convert(fitslist, dtype='float32'):
    if type(fitslist)==FITS: fitslist=[fitslist]
    for fits in fitslist:
        fits.convert_dtype(dtype)



def save(f):
    """this is tmporary, ultimately it will overwrite itself"""
    hdu = fits.PrimaryHDU(data=f.data, header=f.header)
    hdu.writeto('tmp.fits', overwrite=True)

if __name__=='__main__':
    f1 = FITS('../test/frame1.fits')
    f2 = FITS('../test/frame2.fits')
    f3 = FITS('../test/frame3.fits')
    align([f1,f2, f3])
    
