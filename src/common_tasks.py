import os
from fitsclass import FITS
from fileio import *
import numpy as np
from astropy.io import fits

"""
dark frame reduction
>> need to stack multple frames of this together
flat field multiplication
>> need to stack multple frames of this together
"""

def pre_adjustments( fitslist ):
    """
    FUNC: conducts cropping, scaling, naming of fits files beforethey are put though daophot or photutils
    """
    if not os.path.isdir('out'): os.system('mkdir out')    
    if type(fitslist) == FITS: fitslist = [fitslist]
    for fits in fitslist:
        fits.scale()
        fits.crop()
        fits.export("out/%s"%fits.name, overwrite=True)


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
    
    if type(flats)==FITS: 
        flats.normalise()
        return flats
    else:
        for f in flats[1:]:
            flats[0].add(f)
        flats[0].normalise()
        return flats[0]
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
        print(np.nanmedian(frame.data))
    flats[0].add_median(flats[1:])
    flats[0].normalise('median')
        

    return(flats[0])
    


def save(f):
    """this is tmporary, ultimately it will overwrite itself"""
    hdu = fits.PrimaryHDU(f.data)
    hdu.writeto('tmp', overwrite=True)

if __name__=='__main__':
    #pre_adjustments( fitsfromtxt('test/files.txt'))
    f1 = FITS('../test/flat1.fits')
    f2 = FITS('../test/flat2.fits')
    f3 = FITS('../test/flat3.fits')
    flatFrame_build([f1,f2, f3])
    
    frame1=FITS('../test/frame1.fits')
    frame1.divide(f1)
    save(frame1)
