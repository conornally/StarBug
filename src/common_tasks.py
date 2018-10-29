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

def flatFrame_build(flats):
    """
    INPUT: this is a list or single instance of flats
    FUNC: returns a noralised addition of all the frames
    """
    if type(flats)==FITS: 
        flats.normalise()
        return flats
    else:
        for f in flats[1:]:
            flats[0].add(f)
        flats[0].normalise()
        return flats[0]


def save(f):
    """this is tmporary, ultimately it will overwrite itself"""
    hdu = fits.PrimaryHDU(f.data)
    hdu.writeto('tmp', overwrite=True)

if __name__=='__main__':
    pre_adjustments( fitsfromtxt('test/files.txt'))
