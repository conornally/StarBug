from fitsclass import FITS
import numpy as np
from astropy.io import fits

"""
dark frame reduction
>> need to stack multple frames of this together
flat field multiplication
>> need to stack multple frames of this together
"""



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
    d1 = FITS('../test/dark1.fits')
    d2 = FITS('../test/dark2.fits')
    d1 = darkFrame_build([d1,d2])
    
    fl1 = FITS('../test/flat1.fits')
    fl2 = FITS('../test/flat2.fits')
    fl1 = flatFrame_build([fl1,fl2])

    f1 = FITS('../test/frame1.fits')

    f1.subtract(d1)
    f1.divide(fl1)
    save(f1)
