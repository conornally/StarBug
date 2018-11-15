import sys
import numpy as np
from astropy.io import fits

img = fits.open(sys.argv[1])
d = img[0].data
h = img[0].header
img.close()

new = np.array(d, dtype='float32')
out = fits.PrimaryHDU(data = new, header=h)
out.writeto('out32.fits')

