import numpy as np

class Source(object):
    def __init__(self,  ID=np.nan, 
                        x=np.nan, y=np.nan, 
                        ra=np.nan, dec=np.nan, 
                        flux=np.nan, fluxerr=np.nan, 
                        mag=np.nan, magerr=np.nan, 
                        shp=np.nan, rnd=np.nan,
                        bands=1, epochs=1):
        shape = np.ones((epochs, bands))
        self.ID =       ID
        self.x =        x*shape
        self.y =        y*shape
        self.ra=        ra*shape
        self.dec=       dec*shape
        self.flux =     flux*shape
        self.fluxerr=   fluxerr*shape
        self.mag=       mag*shape
        self.magerr=    magerr*shape
        self.shp=       shp*shape
        self.rnd=       rnd*shape

        self.resolved = self.mag / self.mag

        self.colours = np.zeros(2)
        self.colourError = np.zeros(2)
        self.size = np.shape(self.mag)

        self._voidCalcTileMeans()
        self.quality = False
        self.resolved_in = 1
        

    def _voidCalcTileMeans(self):
        self.RA =       np.nanmean( self.ra)
        self.DEC=       np.nanmean( self.dec)
        self.FLUX=      np.nanmean(self.flux,0)
        self.FLUXERR=   np.zeros(np.shape(self.FLUX)[0])
        self.MAG=       np.nanmean( self.mag,0)
        self.MAGERR=    np.zeros(np.shape(self.MAG)[0])
        self.SHARP=     np.nanmean(self.shp,0)
        self.ROUND=     np.nanmean(self.rnd,0)
        self.calc_errors()

    def calc_errors(self):
        #across epochs
        for band in range(self.size[1]):
            N = len(np.where( np.isfinite( self.fluxerr[:,band]))[0])
            self.FLUXERR[band] = np.sqrt( np.nansum( [ (xi/N)**2. for xi in self.fluxerr[:,band]]))
            N = len(np.where( np.isfinite( self.magerr[:,band]))[0])
            self.MAGERR[band] =  np.sqrt( np.nansum( [ (xi/N)**2. for xi in self.magerr[:,band]]))
            


    def get_Value(self, key='flux', epoch=0, band=0):
        returnVal = 0
        if(hasattr(self, key)):
            if( (np.shape(getattr(self, key))[0] > epoch) and (np.shape(getattr(self,key))[1] > band)):
                returnVal = getattr(self, key)[epoch, band]
        return returnVal 

    def construct_colours(self, c1=[0,1], c2=[1,2]):
        self._voidCalcTileMeans()
        self.colours[0] = self.MAG[c1[0]] - self.MAG[c1[1]]
        self.colours[1] = self.MAG[c2[0]] - self.MAG[c2[1]]
        self.colourError[0] = np.sqrt( self.MAGERR[c1[0]]**2. +self.MAGERR[c1[1]]**2. )
        self.colourError[1] = np.sqrt( self.MAGERR[c2[0]]**2. +self.MAGERR[c2[1]]**2. )


    def set_colours(self, colour1=0, colour2=0):
        self.colours[0] = colour1
        self.colours[1] = colour2

    def append_Band(self, source=False, bad=False):
        """
        INPUT: instance of Source()
        FUNC: This will append an epoch of data onto each of the stats in self
        """
        #for now, just ra
        if bad or not source:
            source = Source(bands=1, epochs=self.size[0])
            source.resolved *=False
        else: self.resolved_in +=1

        self.x = np.concatenate((self.x, source.x),1)
        self.y = np.concatenate((self.y, source.y),1)
        self.ra = np.concatenate((self.ra, source.ra),1)
        self.dec = np.concatenate((self.dec, source.dec),1)
        self.flux = np.concatenate((self.flux, source.flux),1)
        self.fluxerr = np.concatenate((self.fluxerr, source.fluxerr),1)
        self.mag = np.concatenate((self.mag, source.mag),1)
        self.magerr = np.concatenate((self.magerr, source.magerr),1)
        self.shp = np.concatenate((self.shp, source.shp),1)
        self.rnd = np.concatenate((self.rnd, source.rnd),1)

        self.resolved = np.concatenate((self.resolved, source.resolved ),1)
        self.size = np.shape(self.flux)
        self._voidCalcTileMeans()

    def append_Epoch(self, source=False, bad=False):
        """
        INPUT: instance of Source()
        FUNC: This will append an band of data onto each of the stats in self
        """

        if bad or not source:
            source = Source(bands=self.size[1], epochs=1)
            source.resolved*=False
        else: self.resolved_in +=1
        self.x = np.concatenate((self.x, source.x),0)
        self.y = np.concatenate((self.y, source.y),0)
        self.ra = np.concatenate((self.ra, source.ra),0) 
        self.dec = np.concatenate((self.dec, source.dec),0) 
        self.flux = np.concatenate((self.flux, source.flux),0) 
        self.fluxerr = np.concatenate((self.fluxerr, source.fluxerr),0) 
        self.mag = np.concatenate((self.mag, source.mag),0) 
        self.magerr = np.concatenate((self.magerr, source.magerr),0) 
        self.shp = np.concatenate((self.shp, source.shp),0) 
        self.rnd = np.concatenate((self.rnd, source.rnd),0) 
    
        self.resolved = np.concatenate((self.resolved, source.resolved),0)
        self.size = np.shape(self.flux)
        self._voidCalcTileMeans()

    def append_Tile(self, source=False, bad=False):
        self.append_Epoch(source, bad) # for now it just does the same, but perhaps in the future ill want to change that
        
    def createExportString(self, style='full'):
        exportstring="%f %f "%(self.RA, self.DEC)
        if(style=='full'):
            for epoch in range(self.size[0]):
                for band in range(self.size[1]):
                    exportstring += "%f "%self.get_Value('flux', epoch, band)
                    exportstring += "%f "%self.get_Value('fluxerr', epoch, band)
                    exportstring += "%f "%self.get_Value('mag', epoch, band)
                    exportstring += "%f "%self.get_Value('magerr', epoch, band)
        zero2nan = lambda x: np.nan if x==0 else x
        if(style=='mean tiles'):
            for band in range(self.size[1]):
                exportstring += "%f "%zero2nan(self.FLUX[band])
                exportstring += "%f "%zero2nan(self.FLUXERR[band])
                exportstring += "%f "%zero2nan(self.MAG[band])
                exportstring += "%f "%zero2nan(self.MAGERR[band])


        return exportstring

    def __getitem__(self, listRef):
        key, epoch, band = listRef
        return self.get_Value(key, epoch, band)
            


    def __repr__(self):
        self._voidCalcTileMeans()
        return("\x1b[1;32m{}\x1b[0m: {} {} {} {}".format(self.ID, self.RA, self.DEC, self.MAG, self.MAGERR))

if __name__=='__main__':
    s1 = Source(ra=0,dec=0,flux=11, fluxerr=1, mag=2, magerr=1)
    s2 = Source(ra=0,dec=0,flux=12, fluxerr=1, mag=1, magerr=1)
    s3 = Source(ra=0,dec=0,flux=21, fluxerr=3, mag=0, magerr=1)

    s = Source( mag=np.array([1,np.nan]), bands=2, epochs=1)
    print(s.resolved)

    s1.append_Band(s2)
    s1.append_Band(s3)
    s1.append_Band(bad=True)
    s1.append_Epoch(bad=True)
    print(np.sum(s1.resolved))

    s1.construct_colours([2,1],[1,0])
    print(s1.colours)

