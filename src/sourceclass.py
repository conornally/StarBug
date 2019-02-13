import numpy as np

class Source(object):
    #def __init__(self, ra=0, dec=0, mag=0, magerr=0, shp=0, rnd=0):
    def __init__(self,  ID=np.nan, 
                        x=np.nan, y=np.nan, 
                        ra=np.nan, dec=np.nan, 
                        flux=np.nan, fluxerr=np.nan, 
                        mag=np.nan, magerr=np.nan, 
                        shp=np.nan, rnd=np.nan,
                        bands=1, epochs=0):
        if not epochs: shape = np.ones((bands))
        else: shape = np.ones((bands, epochs))
        self.ID =       ID *shape
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

        self.numBands = 1
        self.numEpochs= 1

        self.set_means()
        self.quality = False
        self.resolved_in = 1
        

    def set_means(self):
        self.RA =       np.nanmean( self.ra)
        self.DEC=       np.nanmean( self.dec)
        if self.numEpochs > 1:
            self.FLUX=      np.nanmean(self.flux,0)
            #self.FLUXERR=   np.nanmean(self.fluxerr,0)
            self.MAG=       np.nanmean( self.mag,0)
            #self.MAGERR=    np.nanmean( self.magerr,0)
            self.SHARP=     np.nanmean(self.shp,0)
            self.ROUND=     np.nanmean(self.rnd,0)
            self.calc_errors()
        else:
            self.FLUX=      self.flux
            self.FLUXERR=   self.fluxerr
            self.MAG=       self.mag
            self.MAGERR=    self.magerr
            self.SHARP=     self.shp
            self.ROUND=     self.rnd

    def calc_errors(self):
        #across epochs
        df_dxi = 1./self.numEpochs
        dxi = lambda xi: (df_dxi**2.)*(xi**2.)
        for band in range(self.numBands):

            for epoch in range(self.numEpochs):
                pass



            self.FLUXERR[band] = np.sqrt(sum([ dxi(xi) for xi in self.fluxerr[:,band]]))
            self.MAGERR[band] = np.sqrt(sum([ dxi(xi) for xi in self.magerr[:,band]]))




    def append_Band(self, source=False, bad=False):
        """
        INPUT: instance of Source()
        FUNC: This will append an epoch of data onto each of the stats in self
        LIMITATIONS FOR NOW: BAND MUST BE DONE FIRST
        """
        #for now, just ra
        if bad or not source:
            if self.numEpochs==1: source = Source(bands=self.numBands, epochs=0)
            else: source = Source(bands=self.numBands, epochs=self.numEpochs)
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

        self.numBands +=1

    def append_Epoch(self, source=False, bad=False):
        """
        INPUT: instance of Source()
        FUNC: This will append an band of data onto each of the stats in self
        LIMITATIONS FOR NOW: BAND MUST BE DONE FIRST
        """

        if bad or not source:
            if self.numEpochs==1: source = Source(bands=self.numBands, epochs=0)
            else: source = Source(bands=self.numBands, epochs=self.numEpochs)
        else: self.resolved_in +=1
        self.x = np.stack((self.x, source.x))
        self.y = np.stack((self.y, source.y))
        self.ra = np.stack((self.ra, source.ra)) 
        self.dec = np.stack((self.dec, source.dec)) 
        self.flux = np.stack((self.flux, source.flux)) 
        self.fluxerr = np.stack((self.fluxerr, source.fluxerr)) 
        self.mag = np.stack((self.mag, source.mag)) 
        self.magerr = np.stack((self.magerr, source.magerr)) 
        self.shp = np.stack((self.shp, source.shp)) 
        self.rnd = np.stack((self.rnd, source.rnd)) 
    
        self.numEpochs +=1

    def append_Tile(self, source=False, bad=False):
        self.append_Epoch(source, bad) # for now it just does the same, but perhaps in the future ill want to change that
        
    def createExportString(self, full=False):
        self.set_means()
        exportstring="%f %f"%(self.RA, self.DEC)

        if full:
            true = lambda val: val if not np.isnan(val) else 9999.9999
            for epoch in range(self.numEpochs):
                for band in range(self.numBands):
                    if((self.numEpochs>1) and (self.numBands>1)):
                        exportstring += " %f %f"%(true(self.flux[epoch, band]), true(self.fluxerr[epoch, band]))
                    elif(self.numEpochs>1):
                        exportstring += " %f %f"%(true(self.flux[epoch]), true(self.fluxerr[epoch]))
                    elif(self.numBands>1):
                        exportstring += " %f %f"%(true(self.flux[band]), true(self.fluxerr[band]))

            for epoch in range(self.numEpochs):
                for band in range(self.numBands):
                    if((self.numEpochs>1) and (self.numBands>1)):
                        exportstring += " %f %f"%(true(self.mag[epoch, band]), true(self.magerr[epoch, band]))
                    elif(self.numEpochs>1):
                        exportstring += " %f %f"%(true(self.mag[epoch]), true(self.magerr[epoch]))
                    elif(self.numBands>1):
                        exportstring += " %f %f"%(true(self.mag[band]), true(self.magerr[band]))
        else:
            for band in range(self.numBands):
                exportstring += " %f %f"%(self.FLUX[band], self.FLUXERR[band])
            for band in range(self.numBands):
                exportstring += " %f %f"%(self.MAG[band], self.MAGERR[band])
        return exportstring

    def __repr__(self):
        #return("\x1b[1;32m%s\x1b[0m: %.4f %.4f %.4f"%(self.ID, self.RA[0], self.DEC[0], self.MAG[0]))
        return("\x1b[1;32m{}\x1b[0m: {} {} {} {}".format(self.ID, self.RA, self.DEC, self.MAG, self.MAGERR))

if __name__=='__main__':
    s1 = Source(ra=0,dec=0,flux=11, fluxerr=1)
    s2 = Source(ra=0,dec=0,flux=12, fluxerr=1)
    s3 = Source(ra=0,dec=0,flux=21, fluxerr=2)
    s4 = Source(ra=0,dec=0,flux=22, fluxerr=2)
    s1.append_Band(s2)
    print("AfterBand before Epoch", s1.flux)
    print(s1.createExportString(full=True))
    s3.append_Band(s4)
    s1.append_Epoch(s3)
    
    print("AfterBand after Epoch", s1.flux)
    print(s1.createExportString(full=True))
    s1.set_means()
