import numpy as np

class Source(object):
    def __init__(self, ra=0, dec=0, mag=0, magerr=0, shp=0, rnd=0):
        self.ra= np.array([ra])
        self.dec= np.array([dec])
        self.mag= np.array([mag], dtype=np.float64)
        self.magerr= np.array([magerr], dtype=np.float64)
        self.shp= np.array([shp])
        self.rnd= np.array([rnd])

        #self.set_means()
        self.quality = False

    def set_means(self):
        self.RA = np.nanmean( self.ra )
        self.DEC= np.nanmean( self.dec)
        self.MAG= np.nanmean( self.mag)
        self.MAGERR= np.nanmean( self.magerr)
        self.SHARP= np.nanmean(self.shp)
        self.ROUND= np.nanmean(self.rnd)

    def append_Epoch(self, source):
        """
        INPUT: instance of Source()
        FUNC: This will append an epoch of data onto each of the stats in self
        LIMITATIONS FOR NOW: EPOCH MUST BE DONE FIRST
        """
        #for now, just ra
        self.ra = np.concatenate((self.ra, source.ra),0)
        self.dec = np.concatenate((self.dec, source.dec),0)
        self.mag = np.concatenate((self.mag, source.mag),0)
        self.magerr = np.concatenate((self.magerr, source.magerr),0)
        self.shp = np.concatenate((self.shp, source.shp),0)
        self.rnd = np.concatenate((self.rnd, source.rnd),0)

    def append_Band(self, source):
        """
        INPUT: instance of Source()
        FUNC: This will append an band of data onto each of the stats in self
        LIMITATIONS FOR NOW: EPOCH MUST BE DONE FIRST
        """
        self.ra = np.stack((self.ra, source.ra)) 
        self.dec = np.stack((self.dec, source.dec)) 
        self.mag = np.stack((self.mag, source.mag)) 
        self.magerr = np.stack((self.magerr, source.magerr)) 
        self.shp = np.stack((self.shp, source.shp)) 
        self.rnd = np.stack((self.rnd, source.rnd)) 

        
    def set_nan(self):
        pass

    def __repr__(self):
        return("%.4f %.4f %.4f"%(self.RA, self.DEC, self.MAG))

if __name__=='__main__':
    s1 = Source( 1, 2, 3, 4, 5, 6)
    s2 = Source( 2, 3, 4, 5, 6, 7)
    s1.append_Epoch(s2)
    s1.append_Band(s1)
    s1.set_means()
    print(s1.ra)
    print(s1)
