import numpy as np
np.warnings.filterwarnings('ignore')
from src.starbug import StarBug
import argparse

parser = argparse.ArgumentParser(description='parse some args')
parser.add_argument('-test', help='here is some text', type=str)
parser.add_argument('-cat', action='append', dest='CATALOGS', default=[], help='Load in catalog', type=str)
parser.add_argument('-fits', action='append', dest='FITSFILES',  default=[], help='File or glob parsable file root in format "FILE*.fits"', type=str)

args = parser.parse_args()

starbug = StarBug()
if args.CATALOGS:
    for i, cat in enumerate(args.CATALOGS):
        starbug._loadcatAUTO('autocat%d'%i, cat, '', 'custom') #change to custom

if args.FITSFILES: print(args.FITSFILES)

starbug.mainloop()
