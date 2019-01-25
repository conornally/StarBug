import numpy as np
np.warnings.filterwarnings('ignore')
from src.starbugclass import StarBug
import argparse

parser = argparse.ArgumentParser(description='parse some args')
parser.add_argument('-test', help='here is some text', type=str)
parser.add_argument('-cat', help='Load in catalog', type=str)
args = parser.parse_args()

starbug = StarBug()
if args.cat:
    starbug._loadcatAUTO('autocat', args.cat, '', 'sextractor') #change to custom

starbug.mainloop()
