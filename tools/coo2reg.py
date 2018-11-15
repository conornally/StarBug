"""
Converts *.coo output from daophot > find to a region file viewable by ds9
Usage:~$ python coo2reg.py file.coo [circle radius]
"""
from os import system
from sys import argv
from numpy import genfromtxt

if len(argv) <= 1: quit('Usage:~$ python coo2reg.py file.coo [circle radius]')
coo = genfromtxt(argv[1], skip_header=2)[:,1:3]
r = argv[2] if len(argv)==3 else 5
with open(argv[1][:-4]+str('.reg'), 'w') as reg:
    for line in coo:
        reg.write('circle(%f, %f, %f)\n'%(line[0], line[1], r))
system('wc %s'%argv[1][:-4]+str('.reg')) 
