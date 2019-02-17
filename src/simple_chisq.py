import sys, numpy as np
from scipy import optimize


def readfromFile(filename):
    return np.genfromtxt(filename)

def chisqd(coeffs, x,y,yerr, order):
    model = polynomial(coeffs, x, order)
    return np.sum(((model - y)/yerr)**2.)

def polynomial(coeffs, x, order):
    y = 0
    for i in range(order+1):
        y+= coeffs[i]*x**i
    return y


def getCHI(data,order):

    x = data[:,0]
    y = data[:,1]
    dx = data[:,2]
    dy = data[:,3]
    initial=[1]*(order+1)
    return optimize.minimize(chisqd, initial, args=(x, y, dy, order))
    



if __name__=='__main__':
    if(len(sys.argv)>2): order = int(sys.argv[2])
    else: order =1
    out = main(readfromFile(sys.argv[1]), order)
    print("chi: %f"%out['fun'])
    sys.stdout.write("y= ")
    for i in range(len(out['x'])):
        sys.stdout.write("%f*x^%d + "%(out['x'][i],i))
    print('\n')

