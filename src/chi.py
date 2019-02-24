import sys, numpy as np
from scipy import optimize

def chisqd(coeffs, x, y, yerr):
    model = np.polyval(coeffs, x)
    return np.sum(((model - y)/yerr)**2.)

def get_CHI(x, y, dy, order):
    initial=np.ones((order+1))
    out= optimize.minimize(chisqd, initial, args=(x, y, dy))
    return(out['fun'], out['x'])


def get_simpleChierr(coeffs, minX):
    factor=0.01

def get_CHIerr( x, y, dy, coeffs):
    factor=0.01
    sigma = np.sqrt( abs( 2*(len(x)-len(coeffs))))
    chimin=chisqd(coeffs, x, y, dy)
    bestfit = coeffs

    param_errors = np.zeros((len(coeffs)))
    for i, param in enumerate(coeffs):
        chi = chimin
        coeffs = np.copy(bestfit)
        while( chi < (chimin+sigma)):
            coeffs[i]+=factor
            chi = chisqd(coeffs, x, y, dy)
        param_errors[i] = abs(coeffs[i] - bestfit[i])


    #note keeps the final coeff=0, if put through polyval, will not be correct
    deriv_bestfit = list(reversed([i*c for i,c in enumerate(list(reversed(bestfit)))]))
    error = np.sqrt(sum([(bestfit[i]*param_errors[i])**2. for i in range(len(param_errors))]))
    return param_errors,bestfit, error





if __name__=='__main__':
    x = [1, 2, 3, 4, 5]
    y = [1.1, 1.9, 3.01, 4.3, 4.9]
    dy= [0.1, 0.1, 0.05, 0.4, 0.1]
    order=3
    chi, coeffs = get_CHI(x,y,dy,order)
    print(x,coeffs)
    get_CHIerr(x, y, dy, coeffs)

    import matplotlib.pyplot as plt
    plt.scatter(x,y)
    x = np.arange(0,5, 0.1)
    plt.plot(x, np.polyval(coeffs, x))
    plt.show()
