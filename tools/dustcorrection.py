import numpy as np
import matplotlib.pyplot as plt


def dustCorrection(data):
    U = 0.3543 #um
    G = 0.4770 #um
    R = 0.6231 #um

    Rv= 3.1

    Cu = a(U) + b(U)/Rv
    Cg = a(G) + b(G)/Rv
    Cr = a(R) + b(R)/Rv

    ug=  data[:,0]
    dug=  data[:,1]
    gr=  data[:,2]
    dgr=  data[:,3]


    mainsequence = np.genfromtxt("catalogs/pleiades_sloane")
    mask = ( mainsequence[:,1]<0.7)
    mainsequence = mainsequence[mask]
    coeffs = np.polyfit(mainsequence[:,1], mainsequence[:,0],6)
    x=np.arange(-0.3,0.75,0.1)
    #y=np.polyval(coeffs, x)
    deriv_coeffs = list(reversed([ i*c for i, c in enumerate(reversed(coeffs))]))
    
    x=np.arange(-0.3,0.75,0.001)
    y=np.polyval(coeffs, x)
    

    Av = np.arange(0,5.0,0.01)
    Chivals = np.zeros(Av.shape)
    for i, av in enumerate(Av):
        chi=0
        tmpUG = ug-((Cu-Cg)*av)
        tmpGR = ug-((Cg-Cr)*av)
        chi = (tmpUG - np.polyval(coeffs, tmpGR))**2.0/( dug**2.0 + (np.polyval(deriv_coeffs, dgr)*dgr)**2.0  )
        print(chi)
        #if(np.isfinite(schi)): chi+= schi
        Chivals[i] = sum(chi)
    print(Chivals)
    minAv = Av[np.argmin(Chivals)]
    print(minAv)
    print("U-G: %f"%(minAv*(Cu-Cg)))
    print("G-R: %f"%(minAv*(Cg-Cr)))
    fig = plt.figure()
    """
    for s in self.sourcelist:
        s.mag[:,0]-=(Cg*minAv)
        s.mag[:,1]-=(Cr*minAv)
        s.mag[:,2]-=(Cu*minAv)
        s.construct_colours([2,1],[1,0])
        plt.scatter(s.colours[1], -s.colours[0],c='k', s=0.1)
    """
    plt.scatter(data[:,0], data[:,2])
    plt.plot((Cg-Cr)*Av,(Cu-Cg)*Av)
    plt.plot(x,y)
    plt.gca().invert_yaxis()
    plt.plot(Av, np.log(Chivals))
    plt.show()

    
    #plt.show()


def a(x):
    y = x-1.82
    #coeffs = [1, 0.17699, -0.5044, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999]
    coeffs = [0.32999, -0.77530, 0.01979, 0.72085, -0.02427, -0.50447, 0.17699, 1.0]
    return np.polyval(coeffs, y)

def b(x):
    y = x-1.82
    coeffs = [-2.09002, 5.30260, -0.62251, -5.38434, 1.07233, 2.28305, 1.41338, 0]
    return np.polyval(coeffs, y)

if __name__=='__main__':
    data = np.genfromtxt("test/dereddening_teststars.txt", skip_header=1)
    dustCorrection(data)

