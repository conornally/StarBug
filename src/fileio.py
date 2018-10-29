from fitsclass import FITS

def listfromtxt(filename, comment='#'):
    """INPUT: filename of txt file containing list of file names
    NB: needs to be generalised to automatically or manually take out '\n' chars
    """
    with open(filename, 'r') as f:
        names=[]
        for line in f.readlines():
            if line[0] != comment: 
                names.append(line[:-1])
    return names

def fitsfromtxt(filename):
    return([ FITS( name) for name in listfromtxt(filename)])


if __name__=='__main__':
    print(listfromtxt('test/files.txt'))
    print(fitsfromtxt('test/files.txt'))
