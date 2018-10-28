
def load(config='../config'):
    with open(config,'r') as c:
        return(c.readlines())


def get_value(string, lines=False):
    if not lines: lines = load()
    for i in range(len(lines)):
        if string in lines[i]:
            return(parse(lines[i]))



def parse(line):
    return(line[ line.index('=')+1: ])

if __name__=='__main__':
    print('Band1 Zero Point', float(get_value('Band1 Zero Point')))
