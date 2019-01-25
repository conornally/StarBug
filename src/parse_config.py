
def load(config='config'):
    with open(config,'r') as c:
        options = {}
        for line in c.readlines():
            if line != '\n' and line[0] != '#':
                if '\n' in line: line = parse_key(line, '\n')
                options[parse_key(line)] = parse_value(line)

        return(options)


def get_value(string, lines=False, dtype=str, configfile='config'):
    """
    if not lines: lines = load(configfile)
    print(lines)
    for key in lines.keys():
        if string in key:
            return(dtype(parse_value(lines[i])))
    """
    if not lines: lines = load(configfile)
    if string in lines.keys():
        return(dtype(lines[string]))
    else: return False


def parse_value(line, separator='='):
    #if '#' in line: line = line[:line.index('#')]
    if separator not in line: return line
    try: return(float(line[ line.index(separator)+1: ]))
    except: return(str(line[ line.index(separator)+1: ]))

def parse_key(line, separator='='):
    if separator not in line: return line
    key = line[:line.index(separator)]
    if key[-1] == ' ': key = key[:-1]
    return(key)


def display(options):
    if type(options) == dict:
        keys = options.keys()
        print(                 '\x1b[4m##PHOTOMETRY##\x1b[0m')
        if 'FWHM' in keys:      print('FWHM:      %s'%options['FWHM'])
        if 'Threshold' in keys: print('Threshold: %s'%options['Threshold'])
        if 'Sigma' in keys:     print('Sigma:     %s'%options['Sigma'])
        if 'SharpHigh' in keys: print('SharpHigh: %s'%options['SharpHigh'])
        if 'SharpLow' in keys:  print('SharpLow:  %s'%options['SharpLow'])
        if 'RoundHigh' in keys: print('RoundHigh: %s'%options['RoundHigh'])
        if 'RoundLow' in keys:  print('RoundLow:  %s'%options['RoundLow'])


if __name__=='__main__':
    print('Band1 Zero Point', float(get_value('Band1 Zero Point')))
