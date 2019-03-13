import numpy as np

def rebin(a, *args):
    """ integer factor rebin(a, *new_axis_sizes), taken from scipy cookbook. """

    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
        ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
        [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)]

    binArray = eval(''.join(evList))
    return binArray
