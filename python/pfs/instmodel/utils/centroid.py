import numpy as np

def centroid(im):
    """ 1st moment centroids. """

    xFlux = im.sum(axis=0)
    yFlux = im.sum(axis=1)
    flux = yFlux.sum()

    xc1 = (xFlux * np.arange(1, im.shape[1]+1)).sum()
    yc1 = (yFlux * np.arange(1, im.shape[0]+1)).sum()

    return np.array([xc1/flux - 1, yc1/flux - 1])
