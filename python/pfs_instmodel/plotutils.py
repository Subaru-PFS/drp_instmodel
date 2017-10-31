from __future__ import division

import numpy

def asinh(inputArray, scaleMin=None, scaleMax=None, nonLinear=None, doNorm=True):
    """Performs asinh scaling of the input numpy array.

    @type inputArray: numpy array
    @param inputArray: image data array
    @type scaleMin: number or None
    @param scaleMin: minimum data value
    @type scaleMax: number or None
    @param scaleMax: maximum data value
    @type nonLinear: number or None
    @param nonLinear: non-linearity factor
    @rtype: numpy array
    @return: image data array

    """
    
    # If we are normalizing to 0..1, make sure we have floats. Else leave pixel type alone.
    if doNorm:
        imageData = numpy.array(inputArray, copy=True, dtype='f4')
    else:
        imageData = numpy.array(inputArray, copy=True)
    
    if scaleMin is None:
        scaleMin = imageData.min()
    if scaleMax is None:
        scaleMax = imageData.max()
    if nonLinear is None:
        nonLinear = 0.1*(scaleMax - scaleMin) + scaleMin

    if doNorm:
        factor = numpy.arcsinh((scaleMax - scaleMin)/nonLinear)
    else:
        factor = 1
    indices0 = imageData < scaleMin
    indices2 = imageData > scaleMax
    indices1 = ~(indices0 | indices2)
    imageData[indices0] = 0.0
    imageData[indices2] = 1.0
    imageData[indices1] = numpy.arcsinh((imageData[indices1] - scaleMin)/nonLinear)/factor

    return imageData
