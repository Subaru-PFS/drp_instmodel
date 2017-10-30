import numpy as np
from pfs_instmodel.schema.probes import PROBE
    
fiberLim = 313

def fullField():
    global fiberLim
    return range(-fiberLim, fiberLim+1)

def halfField():
    global fiberLim
    return range(-fiberLim, 1)

def keepNth(input, N):
    for i in input:
        if i%N == 0:
            yield i

def keepOdd(input):
    for i in input:
        if i%2 == 1:
            yield i

def keepEven(input):
    for i in input:
        if i%2 == 0:
            yield i

class Fiber(object):
    SCIENCE = 1
    ENGINEERING = 2
    BLANK = 3

class Slit(object):
    block1 = ([Fiber.ENGINEERING] +
              42 * [Fiber.SCIENCE] +
              [Fiber.BLANK])
    block2 = ([Fiber.ENGINEERING] +
              45 * [Fiber.SCIENCE] +
              [Fiber.BLANK])
    block3 = ([Fiber.ENGINEERING] +
              [Fiber.BLANK] +
              42 * [Fiber.SCIENCE] +
              [Fiber.BLANK])
    block4 = ([Fiber.ENGINEERING] +
              42 * [Fiber.SCIENCE] +
              [Fiber.ENGINEERING])
    block5 = ([Fiber.ENGINEERING] +
              45 * [Fiber.SCIENCE] +
              [Fiber.ENGINEERING])
    gap = 19*[Fiber.BLANK]

    slit1 = (block1 +
             block2 +
             block3 +
             block2 +
             block3 +
             block1 +
             block4 +
             gap +
             block5 +
             block1[::-1] +
             block3[::-1] +
             block1[::-1] +
             block3[::-1] +
             block2[::-1] +
             block1[::-1])

    # The second slit type is identical to the first but with
    # three blanked-off science fibers.
    slit2 = list(slit1)
    for f in (280, 309, 359):
        slit2[f-1] = Fiber.BLANK
    slit2 = tuple(slit2)

    slits = slit1, slit1, slit2, slit2

    def __init__(self, slitId):
        self.fiber0 = 1
        scienceFibers = []
        engineeringFibers = []

        for i, fiberType in enumerate(self.slits[slitId]):
            if fiberType is Fiber.SCIENCE:
                scienceFibers.append(i)
            elif fiberType is Fiber.ENGINEERING:
                engineeringFibers.append(i)

        self._scienceFibers = np.array(scienceFibers, dtype='i4')

    @property
    def scienceFibers(self):
        return self._scienceFibers + self.fiber0
    
    def scienceFiberToSlitPos(self, scienceFiberNum):
        """ Return the slit position for the given science fiber.
        
        This is only used for the engineering slits as LAM, which were 
        sent defined in terms of the science fibers.
        """

        if scienceFiberNum < 0:
            indexNum = scienceFiberNum + 300
        else:
            indexNum = scienceFiberNum + 299
            
        return self.scienceFibers[indexNum]

slit1 = Slit(1)

# Science fiber numbers, -300 to +300:

# The illuminated fibers on LAM slit #1
LamSlit1 = [2, 65, 191, 254, 315, 337, 400, 463, 589, 650]

Lam1Flat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in LamSlit1]
Lam1Arcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in LamSlit1]
Lam1Comb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in LamSlit1]
Lam1Slope = [PROBE(i,0.0,1.0,100.0,200.0,'SIMSLOPE', ()) for i in LamSlit1]
Lam1CdArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('CdI',)) for i in LamSlit1]
Lam1HgArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('HgI',)) for i in LamSlit1]
Lam1NeArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('NeI',)) for i in LamSlit1]
Lam1XeArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('XeI',)) for i in LamSlit1]
Lam1KrArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('KrI',)) for i in LamSlit1]

bundledField = slit1.scienceFibers

centerRange = [slit1.scienceFiberToSlitPos(i) for i in -1,1]
edgeRange = [slit1.scienceFiberToSlitPos(i) for i in -300,300]
centerAndEdge = centerRange + edgeRange

combField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in slit1.scienceFibers]
flatField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in slit1.scienceFibers]
skyField =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in slit1.scienceFibers]

combFieldx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepEven(bundledField)]
flatFieldx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepEven(bundledField)]
skyFieldx2 =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepEven(bundledField)]

combFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepNth(bundledField, 10)]
arcFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in keepNth(bundledField, 10)]
flatFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepNth(bundledField, 10)]
skyFieldx10 =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepNth(bundledField, 10)]

combFieldx40 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepNth(bundledField, 40)]
arcFieldx40 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in keepNth(bundledField, 40)]
flatFieldx40 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepNth(bundledField, 40)]
skyFieldx40 =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepNth(bundledField, 40)]

centerAndEdgeFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in centerAndEdge]
centerAndEdgeComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in centerAndEdge]
centerAndEdgeSky = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in centerAndEdge]

quickComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in (-fiberLim, -fiberLim/2, 0, 1, fiberLim+1)]
quickFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in (-fiberLim, -fiberLim/2, 0, 1, fiberLim+1)]

oneArc = (PROBE(200,0.0,1.0,100.0,200.0,'SIMARC', ()),)
oneFlat = (PROBE(200,0.0,1.0,100.0,200.0,'SIMFLAT', ()),)
oneSlope = (PROBE(200,0.0,1.0,100.0,200.0,'SIMSLOPE', ()),)
oneSky = (PROBE(200,0.0,1.0,100.0,200.0,'SKY', ()),)

empty = []
