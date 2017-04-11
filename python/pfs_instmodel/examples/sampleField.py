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
LamSlit1Fibers = [-300, -240, -239, -100, -30, -29, -1, 1, 30, 100, 170, 240, 300]
LamSlit2Fibers = [-299, -290, -246, -194, -148, -112, -56, -10, -2,
                  2, 10, 56, 112, 148, 194, 148, 194, 246, 290]
LamSlit1 = [slit1.scienceFiberToSlitPos(f) for f in LamSlit1Fibers]
LamSlit2 = [slit1.scienceFiberToSlitPos(f) for f in LamSlit2Fibers]

Lam1Flat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in LamSlit1]
Lam1Arcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in LamSlit1]
Lam1CdArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('CdI',)) for i in LamSlit1]
Lam1HgArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('HgI',)) for i in LamSlit1]
Lam1NeArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('NeI',)) for i in LamSlit1]
Lam1XeArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('XeI',)) for i in LamSlit1]
Lam1KrArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('KrI',)) for i in LamSlit1]

bundledField = [i for i in keepInBundle(halfField())]
centerRange = [i for i in keepInBundle(range(-3,2))]
edgeRange = [i for i in keepInBundle(range(-fiberLim, -fiberLim+4))]
centerAndEdge = [-fiberLim, 0] # centerRange + edgeRange
sampledRange = centerRange + edgeRange + [i for i in keepInBundle(range(-fiberLim/2, -fiberLim/2+4))]

combField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in bundledField]
flatField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in bundledField]
skyField =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in bundledField]

combFieldx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepEven(bundledField)]
flatFieldx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepEven(bundledField)]
skyFieldx2 =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepEven(bundledField)]

combFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepNth(bundledField, 10)]
flatFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepNth(bundledField, 10)]
skyFieldx10 =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepNth(bundledField, 10)]

sparseCombField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMCOMB' if i%100 == 0 or i in (-fiberLim, fiberLim+1, ()) else 'UNPLUGGED'), ()) 
                   for i in fullField()]
sparseFlatField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMFLAT' if i%50 == 0 else 'UNPLUGGED'), ()) for i in fullField()]

centerComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in centerRange]
centerFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in centerRange]
centerSky  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in centerRange]

edgeFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in edgeRange]
edgeSky  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in edgeRange]

centerAndEdgeFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in centerAndEdge]
centerAndEdgeComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in centerAndEdge]
centerAndEdgeSky = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in centerAndEdge]

centerFlatx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepEven(centerRange)]
centerSkyx2  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepEven(centerRange)]

edgeFlatx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepEven(edgeRange)]
edgeSkyx2  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepEven(edgeRange)]

quickComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in (-fiberLim, -fiberLim/2, 0, 1, fiberLim+1)]
quickFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in (-fiberLim, -fiberLim/2, 0, 1, fiberLim+1)]

oneArc = (PROBE(0,0.0,1.0,100.0,200.0,'SIMARC', ()),)
oneFlat = (PROBE(0,0.0,1.0,100.0,200.0,'SIMFLAT', ()),)
oneSky = (PROBE(0,0.0,1.0,100.0,200.0,'SKY', ()),)

sampledComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in sampledRange]
sampledFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in sampledRange]
sampledSky  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in sampledRange]
sampledCombx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepEven(sampledRange)]
sampledFlatx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepEven(sampledRange)]
sampledSkyx2  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepEven(sampledRange)]

minFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in (-fiberLim, -fiberLim/2, 0)]
minSky = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in (-fiberLim, -fiberLim/2, 0)]
minComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in (-fiberLim, -fiberLim/2, 0)]

empty = []
