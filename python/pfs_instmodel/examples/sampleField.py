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

    slits = slit1, slit1, slit1, slit1
    
    def __init__(self, slitId):
        self.fiber0 = -325
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
        """ Return the slit position for the given science fiber """

        if scienceFiberNum < 0:
            indexNum = scienceFiberNum + 300
        else:
            indexNum = scienceFiberNum + 299
            
        return self.scienceFibers[indexNum]

slit1 = Slit(1)

# Science fiber numbers, -300 to +300:
#LamSlit1Fibers = [-300, -240, -100, -30, -1, 1, 30, 100, 170, 240, 300]
LamSlit1Fibers = [-300, -240, -239, -100, -30, -29, -1, 1, 30, 100, 170, 240, 300]
LamSlit2Fibers = [-299, -290, -246, -194, -148, -112, -56, -10 -2,
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

bundledField = slit1.scienceFibers
centerRange = [slit1.scienceFiberToSlitPos(i) for i in range(-3,2)]
edgeRange = [i for i in slit1.scienceFibers[:3]]
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

oneArc = (PROBE(0,0.0,1.0,100.0,200.0,'SIMARC', ()),)
oneFlat = (PROBE(0,0.0,1.0,100.0,200.0,'SIMFLAT', ()),)
oneSky = (PROBE(0,0.0,1.0,100.0,200.0,'SKY', ()),)

empty = []
