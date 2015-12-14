from pfs_instmodel.schema.probes import PROBE
    
fiberLim = 313
bundleSpacing = 36

def fullField():
    global fiberLim
    return range(-fiberLim, fiberLim+1)

def halfField():
    global fiberLim
    return range(-fiberLim, 1)

def keepInBundle(input=None):
    global bundleSpacing
    if input is None:
        input = fullField()
    for i in input:
        if i == 0 or i%bundleSpacing != 0:
            yield i

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
