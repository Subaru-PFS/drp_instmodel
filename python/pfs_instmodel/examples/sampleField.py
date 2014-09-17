from pfs_instmodel.schema.probes import PROBE
    
fiberLim = 313
bundleSpacing = 36

def fullField():
    global fiberLim
    return range(-fiberLim, fiberLim+1)

def keepInBundle(input=None):
    global bundleSpacing
    if input is None:
        input = fullField()
    for i in input:
        if i%bundleSpacing != 0:
            yield i

def keepOdd(input):
    for i in input:
        if i%2 == 1:
            yield i

bundledField = [i for i in keepInBundle(fullField())]
centerRange = [i for i in keepInBundle(range(-5,6))]
edgeRange = [i for i in keepInBundle(range(-fiberLim,-fiberLim+10))]
centerAndEdge = centerRange + edgeRange

# Minimal sanity test field.
quickfield = ([PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(-300,-298)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(-298,-21)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(-21,-19)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(-19,-11)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(-11,-9)] +

              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(-11,-2)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(-2,3)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(3,10)] +

              [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(10,12)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(12,20)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(20,22)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(22,299)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(299,301)])

combField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in bundledField]
flatField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in bundledField]
skyField =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in bundledField]

sparseCombField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMCOMB' if i%100 == 0 or i in (-fiberLim, fiberLim+1) else 'UNPLUGGED')) 
                   for i in fullField()]
sparseFlatField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMFLAT' if i%50 == 0 else 'UNPLUGGED')) for i in fullField()]

centerComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in centerRange]
centerFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in centerRange]
centerSky  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in centerRange]

edgeFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in edgeRange]
edgeSky  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in edgeRange]

centerAndEdgeFlat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in centerAndEdge]
centerAndEdgeSky = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in centerAndEdge]

centerFlatx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in keepOdd(centerRange)]
centerSkyx2  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in keepOdd(centerRange)]

edgeFlatx2 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in keepOdd(edgeRange)]
edgeSkyx2  = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in keepOdd(edgeRange)]

quickComb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in (-fiberLim, 0, fiberLim+1)]

empty = []
