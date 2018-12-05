import numpy as np
from pfs_instmodel.schema.probes import PROBE


Sm1Slit = Slit(1).scienceFibers

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
Lam1Sky = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in LamSlit1]

Sm1Flat = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in Sm1Slit]
Sm1Arcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in Sm1Slit]
Sm1Comb = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in Sm1Slit]
Sm1Slope = [PROBE(i,0.0,1.0,100.0,200.0,'SIMSLOPE', ()) for i in Sm1Slit]
Sm1CdArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('CdI',)) for i in Sm1Slit]
Sm1HgArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('HgI',)) for i in Sm1Slit]
Sm1NeArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('NeI',)) for i in Sm1Slit]
Sm1XeArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('XeI',)) for i in Sm1Slit]
Sm1KrArcs = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('KrI',)) for i in Sm1Slit]
Sm1Sky = [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in Sm1Slit]

bundledField = Sm1Slit

combFieldEven = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepEven(bundledField)]
arcFieldEven = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in keepEven(bundledField)]
neArcsEven = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('NeI',)) for i in keepEven(bundledField)]
flatFieldEven = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepEven(bundledField)]
skyFieldEven =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepEven(bundledField)]

combFieldOdd = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepOdd(bundledField)]
arcFieldOdd = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in keepOdd(bundledField)]
neArcsOdd = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ('NeI',)) for i in keepOdd(bundledField)]
flatFieldOdd = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepOdd(bundledField)]
skyFieldOdd =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepOdd(bundledField)]

combFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepNth(bundledField, 10)]
arcFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in keepNth(bundledField, 10)]
flatFieldx10 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepNth(bundledField, 10)]
skyFieldx10 =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepNth(bundledField, 10)]

combFieldx40 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB', ()) for i in keepNth(bundledField, 40)]
arcFieldx40 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMARC', ()) for i in keepNth(bundledField, 40)]
flatFieldx40 = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT', ()) for i in keepNth(bundledField, 40)]
skyFieldx40 =  [PROBE(i,0.0,1.0,100.0,200.0,'SKY', ()) for i in keepNth(bundledField, 40)]

oneArc = (PROBE(315,0.0,1.0,100.0,200.0,'SIMARC', ()),)
oneNe = (PROBE(315,0.0,1.0,100.0,200.0,'SIMARC', ('NeI',)),)
oneHg = (PROBE(315,0.0,1.0,100.0,200.0,'SIMARC', ('HgI',)),)
oneFlat = (PROBE(315,0.0,1.0,100.0,200.0,'SIMFLAT', ()),)
oneSlope = (PROBE(315,0.0,1.0,100.0,200.0,'SIMSLOPE', ()),)
oneSky = (PROBE(315,0.0,1.0,100.0,200.0,'SKY', ()),)

empty = []
