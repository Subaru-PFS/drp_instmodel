from functools import partial
from pfs.instmodel.slit import Slit
from pfs.instmodel.utils.generators import keepOdd, keepEven
from pfs.instmodel.makePfsConfig import (makeArcDesign, makeFlatDesign, makeCombDesign, makeScienceDesign,
                                         makeDarkDesign, makeConstantDesign, Lamps)

__all__ = ["dark",
           "oneArc", "oneNe", "oneHg", "oneFlat", "oneComb", "oneConst", "oneSky", "oneObj",
           "twoArc", "twoNe", "twoHg", "twoFlat", "twoComb", "twoConst", "twoSky", "twoObj",
           "lamArc", "lamNe", "lamHg", "lamFlat", "lamComb", "lamConst", "lamSky", "lamObj",
           "allArc", "allNe", "allHg", "allFlat", "allComb", "allConst", "allSky", "allObj",
           "oddArc", "oddNe", "oddHg", "oddFlat", "oddComb", "oddConst", "oddSky", "oddObj",
           "evenArc", "evenNe", "evenHg", "evenFlat", "evenComb", "evenConst", "evenSky", "evenObj",
           ]

# Fiber identifiers for different slit configurations
singleFiber = [315]
doubleFiber = [311, 315]
lamFibers = [2, 65, 191, 254, 315, 337, 400, 463, 589, 650]
allFibers = Slit(1).scienceFibers
oddFibers = [ii for ii in keepOdd(allFibers)]
evenFibers = [ii for ii in keepEven(allFibers)]

dark = partial(makeDarkDesign, 0, allFibers)

# Configurations with a single fiber
oneArc = partial(makeArcDesign, 0, Lamps.NE | Lamps.HG | Lamps.XE, singleFiber)
oneNe = partial(partial, makeArcDesign, 0, Lamps.NE, singleFiber)
oneHg = partial(makeArcDesign, 0, Lamps.HG, singleFiber)
oneFlat = partial(makeFlatDesign, 0, singleFiber)
oneComb = partial(makeCombDesign, 0, singleFiber, 5)
oneConst = partial(makeConstantDesign, 0, singleFiber, 1.0e4)
oneSky = partial(makeScienceDesign, 0, singleFiber, fracSky=1.0, fracFluxStd=0)
oneObj = partial(makeScienceDesign, 0, singleFiber, fracSky=0.0, fracFluxStd=0.0,
                 minScienceMag=18.0, maxScienceMag=18.0)

# Configurations with two fibers
twoArc = partial(makeArcDesign, 0, Lamps.NE | Lamps.HG | Lamps.XE, doubleFiber)
twoNe = partial(makeArcDesign, 0, Lamps.NE, doubleFiber)
twoHg = partial(makeArcDesign, 0, Lamps.HG, doubleFiber)
twoFlat = partial(makeFlatDesign, 0, doubleFiber)
twoComb = partial(makeCombDesign, 0, doubleFiber, 5)
twoConst = partial(makeConstantDesign, 0, doubleFiber, 1.0e4)
twoSky = partial(makeScienceDesign, 0, doubleFiber, fracSky=1.0, fracFluxStd=0)
twoObj = partial(makeScienceDesign, 0, doubleFiber, fracSky=0.5, fracFluxStd=0.0,
                 minScienceMag=18.0, maxScienceMag=18.0)

# Configurations with the LAM slit
lamArc = partial(makeArcDesign, 0, Lamps.NE | Lamps.HG | Lamps.XE, lamFibers)
lamNe = partial(makeArcDesign, 0, Lamps.NE, lamFibers)
lamHg = partial(makeArcDesign, 0, Lamps.HG, lamFibers)
lamFlat = partial(makeFlatDesign, 0, lamFibers)
lamComb = partial(makeCombDesign, 0, lamFibers, 5)
lamConst = partial(makeConstantDesign, 0, lamFibers, 1.0e4)
lamSky = partial(makeScienceDesign, 0, lamFibers, fracSky=1.0, fracFluxStd=0)
lamObj = partial(makeScienceDesign, 0, lamFibers, fracSky=0.1, fracFluxStd=0.1,
                 minScienceMag=18.0, maxScienceMag=22.0)

# Configurations with all fibers
allArc = partial(makeArcDesign, 0, Lamps.NE | Lamps.HG | Lamps.XE, allFibers)
allNe = partial(makeArcDesign, 0, Lamps.NE, allFibers)
allHg = partial(makeArcDesign, 0, Lamps.HG, allFibers)
allFlat = partial(makeFlatDesign, 0, allFibers)
allComb = partial(makeCombDesign, 0, allFibers, 5)
allConst = partial(makeConstantDesign, 0, allFibers, 1.0e4)
allSky = partial(makeScienceDesign, 0, allFibers, fracSky=1.0, fracFluxStd=0)
allObj = partial(makeScienceDesign, 0, allFibers)

# Configurations with odd fibers only
oddArc = partial(makeArcDesign, 0, Lamps.NE | Lamps.HG | Lamps.XE, oddFibers)
oddNe = partial(makeArcDesign, 0, Lamps.NE, oddFibers)
oddHg = partial(makeArcDesign, 0, Lamps.HG, oddFibers)
oddFlat = partial(makeFlatDesign, 0, oddFibers)
oddComb = partial(makeCombDesign, 0, oddFibers, 5)
oddConst = partial(makeConstantDesign, 0, oddFibers, 1.0e4)
oddSky = partial(makeScienceDesign, 0, oddFibers, fracSky=1.0, fracFluxStd=0)
oddObj = partial(makeScienceDesign, 0, oddFibers)

# Configurations with even fibers only
evenArc = partial(makeArcDesign, 0, Lamps.NE | Lamps.HG | Lamps.XE, evenFibers)
evenNe = partial(makeArcDesign, 0, Lamps.NE, evenFibers)
evenHg = partial(makeArcDesign, 0, Lamps.HG, evenFibers)
evenFlat = partial(makeFlatDesign, 0, evenFibers)
evenComb = partial(makeCombDesign, 0, evenFibers, 5)
evenConst = partial(makeConstantDesign, 0, evenFibers, 1.0e4)
evenSky = partial(makeScienceDesign, 0, evenFibers, fracSky=1.0, fracFluxStd=0)
evenObj = partial(makeScienceDesign, 0, evenFibers)
