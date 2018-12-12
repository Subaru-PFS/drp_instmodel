from pfs_instmodel.slit import Slit
from pfs_instmodel.utils.generators import keepOdd, keepEven
from pfs_instmodel.makePfsConfig import (makeArcConfig, makeFlatConfig, makeCombConfig, makeScienceConfig,
                                         makeConstantConfig, Lamps)

# Fiber identifiers for different slit configurations
singleFiber = [315]
doubleFiber = [311, 315]
lamFibers = [2, 65, 191, 254, 315, 337, 400, 463, 589, 650]
allFibers = Slit(1).scienceFibers
oddFibers = [ii for ii in keepOdd(allFibers)]
evenFibers = [ii for ii in keepEven(allFibers)]

# Configurations with a single fiber
oneArc = makeArcConfig(0, 0, Lamps.NE | Lamps.HG | Lamps.XE, singleFiber)
oneNe = makeArcConfig(0, 0, Lamps.NE, singleFiber)
oneHg = makeArcConfig(0, 0, Lamps.HG, singleFiber)
oneFlat = makeFlatConfig(0, 0, singleFiber)
oneComb = makeCombConfig(0, 0, singleFiber, 5)
oneConst = makeConstantConfig(0, 0, singleFiber, 1.0e4)
oneSky = makeScienceConfig(0, 0, singleFiber, fracSky=1.0, fracFluxStd=0)
oneObj = makeScienceConfig(0, 0, singleFiber, fracSky=0.0, fracFluxStd=0.0,
                           minScienceMag=18.0, maxScienceMag=18.0)

twoArc = makeArcConfig(0, 0, Lamps.NE | Lamps.HG | Lamps.XE, doubleFiber)
twoNe = makeArcConfig(0, 0, Lamps.NE, doubleFiber)
twoHg = makeArcConfig(0, 0, Lamps.HG, doubleFiber)
twoFlat = makeFlatConfig(0, 0, doubleFiber)
twoComb = makeCombConfig(0, 0, doubleFiber, 5)
twoConst = makeConstantConfig(0, 0, doubleFiber, 1.0e4)
twoSky = makeScienceConfig(0, 0, doubleFiber, fracSky=1.0, fracFluxStd=0)
twoObj = makeScienceConfig(0, 0, doubleFiber, fracSky=0.5, fracFluxStd=0.0,
                           minScienceMag=18.0, maxScienceMag=18.0)

# Configurations with the LAM slit
lamArc = makeArcConfig(0, 0, Lamps.NE | Lamps.HG | Lamps.XE, lamFibers)
lamNe = makeArcConfig(0, 0, Lamps.NE, lamFibers)
lamHg = makeArcConfig(0, 0, Lamps.HG, lamFibers)
lamFlat = makeFlatConfig(0, 0, lamFibers)
lamComb = makeCombConfig(0, 0, lamFibers, 5)
lamConst = makeConstantConfig(0, 0, lamFibers, 1.0e4)
lamSky = makeScienceConfig(0, 0, lamFibers, fracSky=1.0, fracFluxStd=0)
lamObj = makeScienceConfig(0, 0, lamFibers, fracSky=0.1, fracFluxStd=0.1,
                           minScienceMag=18.0, maxScienceMag=22.0)

# Configurations with all fibers
allArc = makeArcConfig(0, 0, Lamps.NE | Lamps.HG | Lamps.XE, allFibers)
allNe = makeArcConfig(0, 0, Lamps.NE, allFibers)
allHg = makeArcConfig(0, 0, Lamps.HG, allFibers)
allFlat = makeFlatConfig(0, 0, allFibers)
allComb = makeCombConfig(0, 0, allFibers, 5)
allConst = makeConstantConfig(0, 0, allFibers, 1.0e4)
allSky = makeScienceConfig(0, 0, allFibers, fracSky=1.0, fracFluxStd=0)
allObj = makeScienceConfig(0, 0, allFibers)

# Configurations with odd fibers only
oddArc = makeArcConfig(0, 0, Lamps.NE | Lamps.HG | Lamps.XE, oddFibers)
oddNe = makeArcConfig(0, 0, Lamps.NE, oddFibers)
oddHg = makeArcConfig(0, 0, Lamps.HG, oddFibers)
oddFlat = makeFlatConfig(0, 0, oddFibers)
oddComb = makeCombConfig(0, 0, oddFibers, 5)
oddConst = makeConstantConfig(0, 0, oddFibers, 1.0e4)
oddSky = makeScienceConfig(0, 0, oddFibers, fracSky=1.0, fracFluxStd=0)
oddObj = makeScienceConfig(0, 0, oddFibers)

# Configurations with even fibers only
evenArc = makeArcConfig(0, 0, Lamps.NE | Lamps.HG | Lamps.XE, evenFibers)
evenNe = makeArcConfig(0, 0, Lamps.NE, evenFibers)
evenHg = makeArcConfig(0, 0, Lamps.HG, evenFibers)
evenFlat = makeFlatConfig(0, 0, evenFibers)
evenComb = makeCombConfig(0, 0, evenFibers, 5)
evenConst = makeConstantConfig(0, 0, evenFibers, 1.0e4)
evenSky = makeScienceConfig(0, 0, evenFibers, fracSky=1.0, fracFluxStd=0)
evenObj = makeScienceConfig(0, 0, evenFibers)
