import os

from .makePfsConfig import CatalogId, Lamps
from .spectrum import ArcSpectrum, FlatSpectrum, SlopeSpectrum, CombSpectrum, TextSpectrum, ConstantSpectrum


class SpectrumLibrary:
    def __init__(self, detector, skyModel, skyVarianceOnly=False):
        self.detector = detector
        self.skyModel = skyModel
        self.skyVarianceOnly = skyVarianceOnly

    def getSpectrum(self, catId, objId):
        return {
            CatalogId.ARC: self.getArcSpectrum,
            CatalogId.QUARTZ: self.getQuartzSpectrum,
            CatalogId.FLUXSTD: self.getFluxStdSpectrum,
            CatalogId.SKY: self.getSkySpectrum,
            CatalogId.NULL: self.getNullSpectrum,
            CatalogId.COMB: self.getCombSpectrum,
            CatalogId.SCIENCE: self.getScienceSpectrum,
        }[catId](objId)

    def getArcSpectrum(self, objId):
        lamps = []
        if objId & Lamps.NE:
            lamps.append("NeI")
        if objId & Lamps.HG:
            lamps.append("HgI")
        if objId & Lamps.XE:
            lamps.append("XeI")
        if objId & Lamps.CD:
            lamps.append("CdI")
        if objId & Lamps.KR:
            lamps.append("KrI")
        return ArcSpectrum(lamps)

    def getQuartzSpectrum(self, objId):
        assert objId == 0, "Only one type of quartz spectrum"
        return FlatSpectrum(self.detector, gain=100.0)

    def getNullSpectrum(self, objId):
        return ConstantSpectrum(float(objId))

    def getFluxStdSpectrum(self, objId):
        assert objId == 0, "Currently only one type of fluxStd spectrum"
        return SlopeSpectrum(self.detector, gain=20.0)

    def getSkySpectrum(self, objId):
        assert objId == 0, "Currently only one type of sky spectrum"
        return self.skyModel.getSkyAt(varianceOnly=self.skyVarianceOnly)

    def getCombSpectrum(self, objId):
        return CombSpectrum(spacing=float(objId), gain=200000.0)

    def getScienceSpectrum(self, objId):
        # XXX ignoring objId for now
        filename = os.path.join(os.environ.get("DRP_INSTDATA_DIR", "."),
                                "data", "objects", "ex_gal_sb.dat")
        return TextSpectrum(filename)
