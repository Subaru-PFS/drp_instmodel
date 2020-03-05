import enum


class Arm(enum.Enum):
    """The spectrograph arm"""
    BLUE = "b"
    RED = "r"
    NIR = "n"
    MID = "m"

    @classmethod
    def fromDetectorName(cls, detectorName):
        if (len(detectorName) != 2 or
                detectorName[0] not in {'b', 'r', 'n', 'm'} or
                detectorName[1] not in {'1', '2', '3', '4'}):
            raise ValueError('detectorName (%s) must be in the form "r1"' % (detectorName,))
        return cls(detectorName[0])

    @classmethod
    def fromBand(cls, band):
        menu = {"blue": cls.BLUE,
                "red": cls.RED,
                "nir": cls.NIR,
                "mid": cls.MID,
                }
        return menu[band.lower()]

    def toBand(self):
        cls = self.__class__
        menu = {cls.BLUE: "blue",
                cls.RED: "red",
                cls.NIR: "nir",
                cls.MID: "mid",
                }
        return menu[self]

    def toCcd(self):
        cls = self.__class__
        menu = {cls.BLUE: 0,
                cls.RED: 1,
                cls.NIR: 1,
                cls.MID: 2,
                }
        return menu[self]
