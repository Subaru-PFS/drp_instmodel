import enum
import numpy as np

NUM_FIBERS = 651  # Number of fibers per spectrograph

@enum.unique
class Fiber(enum.Enum):
    SCIENCE = 1
    ENGINEERING = 2
    BLANK = 3


class Slit(object):
    """A definition of a slit.

    This structure clarifies some of logic of the slit layout.

    Parameters
    ----------
    spectrograph : `int`
        Spectrograph number, 1..4.
    allHoles : `bool`, optional
        Return all holes if ``True``, otherwise only the filled holes for that
        spectrograph.
    """

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

    slits = slit1, slit2, slit1, slit2

    def __init__(self, spectrograph, allHoles=False):
        slitId = spectrograph - 1
        self.fiber0 = NUM_FIBERS*slitId + 1
        scienceFibers = []
        engineeringFibers = []

        for i, fiberType in enumerate(self.slit1 if allHoles else self.slits[slitId]):
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

