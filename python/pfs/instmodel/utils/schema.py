""" This should either be sugar around pexConfig or around numpy (rec)arrays. Will do for now. """

from . import configFile

class Schema(object):
    """ Encapsulate (an array of) typed variables. Let them be read from and written to various useful forms.

    Notes
    -----

    We want some sugar here. Specifically, for starters:
      - in/out - a one-line form for human convenience ('par file') 
      - in/out - FITS binary table
      - python subclassing: our variables should be accessible as instance variables.
      - in/out - as numpy arrays.

    We need to choose whether the native form is basically an object, or a row in a table. I'll start
    with the object form, with a mechanism for numpy import/export. 

    Types are often the problem. I'll start with native python or numpy types, which have trivial FITS counterparts.

    The LSST pexConfig world is interesting, particularily the types. We should perhaps extend that instead of 
    exercising extreme NIH.

    """

    pass

def makeParClass(*fields):
    """ Create a very slightly type-checked 'type' class which provides a sugary constructor.
    For example, create a set of points:

    Point = makeParClass(('x',int),('y',int),('label',str))
    points = (Point(3,5,'here'),
              Point(2,2,'there'),
              Point(1,2,'nowhere'))
    """
    schemaNames = [f[0] for f in fields]
    schemaTypes = [f[1] for f in fields]
    # slotNames = schemaNames + ['_schema','__init__','__str__']
    slotNames = schemaNames

    # parClass = type('TESTNAME', ('object',), dict(__slots__=slotNames))

    class _scratchClass(object):
        __slots__ = slotNames
        _schema = dict(names=schemaNames, types=schemaTypes)

        def __init__(self, *parFields):
            """ If we don't use pexConfig, flesh out the checking. """ 

            if len(parFields) != len(self._schema['names']):
                raise RuntimeError("wrong number of args (%s vs. %s)" % 
                                   (parFields, self._schema['names']))

            for i in range(len(self._schema['names'])):
                rawVal = parFields[i]
                slotName = self._schema['names'][i]
                slotType = self._schema['types'][i]
                try:
                    # Handle sets and other things with __contains__
                    ok = rawVal in slotType
                    if not ok:
                        raise ValueError("%s not in %s" % (rawVal, slotType))
                    else:
                        val = rawVal
                except TypeError as e:
                    # Or assume the type can convert or blow up appropriately
                    val = slotType(rawVal)

                setattr(self, slotName, val)

                # And should disable writing?

        def __str__(self):
            return ", ".join([("%s=%s" % (name,
                                          getattr(self,name))) for name in self._schema['names']])

    return _scratchClass

def loadParFile(filename):
    """ load the named par file into a dictionary. """

    d = configFile.readfile(filename)
    return d
