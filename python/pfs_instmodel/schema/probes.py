import pfs_tools.par

""" Schema definitions for probes. """

""" The known targeting or simulation object types. """
PROBETYPE = frozenset(('SKY',
                       'OBJECT', 
                       'UNPLUGGED',
                       'SIMFLAT',         #: sims only: we insert a quartz spectrum
                       'SIMCOMB',         #: sims only: we insert a comb spectrum.
                       ))

""" The external structure defining individual probes. """
PROBE = pfs_tools.par.makeParClass(('fiberId', int),
                                   ('ra', float),
                                   ('dec', float),
                                   ('x', float),
                                   ('y', float),
                                   ('type', PROBETYPE))
