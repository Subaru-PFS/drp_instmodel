import pfs_tools.schema

""" Schema definitions for probes. """

""" The known targeting or simulation object types. """
PROBETYPE = frozenset(('SKY',
                       'OBJECT', 
                       'UNPLUGGED',
                       'SIMFLAT',         #: sims only: we insert a quartz spectrum
                       'SIMCOMB',         #: sims only: we insert a comb spectrum.
                       'SIMARC',          #: sims only: we insert an arc spectrum
                       ))

""" The external structure defining individual probes. We need to add some detail field(s), for
object types, targeting info, etc. Or, say, spacing for SIMCOMBs.
"""
PROBE = pfs_tools.schema.makeParClass(('fiberId', int),
                                      ('ra', float),  #: decimal degrees
                                      ('dec', float), #: decimal degrees
                                      ('x', float),   #: mm from center on focal plane
                                      ('y', float),   #: mm from center on focal plane
                                      ('type', PROBETYPE),
                                      ('args', tuple))
