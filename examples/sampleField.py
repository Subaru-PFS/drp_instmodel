import pfs_tools.par

PROBETYPE = frozenset(('SKY', 'DUMMYFLAT', 'UNPLUGGED',))
PROBE = pfs_tools.par.makeParClass(('fiberId', int),
                         ('ra', float),
                         ('dec', float),
                         ('x', float),
                         ('y', float),
                         ('type', PROBETYPE))
    
field0 = (
    PROBE(0, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(1, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(2, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(3, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(4, 0.0, 1.0, 100.0, 200.0, 'SKY'),
)

field1 = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(300)]
