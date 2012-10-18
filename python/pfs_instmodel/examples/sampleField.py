from pfs_instmodel.schema.probes import PROBE
    
field0 = (
    PROBE(0, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(1, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(2, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(3, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(4, 0.0, 1.0, 100.0, 200.0, 'SKY'),
)

field1 = ([PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(0,2)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(2,100)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(100,102)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(102,298)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(298,300)])

