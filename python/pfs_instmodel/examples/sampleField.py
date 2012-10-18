from pfs_instmodel.schema.probes import PROBE
    
field0 = (
    PROBE(0, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(1, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(2, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(3, 0.0, 1.0, 100.0, 200.0, 'SKY'),
    PROBE(4, 0.0, 1.0, 100.0, 200.0, 'SKY'),
)

field1 = ([PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(0,5)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(5,50)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(50,55)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(55,100)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(100,105)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(105,290)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(290,301)])

