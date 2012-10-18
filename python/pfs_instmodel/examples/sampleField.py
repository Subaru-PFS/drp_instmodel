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
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(105,295)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(295,301)])

quickfield = ([PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(0,2)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(2,10)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(10,12)])

combField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(0,301)]
flatField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(0,301)]
skyField = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(0,301)]

sparseCombField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMCOMB' if i%10 == 0 else 'UNPLUGGED')) for i in range(0,301)]
sparseFlatField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMFLAT' if i%10 == 0 else 'UNPLUGGED')) for i in range(0,301)]
