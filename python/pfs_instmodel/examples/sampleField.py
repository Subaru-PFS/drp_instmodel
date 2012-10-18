from pfs_instmodel.schema.probes import PROBE
    
field1 = ([PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(0,5)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(5,50)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(50,55)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(55,100)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(100,105)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(105,295)] +
          [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(295,301)])

# Minimal sanity test field.
quickfield = ([PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(-300,-298)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(-298,-21)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(-21,-19)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(-19,-11)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(-11,-9)] +

              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(-11,-2)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(-2,3)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(3,10)] +

              [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(10,12)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(12,20)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(20,22)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'UNPLUGGED') for i in range(22,299)] +
              [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(299,301)])

combField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMCOMB') for i in range(-300,301)]
flatField = [PROBE(i,0.0,1.0,100.0,200.0,'SIMFLAT') for i in range(-300,301)]
skyField = [PROBE(i,0.0,1.0,100.0,200.0,'SKY') for i in range(-300,301)]

sparseCombField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMCOMB' if i%10 == 0 else 'UNPLUGGED')) for i in range(-300,301)]
sparseFlatField = [PROBE(i,0.0,1.0,100.0,200.0,('SIMFLAT' if i%10 == 0 else 'UNPLUGGED')) for i in range(-300,301)]
