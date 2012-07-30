
try:
    set_trace
except:
    try:
        from IPython.core.debugger import Tracer
        set_trace = Tracer()
    except:
        import pdb
        set_trace = pdb.set_trace



    
