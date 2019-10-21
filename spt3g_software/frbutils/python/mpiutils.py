def do_collective_analysis(filenames, output_pkl_fn, 
                           anal_class, *anal_args):
    '''
    Anal class is a G3Module with a defined += operator
    anal args are the arguments used to instantiate it
    filenames are a list of files the module applies operations to
    '''

    from spt3g.cpuclusterutils.mpimapreduce import mpi_mapreduce
    def _map_helper(filename, anal_class, anal_args):
        instantiated_class = anal_class(*anal_args)
        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader, filename = filename)
        pipe.Add(instantiated_class)
        pipe.Run()
        return instantiated_class
    def _reduce_helper(x,y):
        x += y
        return x
    data  = mpi_mapreduce(filenames, _map_helper, _reduce_helper, 
                          map_args = (anal_class, anal_args))
    if data != None:
        print('saving')
        pickle.dump(data, open(output_pkl_fn, 'w'))
    else:
        print('not saving')

def do_inputless_collective_analysis(n_runs, output_pkl_fn, 
                                     anal_func, anal_args):
    from spt3g.cpuclusterutils.mpimapreduce import mpi_mapreduce
    def _map_helper(dummy):
        print('maphelper', anal_args)
        return anal_func(*anal_args)
    def _reduce_helper(x,y):
        print('reducehelper')
        x += y
        return x
    data  = mpi_mapreduce(range(n_runs), _map_helper, _reduce_helper)
    if data != None:
        print('saving')
        pickle.dump(data, open(output_pkl_fn, 'w'))
    else:
        print('not saving')



