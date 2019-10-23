from mpi4py import MPI

def subdivide_list( lst, n_bundles):
   n_per_bundle = map(lambda x: len(lst)/n_bundles, range(n_bundles))
   for i in range( len(lst) % n_bundles):
      n_per_bundle[i] += 1
   ind = 0
   out_lst = []
   for nper in n_per_bundle:
      out_lst.append(lst[ind:nper + ind])
      ind += nper
   return out_lst

def mpi_mapreduce(data_lst, 
                  map_operation, 
                  reduce_operation, 
                  comm=None,
                  map_args = (),
                  set_seed = True ):
   if comm == None:
      comm = MPI.COMM_WORLD
   mpisize = comm.Get_size()
   rank = comm.Get_rank()

   if set_seed:
      import numpy, time
      numpy.random.seed(seed = int(time.time() * rank) % 4294967293 )
   if rank == 0:
      data = subdivide_list(data_lst, mpisize)
   else:
      data = None
   #print('hello from rank %d'%rank)
   data = comm.scatter(data, root=0)
   #print('%d '%rank, data)
   data = map(lambda x: map_operation(x, *map_args), data)
   #print('Mapped: %d '%rank, data)
   data = reduce(reduce_operation, data)
   #print('Local reduced: %d '%rank, data)
   reduced_data = comm.reduce(sendobj = data, op = reduce_operation)
   #print('Reduced: %d '%rank, reduced_data)
   if rank == 0:
      return reduced_data
   else:
      return None


if __name__ == '__main__':
   def red(x,y):
      x += y
      return x
      
   def mp(x, v1, v2):
      x.val *= v1 * v2
      return x

   class test(object):
      def __init__(self, val):
         self.val = val
      def __iadd__(self, other):
         self.val += other.val
         return self
   data = [test((i+1)**2) for i in range(2)]
   out_data = mpi_mapreduce(data, mp, red, map_args = (1,1))
   if out_data != None:
      print out_data.val
