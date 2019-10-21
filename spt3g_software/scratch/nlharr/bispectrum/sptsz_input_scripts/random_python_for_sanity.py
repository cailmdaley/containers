import numpy as np
import pylab as pl
import random, os.path



def store_divided_runlist(rlst, dir, n_store):
    f = open(rlst)
    mps = []
    for line in f:
        mps.append(f.strip)
    f.close()
    for i in xrange(n_store):
        outfile_base = dir+'/'+os.path.basename(rlst)+('noise_split_%d'%i)
        random.shuffle(mps)
        lst_a = reduce(lambda x,y: x+'\n'+y, mps[:len(mps)//2])+'\n'
        lst_b = reduce(lambda x,y: x+'\n'+y, mps[len(mps)//2]:)+'\n'
        #open(outfile_base+'_a','w').write(lst_a)
        #open(outfile_base+'_b','w').write(lst_b)
        print "I would write to outfile", outfile_base+'_a', outfile_base+'_b'


if __name__ == '__main__':
    pass
    
