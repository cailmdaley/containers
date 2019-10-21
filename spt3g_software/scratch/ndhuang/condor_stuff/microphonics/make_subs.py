import sys, os
import argparse
import glob
import pickle
import subprocess
import numpy as np
from spt3g import core

def find_arc_for_time(start, stop, arcfile_times):
    i_start = -2
    i_stop = -2
    for i, at in enumerate(arcfile_times):
        if at.time < start.time:
            i_start = i
        if at.time > stop.time:
            i_stop = i + 1
            break
    assert(i_start >= 0 and i_stop >= 0)
    return i_start, i_stop

parser = argparse.ArgumentParser()
parser.add_argument('source')
parser.add_argument('--submit', action = 'store_true')
pargs = parser.parse_args()
source = pargs.source
f = open('/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/condor_base.sub', 'r')
sub = f.read()
f.close()
cpus = 1
disk = '1G'
mem = '2G'
logname = 'microphonics'
exec = '/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/microphonics/run_microphonics.sh'
# outdir = '/spt/user/ndhuang/microphonics'
extras = ''
logdir = '/home/ndhuang/condor/log/{}/'.format(logname)
if not os.path.exists(logdir):
    os.makedirs(logdir)
outlog = os.path.join(logdir, '$(cluster).out')
errorlog = os.path.join(logdir, '$(cluster).err')
condorlog = os.path.join(logdir, 'condor.log')

f = open('/spt/user/ndhuang/obsid_lut.pkl', 'rb')
oblut = pickle.load(f)
f.close()
obsids = sorted(os.listdir(os.path.join('/spt/data/bolodata/downsampled/',
                                        source)), reverse = True)
arcfiles = np.sort(glob.glob('/spt/data/arc/*.dat'))
arcfile_times = np.array([core.G3Time(os.path.basename(af).split('.')[0]) for af in arcfiles])
for id in obsids:
    id = int(id)
    if id not in oblut.keys():
        print(id)
        continue
    i_start, i_stop =find_arc_for_time(oblut[id]['start'], oblut[id]['stop'], arcfile_times)
    infiles = arcfiles[i_start:i_stop]
    outfile = str(id) + '_temps.pkl'
    subfile = os.path.join('/home/ndhuang/condor/submit', 
                           logname + '_' + str(id) + '.sub')
    args = "'{infile}' '{outfile}'".format(
          infile = ' '.join(infiles), outfile = outfile)
    substr = sub.format(exec = exec,
                        args = args,
                        inputfiles = '/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/nhutils.py,/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/microphonics/get_temperatures.py',
                        outputfiles = '', 
                        extras = extras,
                        outlog = outlog,
                        errorlog = errorlog,
                        condorlog = condorlog,
                        cpus = cpus, disk = disk, mem = mem)
    f = open(subfile, 'w')
    f.write(substr)
    f.close()
    if pargs.submit:
        subprocess.check_call(['condor_submit', subfile])
