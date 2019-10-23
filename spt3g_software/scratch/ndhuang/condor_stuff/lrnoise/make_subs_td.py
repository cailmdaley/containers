import sys, os
import argparse
import glob
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('source')
parser.add_argument('obsids', nargs = '+', type = int)
parser.add_argument('--submit', action = 'store_true')
pargs = parser.parse_args()
obsids = pargs.obsids
source = pargs.source
f = open('/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/condor_base.sub', 'r')
sub = f.read()
f.close()
cpus = 1
disk = '64G'
mem = '2G'
logname = 'lrnoise_td'
exec = '/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/lrnoise/run_td.sh'
outdir = '/spt/user/ndhuang/source_noise/'
extras = ''
logdir = '/home/ndhuang/condor/log/{}/'.format(logname)
if not os.path.exists(logdir):
    os.makedirs(logdir)
outlog = os.path.join(logdir, '$(cluster).out')
errorlog = os.path.join(logdir, '$(cluster).err')
condorlog = os.path.join(logdir, 'condor.log')
for id in obsids:
    outdir = os.path.join(outdir, str(id))
    subfile = os.path.join('/home/ndhuang/condor/submit', 
                           logname + '_' + str(id) + '.sub')
    infiles = sorted(glob.glob(os.path.join('/spt/data/bolodata/downsampled/', source, 
                                            str(id), '0*.g3')))
    calfile = os.path.join('/spt/user/production/calibration/calframe', source, 
                           str(id) + '.g3')
    args = "'{calfile} {infile}' '{outfile}'".format(calfile = calfile,
          infile = ' '.join(infiles), outfile = outdir)
    substr = sub.format(exec = exec,
                        args = args,
                        inputfiles = '/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/lrnoise.py',
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
