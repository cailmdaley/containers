import sys, os
import argparse
import glob
import subprocess
sys.path.append('/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang')
import nhutils as nhu

parser = argparse.ArgumentParser()
parser.add_argument('obsids', nargs = '+', type = int)
parser.add_argument('--submit', action = 'store_true')
pargs = parser.parse_args()
obsids = pargs.obsids
with open('/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/condor_base.sub', 'r') as f:
    sub = f.read()
cpus = 1
disk = '2G'
mem = '4G'
logname = 'lr_pointing_map'
exec = '/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/LR_pointing/run.sh'
outdir = '/poleanalysis/ndhuang/LR_pointing/RCW38/'
if not os.path.exists(outdir):
    os.makedirs(outdir)
extras = ''
logdir = '/home/ndhuang/condor/log/{}/'.format(logname)
if not os.path.exists(logdir):
    print(logdir)
    os.makedirs(logdir)
outlog = os.path.join(logdir, '$(cluster).out')
errorlog = os.path.join(logdir, '$(cluster).err')
condorlog = os.path.join(logdir, 'condor.log')
for id in obsids:
    outfile = os.path.join(outdir, str(id) + '.g3')
    subfile = os.path.join('/home/ndhuang/condor/submit', 
                           logname + '_' + str(id) + '.sub')
    rcw38_cal = nhu.find_cal(id, '/poleanalysis/sptdaq/calresult/calibration/RCW38-pixelraster/')
    infiles = sorted(glob.glob(os.path.join('/spt_data/bolodata/fullrate/RCW38-pixelraster',
                                            str(id), '0*.g3')))
    offline_cal = os.path.join('/spt_data/bolodata/fullrate/RCW38-pixelraster',
                               str(id), 'offline_calibration.g3')
    args = "{offline_cal} {rcw38_cal} {infile} -o {outfile}".format(
        offline_cal = offline_cal, rcw38_cal = rcw38_cal, 
        infile = ' '.join(infiles), outfile = outfile)
    substr = sub.format(exec = exec,
                        args = args,
                        # inputfiles = '/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/coadd_LR_sourcemap.py',
                        inputfiles = '',
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
