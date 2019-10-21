import sys, os
import glob
import subprocess

obsids = sys.argv[1:]
f = open('/home/ndhuang/spt_code/scritch/condor_stuff/condor_base.sub', 'r')
sub = f.read()
f.close()
for id in obsids:
    infiles = sorted(glob.glob(os.path.join('/spt/data/bolodata/fullrate/RCW38', 
                                            str(id), '0*.g3')))
    outdir = os.path.join('/spt/user/ndhuang/bolodata/downsampled/RCW38', 
                          str(id))
    args = "'{infile}' '{outfile}'".format(
          infile = ' '.join(infiles), outfile = outdir)
    extras = '\n'.join(['request_disk = 30G',
                        'request_memory = 1G'])
    substr = sub.format(exec = '/home/ndhuang/spt_code/scritch/condor_stuff/run_downsample.sh',
                        args = args,
                        inputfiles = ',/home/ndhuang/spt_code/spt3g_software/pole_processing/downsample_bolodata.py',
                        outputfiles = '', extras = extras)
    f = open('/home/ndhuang/condor/submit/downsample.sub', 'w')
    f.write(substr)
    f.close()
    subprocess.check_call(['condor_submit', 
                           '/home/ndhuang/condor/submit/downsample.sub'])
