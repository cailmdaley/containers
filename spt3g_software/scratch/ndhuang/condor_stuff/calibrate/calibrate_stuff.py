import sys, os
import subprocess

obsids = sys.argv[1:]
f = open('/home/ndhuang/spt_code/scritch/condor_stuff/condor_base.sub', 'r')
sub = f.read()
f.close()
for id in obsids:
    infile = os.path.join('/spt/data/bolodata/fullrate/calibrator/', 
                          str(id), '0000.g3')
    if not os.path.exists(infile):
        continue
    outfile = os.path.join('/spt/user/ndhuang/calibration/calibrator_fullrate', 
                          str(id) + '.g3')
    args = "{infile} {outfile} \'{local_in} {local_out}\'".format(
          infile = infile, outfile = outfile, 
          local_in = os.path.basename(infile), 
          local_out = os.path.basename(outfile))
    substr = sub.format(exec = '/home/ndhuang/spt_code/scritch/condor_stuff/run_calibrate.sh',
                        args = args,
                        inputfiles = ',/home/ndhuang/spt_code/spt3g_software/calibration/scripts/analyze_calibrator.py',
                        outputfiles = '', extras = '\n'.join(['request_memory = 3G']))
    f = open('/home/ndhuang/condor/submit/downsample.sub', 'w')
    f.write(substr)
    f.close()
    subprocess.check_call(['condor_submit', 
                           '/home/ndhuang/condor/submit/downsample.sub'])
