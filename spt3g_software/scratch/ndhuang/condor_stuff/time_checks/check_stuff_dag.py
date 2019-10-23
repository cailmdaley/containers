import sys, os
import glob

arcfiles = glob.glob('/spt/data/arc/*.dat')
arcfiles += glob.glob('/sptpol/2014*.gz')
arcfiles += glob.glob('/sptpol/2015*.gz')
arcfiles += glob.glob('/sptpol/2016*.gz')

f = open('check_times.dag', 'w')
for af in arcfiles:
    bn = os.path.basename(af).split('.')[0]
    f.write('JOB\t{}\t/home/ndhuang/spt_code/spt3g_software/scratch/ndhuang/condor_stuff/time_checks/check.sub\n'.format(bn))
    f.write('VARS\t{}\tINPUT="{}"\n'.format(bn, af))
    f.write('RETRY\t{}\t12\n'.format(bn))
f.close()
