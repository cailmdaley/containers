import glob, os, sys, numpy, subprocess
import cPickle as pickle


obslist = glob.glob('/spt/data/bolodata/downsampled/PMNJ0210-5101-pixelraster/47*')
obslist = [o.split('.')[0].split('/')[-1] for o in obslist]
print "Obslist: ", obslist

halfway = len(obslist)//2
scottkeys = obslist[:halfway]
amundkeys = obslist[halfway:]

pmnj_rt = '/spt/data/bolodata/downsampled/PMNJ0210-5101-pixelraster/'
f0537_rt = '/spt/data/bolodata/downsampled/0537-441-pixelraster/'
exe = 'python sourcemaps_poly_mask_lowvcut.py ' 

pmnj_opts = '-s PMNJ0210-5101 -r 0.2 -x 0.5 -y 0.5 -p 4 -dm'
f0537_opts = '-s 0537-441 -r 0.2 -x 0.5 -y 0.5 -p 4 -dm'

hostname = subprocess.check_output(['hostname']).rstrip()

for obs in sorted(obslist):
    fobs = int(obs)
    data_rt = f0537_rt if fobs < 33000000 else pmnj_rt
    calframe = data_rt + obs +'/offline_calibration.g3 '
    fixedcal = '/home/sguns/mylittlecalframe.g3 '
    #print("WARNING: proceeding with fixed calframe", fixedcal)
    #calframe = fixedcal
    inframes = sorted(glob.glob(data_rt + obs + '/0*.g3'))
    out = '/spt/user/sguns/focusquasar/maps/dynamic/Kcmb/' + obs +'.g3'
    log = '/home/sguns/sourcemaps/dynamicmaps/log/' + obs + '.out'
    err = '/home/sguns/sourcemaps/dynamicmaps/log/' + obs + '.err'
    cmd = exe + calframe
    for fr in inframes:
        cmd += fr + ' '
    cmd += '-o ' + out
    if fobs < 33000000:
        cmd += ' ' + f0537_opts
    else:
        cmd += ' ' + pmnj_opts
    cmd += ' > %s 2> %s &' %(log, err)
#    print("WARNING: launching all jobs on current machine")
    if hostname == 'scott.grid.uchicago.edu' and obs in scottkeys:
        print "launching " + obs
        os.system(cmd)
    if hostname == 'amundsen.grid.uchicago.edu' and obs in amundkeys:
        print "launching " + obs
        os.system(cmd)
