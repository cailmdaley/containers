import glob, os, sys, numpy, subprocess
import cPickle as pickle
import numpy as np
centerpos = pickle.load(open('source_coords.p','rb'))

obslist = [k for k in centerpos if k[0] == '4']
obslist = ['47038378']

halfway = len(obslist)//2
scottkeys = obslist[:halfway]
amundkeys = obslist[halfway:]
scottkeys = obslist

pmnj_rt = '/spt/data/bolodata/downsampled/PMNJ0210-5101-pixelraster/'
f0537_rt = '/spt/data/bolodata/downsampled/0537-441-pixelraster/'

pmnj_opts = '-s PMNJ0210-5101 -r 0.2 -x 0.5 -y 0.5 -p 4 -sm -pw'
f0537_opts = '-s 0537-441 -r 0.2 -x 0.5 -y 0.5 -p 4 -sm -pw'

exe = 'python sourcemaps_pointing.py '

hostname = subprocess.check_output(['hostname']).rstrip()

for obs in sorted(centerpos.keys()):
    fobs = int(obs)
    data_rt = f0537_rt if fobs < 33000000 else pmnj_rt
    calframe = data_rt + obs +'/offline_calibration.g3 '
    if obs[:3] == '354': 
        calframe = '/home/sguns/sourcemaps/mylittlecalframe.g3 '
    inframes = sorted(glob.glob(data_rt + obs + '/0*.g3'))
    out = '/spt/user/sguns/focusquasar/maps/static/Kcmb/pointing_afterevent/' + obs +'.g3'
    log = '/home/sguns/sourcemaps/staticmaps/log/' + obs + '.out'
    err = '/home/sguns/sourcemaps/staticmaps/log/' + obs + '.err'
    cmd = exe + calframe
    for fr in inframes:
        cmd += fr + ' '
    cmd += '-o ' + out
    if fobs < 33000000:
        cmd += ' ' + f0537_opts
    else:
        cmd += ' ' + pmnj_opts
    cmd += ' -ra \'%.5f\' -dec \'%.5f\'' %(centerpos[obs][0], centerpos[obs][1])
    cmd += ' > %s 2> %s &' %(log, err)
    if hostname == 'scott.grid.uchicago.edu' and obs in scottkeys: 
        print "launching", obs
        print cmd
        #os.system(cmd)
    if hostname == 'amundsen.grid.uchicago.edu' and obs in amundkeys: 
        print "launching", obs
        os.system(cmd)
