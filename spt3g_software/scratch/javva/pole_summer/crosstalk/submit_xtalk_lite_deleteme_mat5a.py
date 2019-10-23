import os
import glob
from spt3g.cluster.condor_tools import condor_submit
import argparse
from spt3g import core, mapmaker
import pickle 
import numpy
import time
#import make_rcw38_template

parser = argparse.ArgumentParser()
#parser.add_argument('obsids', nargs = '+', type = str)
parser.add_argument('--submit', action = 'store_true')
parser.add_argument('-s','--source', action ='store', 
                    default='ra0hdec-57.5')
obsnum = '65126904'
pargs = parser.parse_args()
#obsids = pargs.obsids

job_root = 'MAT5a_lstsq_%s'%obsnum
condor_dir = os.path.join('/home/javva/grid/condor_logs', pargs.source, job_root)
out_root = os.path.join('/spt/user/javva/crosstalk/rcw38_singlebolo_coadds_analysis_amps', job_root)
script = '/home/javva/spt3g_software/scratch/javva/pole_summer/crosstalk/rcw38_fit_xtalk_grid_lite_subset_abs.py'

test = True
if pargs.submit:
    test = False

def make_template(obsnum):
    os.system('python make_rcw38_template.py /spt/user/production/calibration/calframe/MAT5A-pixelraster/%s.g3 /spt/user/production/calibration/MAT5A-pixelraster/maps/%s.g3 /spt/data/bolodata/downsampled/MAT5A-pixelraster/%s/0000.g3 -o /spt/user/javva/crosstalk/templates/MAT5A_%s.pkl' %(obsnum,obsnum,obsnum,obsnum))

def make_split_map(obsnum):
    map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename='/spt/user/production/calibration/MAT5A-pixelraster/maps/%s.g3'%obsnum)
    print('got here')
    pipe.Add(lambda fr: 'Id' not in fr or 'Wunpol' in fr or (numpy.asarray(fr['T']) != 0).any())
    pipe.Add(map_extractor)
    pipe.Run()
    print('ran the pipe')
    data = map_extractor.maps
    for b in data.keys():
        if 'Wunpol' in data[b] and b != 'bsmap':
            w = numpy.asarray(data[b]['Wunpol'].TT)[:]
    for row in range(len(w)):
        trans = numpy.where(numpy.diff((w[row] == 0).astype('float')))[0]
        if len(trans) >= 2:
            w[row][:trans[0] + 5] = 0
            w[row][trans[-1] - 5:] = 0
        elif len(trans) == 1:
            if trans[0] > len(w[row])/2:
                w[row][trans[0] - 5:] = 0
            else:
                w[row][:trans[0] + 5] = 0

    for b in data.keys():
        if not os.path.isfile('/spt/user/javva/crosstalk/ind_maps/MAT5A_%s_%s.pkl'%(obsnum, b)):
            maps = {}
            rawmap_b = numpy.asarray(data[b]['T'])[:]
            maps['data'] = rawmap_b
            maps['w'] = w
            with open('/spt/user/javva/crosstalk/ind_maps/MAT5A_%s_%s.pkl'%(obsnum, b), 'wb') as handle:
                pickle.dump(maps, handle, protocol=pickle.HIGHEST_PROTOCOL)

requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX") )'''
baddies = 0

if '/spt/user/javva/crosstalk/templates/MAT5A_'+obsnum+'.pkl' not in glob.glob('/spt/user/javva/crosstalk/templates/*'):
    print("Making the template")
    make_template(obsnum)
    print('finished making the new template')

a = list(core.G3File('/mnt/ceph/srm/spt3g/user/production/calibration/MAT5A-pixelraster/maps/%s.g3'%obsnum))
print('This is the len of a',len(a))
for i in a:
    if 'Id' not in i:
        print('not a map')
        print(i)
        continue

    try:
        main_b = i['Id']
        mbb = main_b
        if not os.path.isfile('/spt/user/javva/crosstalk/ind_maps/MAT5A_%s_%s.pkl'%(obsnum, main_b)):
            print('in this directory, there are %s maps that match the obsid'%len(glob.glob('/spt/user/javva/crosstalk/ind_maps/MAT5A_%s_*'%obsnum)))
            if len(glob.glob('/spt/user/javva/crosstalk/ind_maps/MAT5A_%s_*'%obsnum))<1000:
                   print('making split maps')
                   make_split_map(obsnum)
                   print('made split map')
        print("file is", os.path.isfile('/spt/user/javva/crosstalk/ind_maps/MAT5A_%s_%s.pkl'%(obsnum, main_b)))
        jobname = job_root+'_'+str(mbb)
        infiles = ['/spt/user/javva/crosstalk/templates/MAT5A_%s.pkl'%obsnum, '/spt/user/javva/crosstalk/ind_maps/MAT5A_%s_%s.pkl'%(obsnum, main_b)]
        args_in1 = [os.path.basename(dat) for dat in infiles]
        args_in = ['./'+d for d in args_in1]
        #    extra_args= #'-x 75 -y 50' # --lr'
        args = '{infiles} -mb {main_bolo} -o {outfile}'.format(infiles = ' '.join(args_in),main_bolo = mbb, outfile = obsnum+'_'+jobname+'.pkl')
        print(args)

        condor_submit(script, create_only=test, args = [args],
                      log_root = os.path.join(condor_dir, str(mbb)),
                  output_root = os.path.join(out_root),retry = True, 
                  jobname = jobname,
                  input_files= infiles,
                  clustertools_version='py3-v2',
                  output_files=[obsnum+'_'+jobname+'.pkl'],
                  requirements = requirements,
                  request_disk=4*core.G3Units.GB,
                  request_memory=4*core.G3Units.GB,
                  grid_proxy='/home/javva/javva_proxy')
        time.sleep(2)
    except:
        print('THIS OBSID DIDNT WORK ==',mbb)
        baddies = baddies+1

print("end of processing. didn't complete ",baddies)
