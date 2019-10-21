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
#obsnum = '63643116' #oroginal good one
obsnum = '64146996'
pargs = parser.parse_args()
#obsids = pargs.obsids

main_bolos = pickle.load(open('/home/javva/spt3g_software/scratch/javva/pole_summer/crosstalk/allbolos.pkl', 'rb'))

job_root = 'xtalk_one_map_lite_subset_abs_only_one_from_pixel_lstsq%s'%obsnum
condor_dir = os.path.join('/home/javva/grid/condor_logs', pargs.source, job_root)
out_root = os.path.join('/spt/user/javva/crosstalk/rcw38_singlebolo_coadds_analysis_amps', job_root)
script = '/home/javva/spt3g_software/scratch/javva/pole_summer/crosstalk/rcw38_fit_xtalk_grid_lite_subset_abs.py'

test = True
if pargs.submit:
    test = False

def make_template(obsnum):
    os.system('python make_rcw38_template.py /spt/user/production/calibration/calframe/RCW38-pixelraster/%s.g3 /spt/user/production/calibration/RCW38-pixelraster/maps/%s.g3 /spt/data/bolodata/downsampled/RCW38-pixelraster/%s/0000.g3 -o /spt/user/javva/crosstalk/templates/%s.pkl' %(obsnum,obsnum,obsnum,obsnum))

def make_split_map(obsnum):
    map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
    pipe = core.G3Pipeline()
    ifs = ['/spt/user/production/calibration/calframe/RCW38-pixelraster/'+obsnum+'.g3', '/spt/user/production/calibration/RCW38-pixelraster/maps/'+obsnum+'.g3', '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsnum+'/0000.g3']
    pipe.Add(core.G3Reader, filename=ifs)
    def extract_bp(fr):
    global dfm
    if 'DfMuxHousekeeping' in fr and dfm is None:
        dfm = fr['DfMuxHousekeeping']
        print('found dfm')
    global bp, calresponse, calresponsesn,wm
    if 'BolometerProperties' in fr and bp is None:
        bp = fr['BolometerProperties']
    if 'CalibratorResponse' in fr:
        calresponse = fr['CalibratorResponse']
        calresponsesn = fr['CalibratorResponseSN']
    if 'WiringMap' in fr:
        wm = fr['WiringMap']
        print('found the wm')
    pipe.Add(lambda fr: 'Id' not in fr or 'Wunpol' in fr or (numpy.asarray(fr['T']) != 0).any())
    pipe.Add(map_extractor)
    pipe.Run()
    print('ran the pipe')
    watt_data = map_extractor.maps
    amp_data = {}
    for i in watt_data.keys():
        if i == 'bsmap':
            amp_data[i] = watt_data[i]
            continue
    amp_data[i] = {}
    amp_data[i]['T'] = watt_data[i]['T']/dfmux.HousekeepingForBolo(dfm,wm, i).carrier_frequency
    data = amp_data
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
        if not os.path.isfile('/spt/user/javva/crosstalk/ind_maps/%s_%s.pkl'%(obsnum, b)):
            maps = {}
            rawmap_b = numpy.asarray(data[b]['T'])[:]
            maps['data'] = rawmap_b
            maps['w'] = w
            with open('/spt/user/javva/crosstalk/ind_maps/%s_%s.pkl'%(obsnum, b), 'wb') as handle:
                pickle.dump(maps, handle, protocol=pickle.HIGHEST_PROTOCOL)

requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX") )'''
baddies = 0

if '/spt/user/javva/crosstalk/templates/'+obsnum+'.pkl' not in glob.glob('/spt/user/javva/crosstalk/templates/*'):
    print("Making the template")
    make_template(obsnum)
    print('finished making the new template')
a = list(core.G3File('/mnt/ceph/srm/spt3g/user/production/calibration/RCW38-pixelraster/maps/%s.g3'%obsnum))

for i in a:
    if 'Id' not in i:
        continue
    try:
        main_b = i['Id']
        print('The main bolo is %s'%main_b)
        mbb = main_b
        if not os.path.isfile('/spt/user/javva/crosstalk/ind_maps/%s_%s.pkl'%(obsnum, main_b)):
            if len(glob.glob('/spt/user/javva/crosstalk/ind_maps/%s_*'%obsnum))<1000:
                print('making split maps')
                make_split_map(obsnum)
                print('made split maps')
        jobname = job_root+'_'+str(mbb)
        if os.path.isfile(os.path.join(out_root,obsnum+'_'+jobname+'.pkl')):
            print('file exists')
            continue
        infiles = ['/spt/user/javva/crosstalk/templates/%s.pkl'%obsnum, '/spt/user/javva/crosstalk/ind_maps/%s_%s.pkl'%(obsnum, main_b)]
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
