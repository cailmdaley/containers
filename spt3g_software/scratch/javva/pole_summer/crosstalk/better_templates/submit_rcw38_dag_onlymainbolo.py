import os
import glob
from spt3g.cluster.condor_tools import condor_submit
import argparse
from spt3g import core, mapmaker, dfmux
import pickle 
import numpy
import time


'''
This is a script to submit grid jobs to make the crosstalk matrix off of rcw38
'''

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true')
parser.add_argument('-s','--source', action ='store', 
                    default='ra0hdec-57.5')
parser.add_argument('-on','--obsnum', action ='store', 
                    default=None)

#obsnum = '70792281' #oroginal good one for rcw38
#obsnum = '63643116' #oroginal good one for rcw38
#obsnum = '64146996' #second ok one
pargs = parser.parse_args()
obsnum = pargs.obsnum
#bookkeeping stuff for the output file
job_root = 'injected_source_%s_smset_fasttry'%obsnum
condor_dir = os.path.join('/home/javva/grid/condor_logs', pargs.source, job_root)
out_root = os.path.join('/spt/user/javva/crosstalk/fit_outputs', job_root)
dag_script = os.path.join(condor_dir, job_root+'_atten' + '.dag')
print(dag_script)

bd = '/spt/user/javva'


#script to use to do the fit
script = '/home/javva/spt3g_software/scratch/javva/pole_summer/crosstalk/better_templates/rcw38_fit_better_templates_HARD_CODED_OFFSET_only_main_source_injected_source.py'
test = True
if pargs.submit:
    test = False

#If a template doesn't exist for this observation number, this function submits a job to make it
def make_template(obsnum):
    os.system('python template_to_pkl.py /spt/user/production/calibration/calframe/RCW38-pixelraster/%s.g3 /spt/user/javva/crosstalk/templates/SingleWaferMapsAmps_%s_time_constant_deconvolved_pt5_arcmin.g3 /spt/data/bolodata/downsampled/RCW38-pixelraster/%s/0000.g3 -o /spt/user/javva/crosstalk/templates/%s.pkl' %(obsnum,obsnum,obsnum,obsnum))

def make_fine_coadd(obsnum):
    os.system('python make_fine_coadd.py -on %s'%obsnum)

def make_ind_bolo_g3s(obsnum):
    os.system('python makemaps_vanilla_dc.py /spt/user/production/calibration/calframe/RCW38-pixelraster/%s.g3 /spt/data/bolodata/downsampled/RCW38-pixelraster/%s/0* -o %s/crosstalk/ind_maps/ind_bolo_maps_%s.g3 -s rcw38'%(obsnum,obsnum, bd, obsnum))



bp = None
dfm = None
t_f = None
calresponse = None
calresponsesn = None
if not os.path.exists(condor_dir):
    os.mkdir(condor_dir)
f = open(dag_script,'w')
print('made the dag file')

#This definition takes the autoprocessing .g3 single bolo maps file and splits them into individual maps, so you can run the fit on the grid with less memory use
def make_split_map(obsnum):
    def extract_bp(fr):
        global dfm, t_f
        if 'DfMuxHousekeeping' in fr and dfm is None:
            dfm = fr['DfMuxHousekeeping']
            print('found dfm')
        if 'DfMuxTransferFunction' in fr and t_f is None:
            t_f = fr['DfMuxTransferFunction']
        global bp, calresponse, calresponsesn,wm
        if 'BolometerProperties' in fr and bp is None:
            bp = fr['BolometerProperties']
        if 'CalibratorResponse' in fr:
            calresponse = fr['CalibratorResponse']
            calresponsesn = fr['CalibratorResponseSN']
        if 'WiringMap' in fr:
            wm = fr['WiringMap']
            print('found the wm')
    #Extracts the maps from the large file, and grabs housekeeping information
    map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
    pipe = core.G3Pipeline()
    ifs = ['/spt/user/production/calibration/calframe/RCW38-pixelraster/'+obsnum+'.g3', '/spt/user/javva/crosstalk/ind_maps/ind_bolo_maps_'+obsnum+'.g3', '/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsnum+'/0000.g3']
    pipe.Add(core.G3Reader, filename=ifs)
    pipe.Add(extract_bp)
    pipe.Add(lambda fr: 'Id' not in fr or 'Wunpol' in fr or (numpy.asarray(fr['T']) != 0).any())
    pipe.Add(map_extractor)
    pipe.Run()
    print('ran the pipe')
    #The maps are in watts, so we need to change them to amps 
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
    #Saves the single bolo map and weights map for each bolometer
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

if bd +'/crosstalk/templates/'+obsnum+'.pkl' not in glob.glob(bd+'/crosstalk/templates/*'):
    if not os.path.isfile(bd+'/crosstalk/templates/SingleWaferMapsAmps_%s_time_constant_deconvolved_pt5_arcmin.g3'%obsnum):
        print('Making the fine coadd')
        make_fine_coadd(obsnum)
    print("Making the template")
    make_template(obsnum)
    print('finished making the new template')
if len(glob.glob('%s/crosstalk/ind_maps/%s_*'%(bd,obsnum)))<1000:
    if not os.path.isfile(bd+'/crosstalk/ind_maps/ind_bolo_maps_%s.g3'%obsnum):
            print('making the deconvolved individual maps')
            make_ind_bolo_g3s(obsnum)
    print('making split maps')
    make_split_map(obsnum)
    print('made split maps')


#read in the single bolo maps file to know which bolos are alive in the observation
a = list(core.G3File('/mnt/ceph/srm/spt3g/user/production/calibration/RCW38-pixelraster/maps/%s.g3'%obsnum))

for i in a:
    if 'Id' not in i:
        continue
    try:
        main_b = i['Id']
        print('main b is', main_b)
#        if os.path.isfile(out_root+'/'+obsnum+'_'+job_root+'_'+main_b+'.pkl'):
#            print('File already ran')
#            continue
        print('The main bolo is %s'%main_b)
        mbb = main_b
        # if not os.path.isfile('/spt/user/javva/crosstalk/ind_maps/%s_%s.pkl'%(obsnum, main_b)):
        #     if len(glob.glob('/spt/user/javva/crosstalk/ind_maps/%s_*'%obsnum))<1000:
        #         print('making split maps')
        #         make_split_map(obsnum)
        #         print('made split maps')
        jobname = job_root+'_'+str(mbb)
        # if os.path.isfile(os.path.join(out_root,obsnum+'_'+jobname+'.pkl')):
        #     print('file exists')
        #     continue
        infiles = ['/spt/user/javva/crosstalk/templates/%s.pkl'%obsnum, '/spt/user/javva/crosstalk/ind_maps/%s_%s.pkl'%(obsnum, main_b)]
        args_in1 = [os.path.basename(dat) for dat in infiles]
        args_in = ['./'+d for d in args_in1]
        args = '{infiles} -mb {main_bolo} -o {outfile}'.format(infiles = ' '.join(args_in),main_bolo = mbb, outfile = obsnum+'_'+jobname+'.pkl')
        print(args)

        condor_submit(script, create_only=True, args = [args],
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
        print('got to write bit')
        print(condor_dir, str(mbb), jobname+'.submit')
        f.write('JOB %s %s\n'%(jobname, os.path.join(condor_dir, str(mbb), jobname+'.submit')))
        f.write('RETRY %s 5\n'%jobname)

    except:
        print('THIS OBSID DIDNT WORK ==',mbb)
        baddies = baddies+1

print("end of processing. didn't complete ",baddies)

f.close()
if pargs.submit:
    os.system('condor_submit_dag -update_submit -maxidle 25 %s'%dag_script)
    print('submitted script')

