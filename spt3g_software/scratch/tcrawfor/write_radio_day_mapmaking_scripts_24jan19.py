from spt3g import core, std_processing
import numpy as np
import glob, os

# calls should look like this: 
#python makecoadd_hack2019.py  /poleanalysis/sptdaq/calresult/calibration/elnod/64939562.g3  /poleanalysis/sptdaq/calresult/calibration/calibrator/64939617.g3  /poleanalysis/sptdaq/calresult/calibration/boloproperties/62679602.g3  /poleanalysis/sptdaq/calresult/calibration/RCW38-pixelraster/64599638.g3 /spt_data/bolodata/fullrate/J1337-1257/64939693/0*.g3 -o /poleanalysis/tcrawfor/script_output/64939693_test4.g3 -s J1337-1257 -k

outdir = '/home/tcrawfor/spt_code/spt3g_software/calibration/scripts/fluxpointcal/'
mapdir = '/poleanalysis/tcrawfor/radio_day_2019/'
bdata_dir = '/spt_data/bolodata/fullrate/'
calresult_dir = '/poleanalysis/sptdaq/calresult/calibration/'
pcall_str = 'nohup python makecoadd_hack2019.py '

# range over which to look for files
obsid_min = std_processing.time_to_obsid(core.G3Time('190123 01:00:16'))
obsid_max = std_processing.time_to_obsid(core.G3Time('190123 10:55:46'))
#obsid_min = std_processing.time_to_obsid(core.G3Time('190123 20:58:06'))
#obsid_max = std_processing.time_to_obsid(core.G3Time('190124 01:28:57'))
#obsid_min = 60000000
#obsid_max = 70000000

# cheap way to get source names
rdfiles_2019_all = glob.glob('/spt_data/bolodata/fullrate/J*/6*/0000.g3')
srcnames_all = np.asarray([rdf.split('/')[-3] for rdf in rdfiles_2019_all])
srcnames = np.unique(srcnames_all)
srcdecs = np.asarray([np.float((ntemp.split('-')[-1])[0:2]) for ntemp in srcnames])
srcnames = (srcnames[np.argsort(srcdecs)])[::-1]

# get 2019 elnod and calibrator files
elnod_files = glob.glob(calresult_dir + 'elnod/6*.g3')
elnod_files.sort()
elnod_obsids = np.asarray([np.int((efile.split('/')[-1]).split('.g3')[0]) for efile in elnod_files])
calibrator_files = glob.glob(calresult_dir + 'calibrator/6*.g3')
calibrator_files.sort()
calibrator_obsids = np.asarray([np.int((cfile.split('/')[-1]).split('.g3')[0]) for cfile in calibrator_files])

# loop over sources, find obsids for actual source scans, then find associated elnod and calibrator files
bpfile = calresult_dir + 'boloproperties/62679602.g3'
rcw38pfile = calresult_dir + 'RCW38-pixelraster/64599638.g3'
scrfile = outdir + 'radio_day_' + str(obsid_min) + '.scr'
f1 = open(scrfile,'w')
for source in srcnames:
    srcdir = mapdir + source + '/'
    if len(glob.glob(srcdir)) == 0:
        os.mkdir(srcdir)
    for rdfile in rdfiles_2019_all:
        if source in rdfile:
            thisobsid = np.int(rdfile.split('/')[-2])
            if thisobsid > obsid_min and thisobsid < obsid_max:
                datafiles = rdfile.replace('0000','000*')
                for efile,eobsid in zip(elnod_files,elnod_obsids):
                    if eobsid < thisobsid:
                        elnod_file = efile
                for cfile,cobsid in zip(calibrator_files,calibrator_obsids):
                    if cobsid < thisobsid:
                        calibrator_file = cfile
#                thisstr = pcall_str + elnod_file + ' ' + calibrator_file + ' ' + bpfile + ' ' + rcw38pfile + ' ' + datafiles + ' -o ' + srcdir + str(thisobsid) + '.g3 -s ' + source + ' -r 0.2 -k &'
                thisstr = pcall_str + elnod_file + ' ' + calibrator_file + ' ' + bpfile + ' ' + rcw38pfile + ' ' + datafiles + ' -o ' + srcdir + str(thisobsid) + '.g3 -s ' + source + ' -k &'
                f1.write(thisstr + '\n')

f1.close()


