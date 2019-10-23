from spt3g import core, dfmux, std_processing, util
import numpy as np
import glob, pickle
import get_temperatures
from spt3g.util import tctools

arcdir = '/spt_data/arc/'
logdir = '/spt_data/gcplogs/'
outdir = '/poleanalysis/tcrawfor/script_output/temps_temp/'

arcfiles_all1 = glob.glob(arcdir+'201812*.dat')
arcfiles_all2 = glob.glob(arcdir+'201901*.dat')
arcfiles_all = []
for af in arcfiles_all1:
    arcfiles_all.append(af)
for af in arcfiles_all2:
    arcfiles_all.append(af)
arcfiles_all.sort()
arcmjds_all = np.asarray([core.G3Time(af.split('/')[-1].split('.dat')[0]).mjd for af in arcfiles_all])

logfiles1 = glob.glob(logdir+'201812*.log')
logfiles2 = glob.glob(logdir+'201901*.log')
logfiles = []
for lf in logfiles1:
    logfiles.append(lf)
for lf in logfiles2:
    logfiles.append(lf)
logfiles.sort()

cycle_ends = []
for lf in logfiles:
    f = open(lf,'rt')
    try:
        lines = f.readlines()
        for line in lines:
            if 'cycle_tune' in line and 'Exiting' in line:
                cycle_ends.append(core.G3Time(line[0:15]))
    except:
        pass
    f.close()

data_starts = np.asarray([ttemp+12.*core.G3Units.hours for ttemp in cycle_ends])
data_ends = np.asarray([ttemp+8.*core.G3Units.hours for ttemp in data_starts])

arcfiles = []
for ds,de in zip(data_starts,data_ends):
    for amjd,af in zip(arcmjds_all,arcfiles_all):
        if amjd > ds.mjd and amjd < de.mjd:
            arcfiles.append(af)
arcfiles = np.asarray(arcfiles)
arcfiles.sort()
arcfiles = np.unique(arcfiles)

ucstage = []
icstage = []
uchead = []
ichead = []
ucpump = []
icpump = []
azs = []
els = []
ttimes = []
for file1 in arcfiles:
#for file1 in arcfiles[0:5]:
    print(file1)
    print("")
    outfile = file1.replace(arcdir,outdir)
    outfile = outfile.replace('dat','pkl')
    try:
        dtemp = get_temperatures.get_T(file1,outfile)
        dtemp=pickle.load(open(outfile,'rb'))
        ucstage.append(dtemp['ucstage'])
        icstage.append(dtemp['icstage'])
        uchead.append(dtemp['uchead'])
        ichead.append(dtemp['ichead'])
        ucpump.append(dtemp['ucpump'])
        icpump.append(dtemp['icpump'])
        azs.append(dtemp['az'])
        els.append(dtemp['el'])
        ttimes.append(dtemp['time'])
    except:
        print("Couldn't get temps from file "+file1)
    print("")

ucstage = tctools.list_to_array(ucstage)
icstage = tctools.list_to_array(icstage)
uchead = tctools.list_to_array(uchead)
ichead = tctools.list_to_array(ichead)
ucpump = tctools.list_to_array(ucpump)
icpump = tctools.list_to_array(icpump)
azs = tctools.list_to_array(azs)
els = tctools.list_to_array(els)
ttimes = tctools.list_to_array(ttimes)
mjds = np.asarray([core.G3Time(tt).mjd for tt in ttimes])


