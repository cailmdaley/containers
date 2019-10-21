from spt3g import core, dfmux, std_processing, util
import glob, pickle
import get_temperatures
#from spt3g.util import tctools
import tctools
import numpy as np

#arcdir = '/spt_data/arc/'
arcdir = '/spt/data/arc/'
#arcfiles = []
#for i in np.arange(8):
#    thesearcfiles = glob.glob(arcdir+'20180203_1'+str(i+2)+'*.dat')
#    for taf in thesearcfiles:
#        arcfiles.append(taf)
#thesearcfiles = glob.glob(arcdir+'20180203_20*.dat')
#for taf in thesearcfiles:
#    arcfiles.append(taf)
#arcfiles = glob.glob(arcdir+'20180206_1*.dat')
#arcfiles = glob.glob(arcdir+'20180207_1*.dat')
#arcfiles = glob.glob(arcdir+'20180210_0*.dat')
#arcfiles = glob.glob(arcdir+'20180209_23*.dat')
#arcfiles = glob.glob(arcdir+'20180125*.dat')
#arcfiles = glob.glob(arcdir+'20170824*.dat')
#arcfiles = glob.glob(arcdir+'20180629*.dat')
#arcfiles = glob.glob(arcdir+'20180227*.dat')
arcfiles = glob.glob(arcdir+'20190826*.dat')
arcfiles = np.asarray(arcfiles)
arcfiles.sort()
#arcfiles = arcfiles[0:10]
#arcfiles = arcfiles[65:75]
#arcfiles = arcfiles[7:17]
#arcfiles = arcfiles[14:24]
arcfiles = arcfiles[8:20]
#outdir = '/poleanalysis/tcrawfor/script_output/temps_04feb18/'
outdir = '/big_scratch/tcrawfor/script_output/temps_04feb18/'

ucstage = []
icstage = []
uchead = []
ichead = []
ucpump = []
icpump = []
lctower = []
ttimes = []
for file1 in arcfiles:
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
        lctower.append(dtemp['lctower'])
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
lctower = tctools.list_to_array(lctower)
ttimes = tctools.list_to_array(ttimes)
mjds = np.asarray([core.G3Time(tt).mjd for tt in ttimes])


