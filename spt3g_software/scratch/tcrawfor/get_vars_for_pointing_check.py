import numpy
from spt3g import core, std_processing, gcp
#from spt3g.scratch.tcrawfor import tctools
from spt3g.util import tctools

def grab_boards(f, boardlist = []):
#    keylist = ['TrackerPointing',
#               'ACUStatus',
#               'SourceName',
#               'array',
#               'TrackerStatus',
#               'GCPFeatureBits',
    key = 'antenna0'
#    key = 'TrackerPointing'
    try:
        boardlist.append(f[key])
    except:
        pass

start_time = '24-Oct-2018:20:54:39'
stop_time = '24-Oct-2018:21:00:39'

ra_src = 201.36583 # degrees, J2000, from source.cat
dec_src = -43.018250 # degrees, J2000, from source.cat

#arcdir='/spt_data/arc/'
arcdir='/spt/data/arc/'

data1 = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_boards, boardlist=data1)
#    pipe1.Add(core.G3Writer, filename=fntemp1)
pipe1.Run()

az = []
el = []
az_exp = []
el_exp = []
az_scanoff = []
el_scanoff = []
az_horoff = []
el_horoff = []
ra_geoc_app = []
dec_geoc_app = []
ra2_geoc_app = []
utc = []
for frame in data1:
    az.append(frame['tracker']['actual'][0])
    el.append(frame['tracker']['actual'][1])
    az_exp.append(frame['tracker']['expected'][0])
    el_exp.append(frame['tracker']['expected'][1])
    az_scanoff.append(frame['tracker']['scan_off'][0])
    el_scanoff.append(frame['tracker']['scan_off'][1])
    az_horoff.append(frame['tracker']['horiz_off'][0])
    el_horoff.append(frame['tracker']['horiz_off'][1])
    ra_geoc_app.append(frame['tracker']['equat_geoc'][0])
    dec_geoc_app.append(frame['tracker']['equat_geoc'][1])
    utc.append(frame['tracker']['utc'][0])

naz = tctools.list_to_array(az)/core.G3Units.deg
nel = tctools.list_to_array(el)/core.G3Units.deg
naz_exp = tctools.list_to_array(az_exp)/core.G3Units.deg
nel_exp = tctools.list_to_array(el_exp)/core.G3Units.deg
naz_scanoff = tctools.list_to_array(az_scanoff)/core.G3Units.deg
nel_scanoff = tctools.list_to_array(el_scanoff)/core.G3Units.deg
naz_horoff = tctools.list_to_array(az_horoff)/core.G3Units.deg
nel_horoff = tctools.list_to_array(el_horoff)/core.G3Units.deg
nra_geoc_app = tctools.list_to_array(ra_geoc_app)/core.G3Units.h
ndec_geoc_app = tctools.list_to_array(dec_geoc_app)/core.G3Units.deg
utc2 = []
g3time2 = []
for ttemp in utc:
    for q in np.arange(len(ttemp)):
        utc2.append(ttemp[q].mjd)
        g3time2.append(ttemp[q])
nutc = np.asarray(utc2)

