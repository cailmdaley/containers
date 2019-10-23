from spt3g import core, dfmux
import numpy as np
import pickle

# kludge to get conversion factors
file2 = '/spt_data/bolodata/fullrate/elnod/2668507/0000.g3'
f1 = core.G3File(file2)
f1.next()
wframe = f1.next()
frame2 = f1.next()
converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
converter(wframe)
converter(frame2)
convdict = (converter.convfactors)[core.G3TimestreamUnits.Counts]

# get elnod slopes
f2 = core.G3File('/home/tcrawfor/spt_code/spt3g_software/calibration/scripts/elnod_output_2668507.g3')
eframe = f2.next()
respdict = eframe['ElnodSlopes']
respdict_w_am = {}
for key in respdict.keys():
    respdict_w_am[key] = respdict[key]*convdict[key]

# rough guess at K per airmass
t_atm = 200.
opacity = 0.05
k_am = t_atm*opacity
respdict_w_k = {}
for key in respdict.keys():
    respdict_w_k[key] = respdict_w_am[key]/k_am

