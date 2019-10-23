from spt3g import core, mapmaker, sptpol, todfilter, dfmux, util, timestreamflagging
import spt3g.mapmaker.mapmakerutils as mmu
import spt3g.mapmaker.summingmaps as SM
import spt3g.std_processing as std_processing
from copy import copy
import numpy as np
import pickle


polarized_sum_map_fn = '/home/nlharr/tmp/cena/pol_cena_map.g3'
fit_mask_fn = '/home/nlharr/tmp/cena/fitmask.pkl'

cal_file = '/home/nlharr/tmp/cena/systeminfostarting.g3'
map_sum_file = '/home/nlharr/tmp/cena/summap.g3'

m_frame = [fr for fr in core.G3File(polarized_sum_map_fn)][2]
t = m_frame['T']
q = m_frame['Q']
u = m_frame['U']

fit_mask = copy(t)
np.asarray(fit_mask)[:] = pickle.load(open(fit_mask_fn))

class GetPolAng(object):
    def __init__(self):
        self.pol_ang = {}
        self.pol_eff = {}
        self.fit_params = {}
        self.bprops = None
        self.angle = {}
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            self.bprops = frame['BolometerProperties']
        if frame.type != core.G3FrameType.Map:
            return
        mp = frame['T'] - t
        mask = copy(fit_mask)

        np.asarray(mask)[:] = np.abs(np.sign(frame['Wunpol'].TT)) * fit_mask
        fit_params = np.array(todfilter.get_polarization_fit_params(q,u,mask,mp))
        self.fit_params[frame['Id']] = fit_params

        #import pdb, rlcompleter
        #pdb.Pdb.complete = rlcompleter.Completer(locals()).complete
        #pdb.set_trace()

        fit_params /=  np.sum( fit_params**2 )**.5
        #self.angle[frame['Id']] = np.arctan2(fit_params[1], fit_params[0]) / 2.0


        self.angle[frame['Id']] = np.arccos(fit_params[0]) / 2.0
        if fit_params[1] < 0:
            self.angle[frame['Id']] = np.pi - self.angle[frame['Id']]

pa = GetPolAng()

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = [cal_file, map_sum_file])
pipe.Add(SM.FilterMapsNotInBoloProps)
pipe.Add(SM.FilterMapToBand)
pipe.Add(pa)
pipe.Run()

import pylab as pl


if True:
    for k in pa.angle:
        pl.plot( pa.bprops[k].pol_angle, pa.angle[k] - pa.bprops[k].pol_angle, '.')

    pl.title("Cen A Polarization Fit")
    pl.xlabel("Sptpol Assumed Polarization Angle (rad)")
    pl.ylabel("Polarization Angle Inferred from Centaurus A (rad)")
    pl.show()
