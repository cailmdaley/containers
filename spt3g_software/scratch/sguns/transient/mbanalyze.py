from spt3g import core, calibration, mapmaker, std_processing
import numpy as np
import scipy.stats, scipy.optimize, os, sys
import argparse as ap
import fitlcmb
import pickle
import functools
import astropy.time

def obsid_to_mjd(obsid):
    return std_processing.utils.obsid_to_g3time(obsid).mjd
def mjd_to_obsid(mjd):
    return std_processing.utils.time_to_obsid(astropy.time.Time(mjd,format='mjd').iso.replace(' ','T'))
def dmjd_to_dobsid(dmjd):
    return dmjd*86400

class FitLightCurves(object):
    def __init__(self, bands, map_key, rms_key, debug=False):
        self.bands = bands if type(bands) != str else [bands]
        self.map_key = map_key
        self.rms_key = rms_key
        self.coadd_map = {}
        self.coadd_rms = {}
        self.mask = {}
        self.maps = {b: [] for b in self.bands}
        self.rms = {b: [] for b in self.bands}
        self.obs_ids = {b: [] for b in self.bands}
        self.common_obs_ids = []
        self.llhdiff = []
        self.n_maps = {b: 0 for b in self.bands}
        self.map_res = 0
        self.x_offset = 0
        self.y_offset = 0
        self.patch_shape = ()
        self.starttime = 0
        self.debug = debug

    def _makecoadd(self, band):
        core.log_fatal("This function is not yet supported")
        self.coadd_map[band] = 1./np.nansum(1./self.rms[band]**2,axis=0) * np.nansum(self.maps[band]/self.rms[band]**2,axis=0)
        if self.debug:
            print("WARNING: COADD SHENANIGANS")
            self.coadd_map[band] /= self.coadd_map[band].max()
        self.coadd_rms[band] = np.sqrt(1./np.nansum(1./self.rms[band]**2,axis=0))

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Map:
            if self.map_res == 0:
                self.map_res = frame['res']
                self.patch_shape = frame['patch shape']
                self.x_offset = frame['x offset']
                self.y_offset = frame['y offset']

            if 'Obs' in frame and frame['Obs'] == 'coadd':

                if self.coadd_map:
                    core.log_error("Found an additional coadd frame, skipping")
                    return
                
                if frame['Band'] not in self.bands:
                    core.log_warn("%s not a recognized band. Skipping"%frame['Band'])
                    return 

                band = frame['Band']
                self.coadd_map[band] = np.array(frame[self.map_key][band]).reshape(self.patch_shape)
                self.coadd_rms[band] = np.array(frame[self.rms_key][band]).reshape(self.patch_shape)
                self.coadd_sn[band] = self.coadd_map[band] / self.coadd_rms[band]
                return

            if 'Obs' in frame and frame['Obs'] != 'coadd':

                if frame['Band'] not in self.bands:
                    core.log_warn("%s not a recognized band. Skipping"%frame['Band'])
                    return 
                
                band = frame['Band']
                if not frame['Obs'] in self.obs_ids[band]:
                    self.obs_ids[band].append(frame['Obs'])
                else:
                    #print(self.obs_ids[band])
                    core.log_warn('Duplicate obsid / band: %s / %s. Skipping'%(frame['Obs'], frame['Band']))
                    return 
                self.n_maps[band] += 1
                self.maps[band].append(np.array(frame[self.map_key]).reshape(self.patch_shape))
                self.rms[band].append(np.array(frame[self.rms_key]).reshape(self.patch_shape))
                return

        if frame.type == core.G3FrameType.EndProcessing:


            
            self.obs_ids = {band: np.array(self.obs_ids[band], dtype=int) for band in self.bands}
            self.common_obs_ids = functools.reduce(np.intersect1d, [self.obs_ids[band] for band in self.bands])
            core.log_info("Analyzing lightcurves for %d maps"%self.common_obs_ids.size)
            for band in self.bands:

                common = np.isin(self.obs_ids[band], self.common_obs_ids)
                self.obs_ids[band] = self.obs_ids[band][common]
                self.maps[band] = np.array(self.maps[band])[common][self.obs_ids[band].argsort()]
                self.rms[band] = np.array(self.rms[band])[common][self.obs_ids[band].argsort()]
                if not band in self.coadd_map:
                    self._makecoadd(band)

                self.mask[band] = np.ones(self.coadd_map[band].shape)
                coadd_sn = self.coadd_map[band] / self.coadd_rms[band]
                self.mask[band][coadd_sn > 8] = 0
                for i in np.transpose(np.where(coadd_sn > 20)):
                    minind = np.ones(2,dtype=int)
                    minind[0] = max(0, i[0] - 30)
                    minind[1] = max(0, i[1] - 30)
                    self.mask[band][minind[0]:i[0]+31,minind[1]:i[1]+31] = 0 #XXX: convert this to a physical distance
                self.coadd_map[band] *= self.mask[band]


                for i in range(self.maps[band].shape[0]):
                    self.maps[band][i] = self.maps[band][i]*self.mask[band] - self.coadd_map[band] 

                self.maps[band][~np.isfinite(self.maps[band])] = 0
                self.rms[band][self.rms[band] > 6*np.median(self.rms[band][np.isfinite(self.rms[band])])] = np.infty
                
                finitemask = np.logical_and(np.isfinite(self.maps[band]), np.isfinite(self.rms[band]))
                finitemask = np.logical_and(finitemask, self.rms[band] != 0)
                finitemask = np.logical_and(finitemask, self.maps[band] != 0)
                self.rms[band][np.logical_not(finitemask)] = np.infty


            self.times = [obsid_to_mjd(coi) for coi in self.common_obs_ids]
            self.times = np.repeat(self.times, len(self.bands))
            self.starttime = self.times.min()
            self.times -= self.starttime
            
            for c, run in np.ndenumerate(np.zeros(self.patch_shape, dtype=bool)):
                lightcurve = []
                lightcurve_errors = []
                for band in self.bands:
                    run = run or np.sum(self.maps[band][:,c[0],c[1]] != 0) > 0.2*self.maps[band].shape[0]
                    lightcurve.append(self.maps[band][:,c[0],c[1]])
                    lightcurve_errors.append(self.rms[band][:,c[0],c[1]])
                if not run:
                    continue
                lightcurve = np.column_stack(lightcurve).flatten()
                lightcurve_errors = np.column_stack(lightcurve_errors).flatten()
                r = list(fitlcmb.fitlightcurve(lightcurve, lightcurve_errors, self.times, len(self.bands)))
                r[2] += self.starttime
                r[2] = mjd_to_obsid(r[2])
                r[3] = dmjd_to_dobsid(r[3])
                self.llhdiff.append((c,r))

if __name__ == "__main__":
    P = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
    P.add_argument('-d', '--delete', action='store_true', help = 'delete intermediate maps') 
    P.add_argument('-o', '--output', action='store', help = 'output file')
    P.add_argument('inputs', nargs='+', default = [], help = 'input map files')
    args = P.parse_args()

    fc = FitLightCurves(bands=['90GHz','150GHz','220GHz'], map_key = 'Tfilt', rms_key = 'Sigmafilt')
#    fc = FitLightCurves(bands=['radio','infra','gamma'], map_key = 'Tfilt', rms_key = 'Sigmafilt', debug=True)
#    print("DEBUG SESSION")

    p = core.G3Pipeline()
    p.Add(core.G3Reader, filename = args.inputs)
    p.Add(fc)
    p.Run()
    if args.delete:
        fc.maps = 0
        fc.rms = 0
        fc.coadd_map = 0
        fc.coadd_rms = 0
    pickle.dump(fc, open(args.output,'wb'))

