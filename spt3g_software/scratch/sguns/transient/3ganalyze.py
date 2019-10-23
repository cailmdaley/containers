from spt3g import core, calibration, mapmaker
import numpy as np
import scipy.stats, scipy.optimize, os, sys
import argparse as ap
import fitlc
import pickle


class FitLightCurves(object):
    def __init__(self, band, map_key, rms_key):
        self.band = band
        self.map_key = map_key
        self.rms_key = rms_key
        self.coadd_map = None
        self.coadd_rms = None
        self.mask = None
        self.maps = []
        self.rms = []
        self.obs_ids = []
        self.llhdiff = []
        self.n_maps = 0
        self.map_res = 0
        self.x_offset = 0
        self.y_offset = 0
        self.x_size = 0
        self.y_size = 0

    def _gaussian(self, height, center_x, center_y, width_x, *args):
            """Returns a gaussian function with the given parameters"""
            width_x = float(width_x)
            width_y = float(width_x)
            rotation = 0
            if len(args) > 0:
                    width_y = float(args[0])
            if len(args) > 1:
                    rotation = float(args[1])
            return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2 + 2*rotation*(center_x-x)*(center_y-y))/2)

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Map:
            if self.band in frame['Id'] and self.map_res == 0:
                print("setting properties")
                self.map_res = frame['res']
                self.chunk_size = frame['chunk size']
                self.x_offset = frame['x start']
                self.y_offset = frame['y start']
            if self.band in frame['Id'] and frame['Obs'] == 'coadd':
                assert(frame['res'] == self.map_res)
                self.coadd_map = np.array(frame[self.map_key]).reshape(self.chunk_size)
                self.coadd_rms = np.array(frame[self.rms_key]).reshape(self.chunk_size)
                self.coadd_sn = self.coadd_map / self.coadd_rms
                self.mask = np.ones(self.coadd_map.shape)
                self.mask[self.coadd_sn > 8] = 0
                for i in np.transpose(np.where(self.coadd_sn > 20)):
                    #mask[i[0]-15:i[0]+15,i[1]-15:i[1]+15] = 0
                    minind = np.ones(2,dtype=int)
                    minind[0] = max(0, i[0] - 30)
                    minind[1] = max(0, i[1] - 30)
                    self.mask[minind[0]:i[0]+31,minind[1]:i[1]+31] = 0 #SG: made symmetric
                self.coadd_map *= self.mask 
                return
            if self.band in frame['Id']:
                self.obs_ids.append(frame['Obs'])
                self.n_maps += 1
                assert(frame['res'] == self.map_res)
                self.maps.append(self.mask*np.array(frame[self.map_key]).reshape(self.chunk_size) - self.coadd_map)
                self.rms.append(np.array(frame[self.rms_key]).reshape(self.chunk_size))
                return

        if frame.type == core.G3FrameType.EndProcessing:
            print("Analyzing lightcurves for %d maps"%self.n_maps)
            self.obs_ids = np.array(self.obs_ids, dtype=float)
            self.maps = np.array(self.maps)[self.obs_ids.argsort()]
            self.rms = np.array(self.rms)[self.obs_ids.argsort()]
            self.maps[~np.isfinite(self.maps)] = 0
            self.rms[self.rms > 6*np.median(self.rms[np.isfinite(self.rms)])] = np.infty

            finitemask = np.logical_and(np.isfinite(self.maps), np.isfinite(self.rms))
            finitemask = np.logical_and(finitemask, self.rms != 0)
            finitemask = np.logical_and(finitemask, self.maps != 0)
            self.rms[np.logical_not(finitemask)] = np.infty

            #beamsize = 1.1/2.35482/ (self.map_res/core.G3Units.arcmin) #beamwidth / fwhm to sigma / pixel res
            #maskwidth = int(np.ceil(2*beamsize))
            #template = self._gaussian(1, maskwidth, maskwidth, beamsize)(*np.indices((2*maskwidth+1, 2*maskwidth+1)))

            maskwidth = 0
            template = np.ones(1)

            times = self.obs_ids[self.obs_ids.argsort()] #SG: convert to days
            times = np.column_stack((times,)*template.size).flatten() #SG: repeat time index once for each pixel in template)
            
            coarseness = 2
            samplegrid = self.maps[:,::coarseness,::coarseness]
            print("samplegrid is size:",samplegrid.size,"shape:",samplegrid.shape)
            for c in np.column_stack(np.where(np.sum(samplegrid != 0, axis = 0) > 0.2*samplegrid.shape[0])): #SG: for c in list of coordinates where condition is met
                if c[0]*coarseness < maskwidth or c[1]*coarseness < maskwidth:
                    continue
                if c[0]*coarseness+maskwidth+1 > self.maps.shape[1] or c[1]*coarseness+maskwidth+1 > self.maps.shape[2]: #SG doublecheck def of rawmaps
                    continue
                pixmask = (slice(None), slice(c[0]*coarseness-maskwidth, c[0]*coarseness+maskwidth+1),
                                        slice(c[1]*coarseness-maskwidth, c[1]*coarseness+maskwidth+1))
                lightcurve = self.maps[pixmask]
                lightcurve_errors = self.rms[pixmask].copy()

                for k in lightcurve_errors:
                    if False and not np.isfinite(k).all():
                        k[:] = np.infty
                r = fitlc.fitlightcurve(lightcurve.flatten(), lightcurve_errors.flatten(), times, template.flatten())
                self.llhdiff.append((c,r))
                #print('LLH Difference', r[0])
                #print('Params', r[1:])

if __name__ == "__main__":
    P = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
    P.add_argument('-d', '--delete', action='store_true', help = 'delete intermediate maps') 
    P.add_argument('-o', '--output', action='store', help = 'output file')
    P.add_argument('inputs', nargs='+', default = [], help = 'input map files')
    args = P.parse_args()
    fc = FitLightCurves(band='150GHz', map_key = 'PSFiltered', rms_key = 'PSFiltered Sigma')
    p = core.G3Pipeline()
    p.Add(core.G3Reader, filename = args.inputs)
    p.Add(fc)
    p.Run()
    if args.delete:
        fc.maps = 0
        fc.rms = 0
        fc.coadd_map = 0
        fc.coadd_rms = 0
    #pickle.dump(fc, open('fc_150GHz_2018_dec_52.25_x_%d_y_%d.p'%(args.x_offset,args.y_offset),'wb'))
    pickle.dump(fc, open(args.output,'wb'))

