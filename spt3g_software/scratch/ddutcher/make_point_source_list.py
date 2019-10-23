# Module I used to create the point source list for 1500d

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from copy import copy

sys.path.append('/home/ddutcher/code/spt3g_software/build')
from spt3g import core, coordinateutils, mapspectra, mapmaker
from spt3g.sources import source_utils
from spt3g.pointing import check_astrometry_at20 as ca

def find_and_write_ptsrc(outdir = '/home/ddutcher'):
    map_to_use = '/spt/user/ddutcher/coadds/lpf9000_tonly_0.5res.g3'
    false_detections = np.array((
        (2020, 3143),
        (2026, 3146),
        (2028, 3147),
        (2030, 3148),
        (2211, 3243),
        (2228, 3251),
        (5727, 1247)))
    
    three_bands = list(core.G3File(map_to_use))

    Tmap90 = mapmaker.mapmakerutils.remove_weight_t(three_bands[0]['T'], three_bands[0]['Wunpol'])
    apod90 = mapspectra.apodmask.make_border_apodization(
                three_bands[0]['Wunpol'],apod_type='cos',radius_arcmin=5,
                weight_threshold = 0.3)
    out_src90 = source_utils.find_sources_in_map(Tmap90, pixel_mask = apod90, nsigma=10, threshold_mjy=0,
                                                 plot_sources = False, beamsize=1.77*core.G3Units.arcmin)

    Tmap150 = mapmaker.mapmakerutils.remove_weight_t(three_bands[1]['T'], three_bands[1]['Wunpol'])
    apod150 = mapspectra.apodmask.make_border_apodization(
                three_bands[1]['Wunpol'],apod_type='cos',radius_arcmin=5,
                weight_threshold = 0.3)
    out_src150 = source_utils.find_sources_in_map(Tmap150, pixel_mask = apod150, nsigma=10, threshold_mjy=0,
                                                  plot_sources = False, beamsize=1.51*core.G3Units.arcmin)

    Tmap220 = mapmaker.mapmakerutils.remove_weight_t(three_bands[2]['T'], three_bands[2]['Wunpol'])
    apod220 = mapspectra.apodmask.make_border_apodization(
                three_bands[2]['Wunpol'],apod_type='cos',radius_arcmin=45,
                weight_threshold = 0.3)
    out_src220 = source_utils.find_sources_in_map(Tmap220, pixel_mask = apod220, nsigma=10, threshold_mjy=0,
                                                  plot_sources = False, beamsize=1.40*core.G3Units.arcmin)
      
    for band, results in [('90',out_src90),
                          ('150',out_src150),
                          ('220',out_src220)]:
        fluxes=[]
        ras, decs = [], []
        for src, dat in results.items():
            real = True
            for x,y in false_detections:
                if abs(dat['xpeak']-x)<2 and abs(dat['ypeak'] - y)<2:
                    real = False
            if not real:
                continue
            ra, dec = Tmap150.pixel_to_angle(int(round(dat['xpeak'])),
                                             int(round(dat['ypeak'])))
            ra /= core.G3Units.deg
            dec /= core.G3Units.deg
            if ra < 0:
                ra += 360.
            ras.append(ra)
            decs.append(dec)
            fluxes.append(dat['flux_mjy'])

        include_flux = True
        f = open(os.path.join(outdir,'1500d_%sGHz_10sigma_ptsrc.txt'%band),'w')
        header = "# Index\tRA\tDEC\tRadius"
        units="# \t(deg)\t(deg)\t(deg)"

        if include_flux:
            header += "\tFlux"
            units += "\t(mJy)"
        f.write(header+"\n")
        f.write(units+"\n")

        for idx, ra in enumerate(ras):
            line = "%s \t%.5f \t%.5f \t%.4f"%(
                idx+1, ra, decs[idx] , 5*core.G3Units.arcmin/core.G3Units.deg)
            if include_flux:
                line += " \t%.4f"%fluxes[idx]
            f.write(line+"\n")
        f.close()

def combine_ptsrc_lists(outdir = '/home/ddutcher'):
    dd_90 = np.loadtxt(os.path.join(outdir,'1500d_90GHz_10sigma_ptsrc.txt'))
    dd_150 = np.loadtxt(os.path.join(outdir,'1500d_150GHz_10sigma_ptsrc.txt'))
    dd_220 = np.loadtxt(os.path.join(outdir,'1500d_220GHz_10sigma_ptsrc.txt'))
    
    dict90 = {}
    dict150 = {}
    dict220 = {}
    for idx in range(len(dd_90)):#enumerate(['Index','RA','DEC', 'Radius','Flux']):
        dict90[idx] = {'ra':dd_90[idx,1],
                       'dec':dd_90[idx,2],
                       'flux':dd_90[idx,4]}

    for idx in range(len(dd_150)):#enumerate(['Index','RA','DEC', 'Radius','Flux']):
        dict150[idx] = {'ra':dd_150[idx,1],
                       'dec':dd_150[idx,2],
                       'flux':dd_150[idx,4]}

    for idx in range(len(dd_220)):#enumerate(['Index','RA','DEC', 'Radius','Flux']):
        dict220[idx] = {'ra':dd_220[idx,1],
                       'dec':dd_220[idx,2],
                       'flux':dd_220[idx,4]}
        
    three_band = {}

    ra150s = np.array([dat['ra'] for src, dat in dict150.items()])
    dec150s = np.array([dat['dec'] for src, dat in dict150.items()])

    ra220s = np.array([dat['ra'] for src, dat in dict220.items()])
    dec220s = np.array([dat['dec'] for src, dat in dict220.items()])

    count = max(dict90.keys())

    for src, dat in dict90.items():
        three_band[src] = {'ra':dat['ra'],
                           'dec':dat['dec'],
                           'flux90':dat['flux'],
                           'flux150':np.nan,
                           'flux220':np.nan}
        # match 150s
        matched = False
        idx = np.where(np.abs(ra150s - dat['ra']) < 2/60.)[0]
        if len(idx) != 0:
            for i in idx:
                if np.abs(dec150s[i] - dat['dec']) < 2/60.:
                    matched = True
                    break
        else:
            matched = False        
        if matched:
            three_band[src]['flux150'] = dict150[i]['flux']
            dict150.pop(i,None)
        else:
            three_band[src]['flux150'] = np.nan

        # match 220s
        matched = False
        idx = np.where(np.abs(ra220s - dat['ra']) < 2/60.)[0]
        if len(idx) != 0:
            for i in idx:
                if np.abs(dec220s[i] - dat['dec']) < 2/60.:
                    matched = True
                    break
        else:
            matched = False        
        if matched:
            three_band[src]['flux220'] = dict220[i]['flux']
            dict220.pop(i,None)
        else:
            three_band[src]['flux220'] = np.nan

    # Loop through remaining 150s
    for src, dat in dict150.items():
        count += 1
        three_band[count] = {'ra':dat['ra'],
                             'dec':dat['dec'],
                             'flux90':np.nan,
                             'flux150':dat['flux'],
                             'flux220':np.nan}

        # match 220s
        matched = False
        idx = np.where(np.abs(ra220s - dat['ra']) < 2/60.)[0]
        if len(idx) != 0:
            for i in idx:
                if np.abs(dec220s[i] - dat['dec']) < 2/60.:
                    matched = True
                    break
        else:
            matched = False        
        if matched:
            three_band[src]['flux220'] = dict220[i]['flux']
            dict220.pop(i,None)
        else:
            three_band[src]['flux220'] = np.nan

    # Loop through remaining 220s
    for src, dat in dict220.items():
        count += 1
        three_band[count] = {'ra':dat['ra'],
                             'dec':dat['dec'],
                             'flux90':np.nan,
                             'flux150':np.nan,
                             'flux220':dat['flux']}
    return three_band

if __name__=="__main__":
    find_and_write_ptsrc()