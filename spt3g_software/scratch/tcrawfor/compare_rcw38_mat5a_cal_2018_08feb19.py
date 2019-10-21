from spt3g import core, dfmux, std_processing, calibration   
import glob, pickle

# magic numbers
kcmb_conversion_factors = {
  'RCW38': {
    90.0*core.G3Units.GHz: 4.0549662e-07*core.G3Units.K,
    150.0*core.G3Units.GHz: 2.5601153e-07*core.G3Units.K,
    220.0*core.G3Units.GHz: 2.8025804e-07*core.G3Units.K,
  },
  'MAT5A': {
    90.0*core.G3Units.GHz: 2.5738063e-07*core.G3Units.K, # center (608, 555)
    150.0*core.G3Units.GHz: 1.7319235e-07*core.G3Units.K,
    220.0*core.G3Units.GHz: 2.145164e-07*core.G3Units.K,
  },
}

# W/K for a perfect dual-pol bolo with assumed bandwidths
kb = 1.38e-23 * core.G3Units.W / core.G3Units.K
wpk_perfect = {}
wpk_perfect[90.0*core.G3Units.GHz] = 2.*kb*28.2e9 # new from Brad
wpk_perfect[150.0*core.G3Units.GHz] = 2.*kb*42.6e9
wpk_perfect[220.0*core.G3Units.GHz] = 2.*kb*51.9e9

# rough T_RJ to T_CMB conversion
rj2cmb = {}
rj2cmb[90.0*core.G3Units.GHz] = 1.23
rj2cmb[150.0*core.G3Units.GHz] = 1.73
rj2cmb[220.0*core.G3Units.GHz] = 3.07

outd = {}
#cfdir = '/poleanalysis/sptdaq/calresult/calibration/calframe/'
cfdir = '/spt/user/production/calibration/calframe/'

rcw38files = glob.glob(cfdir + 'RCW38-pixelraster/5???????.g3')
rcw38files.sort()
rcw38files2 = glob.glob(cfdir + 'RCW38/5???????.g3')
rcw38files2.sort()
for rfile in rcw38files:
    obsid = rfile.split('/')[-1].split('.')[0] 
    if np.int(obsid) > 53000000:
        continue
    for rfile2 in rcw38files2:
        obsid2 = rfile2.split('/')[-1].split('.')[0] 
        if np.int(obsid2) > np.int(obsid):
            thisrfile2 = rfile2
            break
    try:
        cframe = core.G3File(thisrfile2).next()
    except:
        continue
    bprops = cframe['BolometerProperties']
    for bolo in bprops.keys():
        if bolo not in outd:
            outd[bolo] = {}
        if 'RCW38' not in outd[bolo]:
            outd[bolo]['RCW38'] = {}
        fband = bprops[bolo].band
        outd[bolo]['band'] = fband
        try:
            outd[bolo]['RCW38'][obsid] = -cframe['CalibratorResponse'][bolo] * cframe['RCW38FluxCalibration'][bolo] * cframe['RCW38IntegralFlux'][bolo] / kcmb_conversion_factors['RCW38'][fband] * rj2cmb[fband]
        except:
            pass

mat5afiles = glob.glob(cfdir + 'MAT5A-pixelraster/5???????.g3')
mat5afiles.sort()
mat5afiles2 = glob.glob(cfdir + 'MAT5A/5???????.g3')
mat5afiles2.sort()
for mfile in mat5afiles:
    obsid = mfile.split('/')[-1].split('.')[0] 
    if np.int(obsid) > 53000000:
        continue
    for mfile2 in mat5afiles2:
        obsid2 = mfile2.split('/')[-1].split('.')[0] 
        if np.int(obsid2) > np.int(obsid):
            thismfile2 = mfile2
            break
    try:
        cframe = core.G3File(thismfile2).next()
    except:
        continue
    bprops = cframe['BolometerProperties']
    for bolo in bprops.keys():
        if bolo not in outd:
            outd[bolo] = {}
        if 'MAT5A' not in outd[bolo]:
            outd[bolo]['MAT5A'] = {}
        fband = bprops[bolo].band
        outd[bolo]['band'] = fband
        try:
            outd[bolo]['MAT5A'][obsid] = -cframe['CalibratorResponse'][bolo] * cframe['MAT5AFluxCalibration'][bolo] * cframe['MAT5AIntegralFlux'][bolo] / kcmb_conversion_factors['MAT5A'][fband] * rj2cmb[fband]
        except:
            pass

for bolo in outd.keys():
    if 'RCW38' in outd[bolo]:
        caltemps = np.asarray(list(outd[bolo]['RCW38'].values()))
        try:
            outd[bolo]['RCW38']['median'] = np.nanmedian(caltemps)
            outd[bolo]['RCW38']['opteff'] = outd[bolo]['RCW38']['median'] / wpk_perfect[outd[bolo]['band']]
        except:
            pass
    if 'MAT5A' in outd[bolo]:
        caltemps = np.asarray(list(outd[bolo]['MAT5A'].values()))
        try:
            outd[bolo]['MAT5A']['median'] = np.nanmedian(caltemps)
            outd[bolo]['MAT5A']['opteff'] = outd[bolo]['MAT5A']['median'] / wpk_perfect[outd[bolo]['band']]
        except:
            pass

pickle.dump(outd,open('/spt/user/tcrawfor/public/cals_rcw38_mat5a_2018_08mar19.pkl','w'))

medopteffs = {}
for band in [900.,1500.,2200.]:
    medopteffs[band] = {}
    for source in ['RCW38','MAT5A']:
        oetemp = []
        for bolo in outd.keys():
            if outd[bolo]['band'] == band:
                try:
                    oetemp.append(outd[bolo][source]['opteff'])
                except:
                    pass
        medopteffs[band][source] = np.nanmedian(np.asarray(oetemp))

