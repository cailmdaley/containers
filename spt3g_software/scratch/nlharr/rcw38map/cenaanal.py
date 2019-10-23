from spt3g import core, calibration, mapmaker
import numpy as np

f = core.G3File('cenamaps.g3')
bprops = None
maps = []
ids = []
i=0
for fr in f:
    print i, fr
    i+=1
    if fr.type == core.G3FrameType.Map:
        ids.append(fr["Id"])
        maps.append(fr['T'] / fr['Wunpol'].TT)
    elif fr.type == core.G3FrameType.Calibration:
        bprops = fr['BolometerProperties']

#signal_band = m[160:200, 130:170]
#noise_band = m[170:190, 165:180]

blobes = []
ulobes = []
bids = []
for i, m in enumerate(maps):
    signal_section = np.asarray(m)[160:200, 130:170]
    noise_section = np.asarray(m)[170:190, 165:180]

    noise = np.nanstd(noise_section.flatten())
    max_ind = np.where(np.abs(signal_section) == np.max(np.abs(signal_section)))
    if np.size(max_ind) != 2:
        continue
    sig = signal_section[max_ind][0]
    stn = abs(float(sig/noise))
    #import pdb; pdb.set_trace()
    print 'sig', sig, 'noise', noise, 'stn', stn
    if stn > 20:
        nmap = np.asarray(m) / sig
        blobe = np.sum(nmap[181:186, 141:145])
        ulobe = np.sum(nmap[170:174, 152:156])
        blobes.append(blobe)
        ulobes.append(ulobe)
        bids.append(ids[i])

pol_angs = map(lambda s: bprops[s].pol_angle, bids)



