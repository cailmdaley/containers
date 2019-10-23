from spt3g import core, calibration,std_processing
import matplotlib.pyplot as plt
import numpy as np
import os, sys

# IMPORTANT: these tests are supposed to be carried out on calibrator and elnod data that was taken close together. Use nearly adjacent observation ids, or check the autoprocessing logs to see which elnod a calibration stare was originally processed with. #

def nanequals(a,b):
    if np.isnan(a) and np.isnan(b):
        return True
    else:
        return a==b

def compare_data_rotations_from_elnod_cal(calfile, elnodfile, some_data, bp = '/poleanalysis/sptdaq/calresult/calibration/boloproperties/62679602.g3', 
                                          n_scan_frames = 2, plot = True, n_bolos_to_plot = 100):
    #FrameGrabber stolen from /scratch/nlharr/secretcode.py 
    class FrameGrabber(object):
        def __init__(self):
            self.frames = []
        def __call__(self, frame):
            if frame.type != core.G3FrameType.EndProcessing:
                self.frames.append(frame)

    fg = FrameGrabber()
    p = core.G3Pipeline()
    p.Add(core.G3Reader, filename = [bp, calfile, elnodfile, some_data], n_frames_to_read = 5 + n_scan_frames)
    p.Add(std_processing.CalibrateRawTimestreams, units = core.G3TimestreamUnits.Power, output = 'RotatedWithElnod', keep_original = True, calibrator_iq_rotation = False)
    p.Add(std_processing.CalibrateRawTimestreams, units = core.G3TimestreamUnits.Power, output = 'RotatedWithCalibrator', keep_original = True, calibrator_iq_rotation = False)
    p.Add(lambda fr: fr.type == core.G3FrameType.Scan)
    p.Add(fg)
    p.Run() 
    
    tsmap_calib = core.G3TimestreamMap.concatenate([frame['RotatedWithCalibrator'] for frame in fg.frames])
    tsmap_elnod = core.G3TimestreamMap.concatenate([frame['RotatedWithElnod'] for frame in fg.frames])
    print("Bolometer keys in calib and elnod timestream maps are the same: %r"%(set(tsmap_calib.keys()) == set(tsmap_elnod.keys()))) 
    tsmap_sub = core.G3TimestreamMap()
    for bolo in tsmap_calib.keys():
        ts_calib = tsmap_calib[bolo]
        ts_elnod = tsmap_elnod[bolo]
        tsmap_sub[bolo] = ts_calib - ts_elnod
        # Haven't thought of a good way to condense the difference into a single number yet (could just do np.mean(..) but I don't care about DC offsets between the two)
        # So just go ahead and plot some bolos and convince yourself the difference is 0 everywhere.
    if plot:
        fig,ax = plt.subplots(1,1)
        n = 0
        while (n < n_bolos_to_plot) and len(tsmap_sub.keys()) > 0:
            b = np.random.choice(tsmap_sub.keys())
            if not np.isnan(tsmap_sub[b]).any():
                ax.plot(tsmap_sub[b], label=b)
                n += 1
            tsmap_sub.pop(b, 0)
        ax.set_title('%d differenced elnod-rotated and calib-rotated timestreams'%n)
        ax.set_ylabel('Watts')
        #ax.legend()
        fig.show()
    

def compare_elnod_cal_phase_angles(calfile, elnodfile, plot = True, cal_SN_cut = False):

    cal = core.G3File(calfile).next()
    eln = core.G3File(elnodfile).next()

    cvec_keys = ['CalibratorResponseI', 'CalibratorResponseQ']
    evec_keys = ['ElnodEigenvalueDominantVectorI', 'ElnodEigenvalueDominantVectorQ'] 
   
    print("Identical bolometer keys in cal and elnod angles: %r"%(set(cal[cvec_keys[0]].keys()) == set(eln[evec_keys[0]].keys())))
    asym_nans = 0
    diffs,c_angles,e_angles = [],[],[]
    for b in set(cal[cvec_keys[0]].keys()) & set(eln[evec_keys[0]].keys()):
        c0,c1 = cal[cvec_keys[0]][b], cal[cvec_keys[1]][b]
        e0,e1 = eln[evec_keys[0]][b], eln[evec_keys[1]][b]
        c_angle = np.arctan2(c1,c0)
        e_angle = np.arctan2(e1,e0)
        if np.isnan(c_angle) != np.isnan(e_angle):
            asym_nans += 1
        passed_cal_cut = not (cal_SN_cut and cal['CalibratorResponseSN'][b] < 20.)
        if not (np.isnan(c_angle) or np.isnan(e_angle)) and passed_cal_cut:
            c_angle = c_angle - 2*np.pi*core.G3Units.rad if c_angle > np.pi*core.G3Units.rad else c_angle
            e_angle = e_angle - 2*np.pi*core.G3Units.rad if e_angle > np.pi*core.G3Units.rad else e_angle
            diffs.append(np.abs(c_angle - e_angle))
            c_angles.append(c_angle)
            e_angles.append(e_angle)
    print("Found %d bolos that have nan in one frame but not the other"%(asym_nans))
    dm = np.mean(diffs)
    print("average (absolute) difference in phase angle between cal en elnod phasing: %.3f degrees (and cos(%.1f) = %.2f)"%(dm/core.G3Units.deg,dm/core.G3Units.deg,np.cos(dm/core.G3Units.rad)))
    if plot:
        fig,ax = plt.subplots(1,1)
        bins = np.linspace(-180,180,24)
        ax.hist(np.array(c_angles)/core.G3Units.deg, color='blue', alpha=0.5, label='Calibrator Angles',bins=bins)
        ax.hist(np.array(e_angles)/core.G3Units.deg, color='red', alpha=0.5, label='Elnod Angles',bins=bins)
        ax.set_title('Elnod and Calibrator Phase Angles')
        ax.legend()
        fig.show()
    


def do_two_calframes_look_similar(file1, file2, plot = True):

    cal1 = core.G3File(file1).next()
    cal2 = core.G3File(file2).next()

    phasediff = cal1['CalibratorResponsePhase'] - cal2['CalibratorResponsePhase']
    freqdiff = cal1['CalibratorResponseFrequency'] - cal2['CalibratorResponseFrequency']
    print("Phase difference: %d samples"%(phasediff))
    print("Frequency difference: %.9f Hz "%(phasediff/core.G3Units.Hz))

    mapkeys = ['CalibratorResponse', 'CalibratorResponseSN']
    if plot:
        fig, ax = plt.subplots(1, len(mapkeys))
        nplot = 0
    for k in mapkeys:
        asym_nans = 0
        diffs,vals1,vals2 = [],[],[]
        print("Identical bolometer keys in %s: %r"%(k, set(cal1[k].keys()) == set(cal2[k].keys())))
        for b in set(cal1[k].keys()) & set(cal2[k].keys()):
            if not (np.isnan(cal1[k][b]) or np.isnan(cal2[k][b])):        
                diffs.append(cal1[k][b] - cal2[k][b])
                vals1.append(cal1[k][b])
                vals2.append(cal2[k][b])
            if np.isnan(cal1[k][b]) != np.isnan(cal2[k][b]):
                asym_nans += 0
        if k == 'CalibratorResponse':
            u = 1e-6 * core.G3Units.pW
            U = 'attoWatts'
        else:
            u = 1.
            U = ''
        print("Average difference in %s: %.9f %s"%(k, np.mean(diffs)/u, U))
        print("Found %d bolos that have nan for %s in one frame but not the other"%(asym_nans, k))
        if plot:
            binmin = min(min(vals1),min(vals2))/u
            binmax = max(max(vals1),max(vals2))/u
            bins = np.linspace(binmin, binmax, 21)
            ax[nplot].hist(np.array(vals1)/u, color='blue', alpha=0.5, label = os.path.basename(file1), bins=bins) 
            ax[nplot].hist(np.array(vals2)/u, color='red', alpha=0.5, label = os.path.basename(file2), bins=bins)
            ax[nplot].set_title(k)
            ax[nplot].set_xlabel(U)
            ax[nplot].legend()
            nplot += 1
            fig.show()

def are_two_calframes_identical(file1, file2):

    cal1 = core.G3File(file1).next()
    cal2 = core.G3File(file2).next()

    keys = ['CalibratorResponsePhase', 'CalibratorResponseFrequency'] 
    mapkeys = ['CalibratorResponse', 'CalibratorResponseSN']

    same = True

    for k in keys:
        if k in cal1 and k in cal2:
            same = same and nanequals(cal1[k],cal2[k])
    for k in mapkeys:
        if k in cal1 and k in cal2:
            same = same and (set(cal1[k].keys()) == set(cal2[k].keys()))
            for b in cal1[k].keys():
                same = same and nanequals(cal1[k][b],cal1[k][b])
    if same:
        print("They're the same! Your code has no obvious bugs! Hurray!")
        return True
    else:
        print("They're different. Your code sucks. Back to the drawing board!")
        return False
