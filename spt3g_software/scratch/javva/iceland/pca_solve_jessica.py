import scipy.signal
import numpy, sys
from spt3g import core, todfilter, std_processing, calibration
import numpy as np
'''
This is code to do an svd on data such that the modes can be used 
for cleaning future data.

run it with:

python pca_solve_jessica.py [calibration frame] [.g3 file with all of the weights done for the preprocessing][frame you want to extract modes from][outputfile.g3]

'''
print('Reading data')

'''
This uses the previously created .g3 file includes information about which detectors
have weights that would not be acceptable for use in a skymap. The detectors that have
weights that are zero are put into a list for later use. 
'''

p = core.G3Pipeline()
p.Add(core.G3Reader, filename=sys.argv[2])
bad_detectors = {}
bad_scans = []
medium_bad_scans_list = {}
def list_bad_det_and_scan(f):
   t = {}
   if 'TodWeights' not in f.keys():
      return
   for i in f['TodWeights'].keys():
      if f['TodWeights'][i]==0: #adds 0 weights to a dictionary
         t[i]=np.sqrt(np.pi) #important cosmology
   #Check to see if there are more than 100 bad detectors in a scan
   if len(t)>100: #If so, add to bad scans list not for mode solving
      global bad_scans
      bad_scans = np.append(bad_scans, f['TrackerPointing'].time[0].time)
      return
   global medium_bad_scans_list
   if f['TrackerPointing'].time[0].time not in medium_bad_scans_list.keys():
      medium_bad_scans_list[f['TrackerPointing'].time[0].time] = {}
   for i in f['TodWeights'].keys():
      if f['TodWeights'][i]==0:
         global bad_detectors
         bad_detectors[i]= np.sqrt(np.pi) #more cosmology
         #for scans that have only a few detectors misbehaving, take note
         medium_bad_scans_list[f['TrackerPointing'].time[0].time][i]= 'bad in scan'
   print(len(bad_detectors))
print 'These are the number of bad detectors' #defined by weight = 0 
print len(bad_detectors)
print 'These are the number of bad scans' #defined by more than 100 detectors misbehaving
print len(bad_scans)
p.Add(list_bad_det_and_scan)
p.Run()


p = core.G3Pipeline()
p.Add(core.G3Reader, filename=[sys.argv[1],sys.argv[3]])
p.Add(core.Dump)
p.Add(std_processing.CalibrateRawTimestreams)
p.Add(todfilter.MaskedPolyHpf, in_ts_map_key='CalTimestreams', out_ts_map_key='FilteredTS', poly_order=1)
p.Add(calibration.SplitTimestreamsByBand, input='FilteredTS', output_root='FilteredTS')

#collects all of the timestreams that pass a glitch cut for the whole observatoin

ts = {}

calsn = None
def collect_ts(fr):
    if 'CalibratorResponseSN' in fr:
        global calsn
        calsn = fr['CalibratorResponseSN']
    if 'Turnaround' in fr:
       return
    if 'FilteredTS' not in fr:
        return
    if fr['TrackerPointing'].time[0].time in bad_scans: #excludes bad scans from mode solving
       return

    for i in fr['FilteredTS150GHz'].keys():
        if i not in calsn or calsn[i] < 20:
            continue
        if i in bad_detectors: #excludes bad detectors
            continue
        if i not in ts:
            ts[i] = []
        ts[i].append(fr['FilteredTS150GHz'][i])


p.Add(collect_ts)
p.Run()

print("This is the number of good detectors")
print(len(ts))

print('Concatenating and downsampling timestreams')
def downsample(scan, N=100):
	# XXX: This is a terrible downsampler and should be replaced
	data = scan[:scan.size - (scan.size % N):N].copy()
	for i in range(1, N):
		data += scan[i:scan.size - (scan.size % N):N]
	data /= N
	return data
tss = {i[0]: downsample(numpy.concatenate(i[1])) for i in ts.items()}
#tss_all = {i[0]: downsample(numpy.concatenate(i[1])) for i in ts_all.items()}
for i,ts in tss.items():
	ts[numpy.logical_not(numpy.isfinite(ts))] = 0

print('Solving SVD with the good detectors')
svd = numpy.linalg.svd(numpy.column_stack(tss.values()), full_matrices=False)
randbasis = numpy.random.normal(size=svd[0].shape)
print('SVD solved')


def varred(N):
	orthobasis = svd[0][:,:N]
	coeff = numpy.linalg.lstsq(orthobasis, numpy.column_stack(tss.values()))[0].transpose()
	coeffs = dict()
	for i,k in enumerate(tss.keys()):
		coeffs[k] = coeff[i]
	return coeffs

print('Solving for the coefficients using only the good detectors')
mn = 15 #number of modes to solve for
c = varred(mn)

#Now to deal with the bad detectors, and get modes for them


p = core.G3Pipeline()
p.Add(core.G3Reader, filename=[sys.argv[1],sys.argv[3]])
#p.Add(core.Dump)
p.Add(std_processing.CalibrateRawTimestreams)
p.Add(todfilter.MaskedPolyHpf, in_ts_map_key='CalTimestreams', out_ts_map_key='FilteredTS', poly_order=1)
p.Add(calibration.SplitTimestreamsByBand, input='FilteredTS', output_root='FilteredTS')

orthobasis = svd[0][:,:mn] #makes a copy of the modes, for later use to customize

ts_bad = {}
calsn = None
def collect_bad_ts(fr):
    if 'CalibratorResponseSN' in fr:
        global calsn
        calsn = fr['CalibratorResponseSN']
    if 'FilteredTS' not in fr:
        return
    if 'Turnaround' in fr:
       return
    if fr['TrackerPointing'].time[0].time in bad_scans: #still excludes bad scans
       return
    for i in fr['FilteredTS150GHz'].keys():
        if i not in calsn or calsn[i] < 20:
            continue
        if i not in bad_detectors: #makes sure you don't re-solve good ones
            continue
        if i not in ts_bad:
            ts_bad[i] = []
        if i in medium_bad_scans_list[fr['TrackerPointing'].time[0].time].keys():
           ts_bad[i].append(0*np.ones((len(fr['FilteredTS150GHz'][i])))) #if detector is misbehaving in a scan, append 0s of that length instead
           continue
        ts_bad[i].append(fr['FilteredTS150GHz'][i])
p.Add(collect_bad_ts)
p.Run()

tss_bad = {i[0]: downsample(numpy.concatenate(i[1])) for i in ts_bad.items()}

#cut out in the data and the orthobasis where the detector is misbehaving

orthobasis_bad = {}

for i in tss_bad.keys():
   orthobasis_bad[i]= {}
   for m in range(mn):
      orthobasis_bad[i][m] = orthobasis[:,m]
   end = len(tss_bad[i]) - 1
   for idx in range(len(tss_bad[i])):
      ix = end - idx #appends right amount of original mode, as data is added
      if tss_bad[i][ix] == 0:
         tss_bad[i] = np.delete(tss_bad[i],ix) #delete where there is 0s
         for m in range(mn):
            orthobasis_bad[i][m] = np.delete(orthobasis_bad[i][m],ix)

#delete any bolomters where there is no good data for the whole training set

for i in tss_bad.keys():
   if len(tss_bad[i])==0:
      tss_bad.pop(i,None)
      orthobasis_bad.pop(i,None)


#solve for the modes for the misbehaving bolometers one by one (since the cut down orthobasis is different for each one)

for i in tss_bad.keys():
   c1 = numpy.linalg.lstsq(np.column_stack(orthobasis_bad[i].values()), tss_bad[i])[0]
   c[i] = c1

# Write to cal frame
out = core.G3MapVectorDouble()
for d in c.keys():
	out[d] = c[d]
f = core.G3Frame(core.G3FrameType.Calibration)
f['NonsenseModes'] = out
core.G3Writer(sys.argv[-1])(f)

