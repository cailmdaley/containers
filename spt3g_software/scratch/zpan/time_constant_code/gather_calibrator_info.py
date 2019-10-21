# This code stores all calibrator scans' frequencies, elevation, and time. 
from spt3g import core,gcp,dfmux
import pylab as pl
import numpy as np
import glob, os
import pickle as pkl

# location to save data, you should copy the default file to a location you can write
# it will only analyze data which is not in the archive
out_fname = '/home/panz/spt3g/calibrator/calibrator.pkl'
exists = os.path.isfile(out_fname)
if exists: 
    out_data = pkl.load(open(out_fname, 'rb'))
else:
    out_data = {}

# all calibrator data
full_rate_dir = '/spt/data/bolodata/fullrate/calibrator/'
obsids  = [ name for name in os.listdir(full_rate_dir)\
          if os.path.isdir(os.path.join(full_rate_dir, name)) ]
raw_data=[]
for obsid in obsids:
    raw_data.append(full_rate_dir + obsid + '/0000.g3')

# new calibrator data to be analyzed
raw_data_new = []

# data to be analyzed 
for data in sorted(raw_data):
    obsid= int(data.split('/')[-2])
    if obsid not in list(out_data.keys()) and obsid> 40000000:
        raw_data_new.append(data)

# print the length of data to be analyzed
print('length of data to be analyzed',len(raw_data_new))


class gather_cal_info(object):
        def __init__(self):
            self.calibrator = {}
            self.id='30000000'
            self.time=''
        def __call__(self, fr):
            if fr.type == core.G3FrameType.Observation:
                self.id= fr['ObservationID']
                self.time= str(fr['ObservationStart'])
                print(self.id)
            if fr.type == core.G3FrameType.Scan:
                try: 
                    minaz=min(fr['RawBoresightAz'])/core.G3Units.deg
                    len_data= len(fr['RawBoresightAz'])
                    maxaz=max(fr['RawBoresightAz'])/core.G3Units.deg
                    minel=min(fr['RawBoresightEl'])/core.G3Units.deg
                    maxel=max(fr['RawBoresightEl'])/core.G3Units.deg
                    freq=np.fft.fftfreq(len(fr['CalibratorOn']), d=1/152.6)
                    ts= fr['CalibratorOn']  
                    chopper_freq=abs(freq[np.where(abs(np.fft.fft(np.array(ts)-np.mean(ts)))\
                                 ==max(abs(np.fft.fft(np.array(ts)-np.mean(ts)))))][0])
                    print('starting time is %s'%self.time)
                    print('min and max az are %10.2f, %10.2f'%(minaz, maxaz))
                    print('min and max el are %10.2f, %10.2f'%(minel, maxel))
                    print('chopper frequency is %10.5f'%chopper_freq)
                    print('length of data is %d'%len_data)
                    self.calibrator[self.id]={}
                    self.calibrator[self.id]['time']=self.time
                    self.calibrator[self.id]['az_min']=minaz
                    self.calibrator[self.id]['az_max']=maxaz
                    self.calibrator[self.id]['el']=(minel+maxel)/2.0
                    self.calibrator[self.id]['freq']=chopper_freq
                    self.calibrator[self.id]['length']= len_data
                except:
                    print('This obsid is broken: %d'%self.id)
            if fr.type == core.G3FrameType.EndProcessing:
                out_data.update(self.calibrator)
                pkl.dump(out_data, open(out_fname,'wb'))
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename= raw_data_new)
pipe.Add(gather_cal_info)
pipe.Add(lambda fr: fr.type != core.G3FrameType.EndProcessing)
pipe.Run(profile=True)

