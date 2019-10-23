import copy
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core


#-------------------------
# Get info saved in file -
#-------------------------
rising_edge_times = []
irig_info = []
synch_clk_cnts = []
def get_irig_data(frame):
    global rising_edge_times
    global irig_info
    global synch_clk_cnts
    if frame.type == core.G3FrameType.Timepoint and 'whwp_irig_rising_edge_time' in frame.keys():
        rising_edge_times.append(int(frame['whwp_irig_rising_edge_time']))
        irig_info.append(np.array(frame['whwp_irig_info']))
        synch_clk_cnts.append(np.array(frame['whwp_irig_synch_clk_cnts']))

length_of_encoder_queue = []
encoder_cnts = []
clk_cnts = []
def get_encoder_data(frame):
    global encoder_cnts
    global clk_cnts
    if frame.type == core.G3FrameType.Timepoint and 'whwp_encoder_cnts' in frame.keys():
        #print frame['whwp_encoder_cnts']
        length_of_encoder_queue.append(len(frame['whwp_encoder_cnts']))
        for cnt in frame['whwp_encoder_cnts']:
            encoder_cnts.append(cnt)
        for cnt in frame['whwp_encoder_clk_cnts']:
            clk_cnts.append(cnt)

enc_cnts_at_ref = []
def get_ref_data(frame):
    global enc_cnts_at_ref
    if frame.type == core.G3FrameType.Timepoint and 'whwp_encoder_cnts_at_ref' in frame.keys():
        enc_cnts_at_ref.append(frame['whwp_encoder_cnts_at_ref'])

            
#-------------------
# run the pipeline -
#-------------------
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename='whwp_consumer_test.g3')
pipe.Add(get_irig_data)
pipe.Add(get_encoder_data)
pipe.Add(get_ref_data)
pipe.Run()

#----------------------------------------------
# account for counter overflow on the arduino -
#----------------------------------------------
def account_for_wrapping(input_array, nbits):
    """
    Inputs:
       input_array: array over data taking period where each entry is the encoder count
       nbits: number of bits used for storage of encoder count
    Outputs:
       out_arr: array as the same length of input_array, but wraps back to the start to
                account for reaching 2pi radians
    """
    add_val = int(2.0**nbits)
    out_arr = copy.copy(input_array)
    arr = copy.copy(input_array)
    arr[0] = 1
    arr[1:] = arr[1:] - arr[:-1]
    for ind in np.where(arr < 0)[0]:
        out_arr[ind:] += add_val
    return out_arr

def convert_to_radians(raw_counts):
    out = []
    for i in range(len(raw_counts)):
        out.append(raw_counts[i] * 2 * np.pi / 1.e5)
        if out[i] > 2*np.pi:
            out[i] = out[i] % (2*np.pi)
    return out
            
encoder_cnts_monotonic = np.array(account_for_wrapping(np.array(encoder_cnts), 16))
hwp_angle = convert_to_radians(encoder_cnts_monotonic)
enc_cnts_at_ref_monotonic = np.array(account_for_wrapping(np.array(enc_cnts_at_ref), 16))

clk_cnts = np.array(account_for_wrapping(np.array(clk_cnts), 32))
time = (clk_cnts - clk_cnts[0]) / 16.e6 # Arduino clock is 16 MHz or so


#plt.figure(1)
#plt.plot(time, encoder_cnts,'r-')
#plt.xlabel('Time [s]')
#plt.ylabel('Encoder Counts')
#plt.savefig('counts_vs_time.png')

plt.figure(1)
deriv = np.diff(encoder_cnts_monotonic)
plt.plot(time[1:], deriv,'r.')
plt.xlabel('Time [s]')
plt.ylabel('Encoder counts passed between samples')
plt.savefig('counts_passed_vs_time.png')

plt.figure(2)
plt.plot(time, hwp_angle,'b.')
plt.xlabel('Time according to the Arduino clock [s]')
plt.ylabel('HWP angle [rad]')


#plt.figure(3, figsize=(15,5))
#plt.hist(deriv, bins=1000)
#plt.xlabel('Encoder counts passed between samples')
#plt.yscale('log')
#plt.ylim(0.9,1e5)
#plt.savefig('counts_passed_hist.png')

#plt.figure(4)
#plt.hist(length_of_encoder_queue)
#plt.xlim(0,500)
#plt.xlabel('Number of encoder samples taken between bolo samples')
#plt.savefig('encoder_samples_per_bolo_sample')

plt.show()


