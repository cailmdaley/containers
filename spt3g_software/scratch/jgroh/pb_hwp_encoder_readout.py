"""
Python consumer test script that parses UDP packets from the HWP encoder Arduino Leonardo

Note that the arduino is hardcoded to require the listener to have ip address 192.168.1.42
To do this on linux:
ip a add 192.168.1.42/255.255.255.0 dev eth0

Based on code written by Nick Harrington
John Groh, August 2016
"""

import socket
import struct
import copy
import time
import numpy as np
from collections import deque
import matplotlib.pyplot as plt
import pickle

'''
I think the following comment is something to do with debugging over the serial port?

loopback command is: 
while true; do cat /dev/ttyACM0 | nc -l 9634; done;
'''

########### CAREFUL!!!!
#COUNTER_INFO_LENGTH = 100
COUNTER_INFO_LENGTH = 150
COUNTER_PACKET_SIZE = 2 + 6 * COUNTER_INFO_LENGTH
IRIG_PACKET_SIZE = 66

def de_irig(val, base_shift=0):
    """
    Helper function for parsing IRIG information
    """
    return (((val >> (0+base_shift)) & 1) + 
            ((val >> (1+base_shift)) & 1) * 2 + 
            ((val >> (2+base_shift)) & 1) * 4 + 
            ((val >> (3+base_shift)) & 1) * 8 + 
            ((val >> (5+base_shift)) & 1) * 10 + 
            ((val >> (6+base_shift)) & 1) * 20 + 
            ((val >> (7+base_shift)) & 1) * 40 )

def pretty_print_irig_info(v):
    """
    Print IRIG information in hh:mm:ss format to the console
    """
    secs = de_irig(v[0], 1)
    mins = de_irig(v[1], 0)
    hours = de_irig(v[2], 0)
    print("%d:%d:%d"%(hours, mins, secs))

class EncoderParser(object):
    _read_chunk_size = 8192
    def __init__(self, 
                 arduino_port = 8888, 
                 read_chunk_size = 4096  ):
        self.counts = 0
        self.error_queue = deque()
        self.counter_queue = deque()
        self.irig_queue = deque()

        self._is_aligned = False

        self._data = ''
        self.read_chunk_size_ = read_chunk_size

        sock_protocol = socket.SOCK_DGRAM
        self._s = socket.socket(socket.AF_INET, sock_protocol)
        self._s.bind(("", arduino_port))
        
        self.hd1index = 0
        self.trueIndex = 0
        self._index_sec = 0


    def _check_data_length(self, start_index, size_of_read):
        if start_index + size_of_read > len(self._data):
            self._data = self._data[start_index:]
            return False
        else:
            return True

    def grab_and_parse_data(self):

        # DEBUG
        time1 = time.time()
        self._data += self._s.recv(self._read_chunk_size)


        parse_index = 0
        
        prev_header = ''
        self.trueIndex = 0
        while True:
            if not self._check_data_length(parse_index, 2): 
                #print "BROKEN 0!"
                break

            #Check header for guidacne on action
            header = self._data[parse_index : parse_index + 2]
            if header == '\xaf\x1e':
                if not self._check_data_length(parse_index, COUNTER_PACKET_SIZE): 
                    print "BROKEN 1!"
                    break
                self._parse_counter_info(self._data[parse_index + 2 : parse_index + COUNTER_PACKET_SIZE ]) 
                parse_index += COUNTER_PACKET_SIZE
                self.hd1index += 1
            elif header == '\xfe\xca':
                if not self._check_data_length(parse_index, IRIG_PACKET_SIZE): 
                    print "BROKEN 2!"
                    break
                #self._parse_counter_info(self._data[parse_index + 2 : parse_index + COUNTER_PACKET_SIZE ]) 
                #parse_index += COUNTER_PACKET_SIZE
                self._parse_irig_info(self._data[parse_index + 2 : parse_index + IRIG_PACKET_SIZE ]) 
                parse_index += IRIG_PACKET_SIZE
                #rint 'DEBUG: self.hd1index is currently {:s}'.format(str(self.hd1index))
                self.hd1index = 0
                self._index_sec += 1
            elif header == '\x2a\xe1':
                print "header 3"
                if not self._check_data_length(parse_index, 4): 
                    print "BROKEN 3!"
                    break
                self._parse_error_info(self._data[parse_index + 2 : parse_index + 4 ]) 
                parse_index += 4
            elif header == '\xcd\xab':
                print "\n\n\n------\nRef Pulse Header Found\n------\n\n\n"
            else:
                print("bad header", header, prev_header)
                self._align_stream(parse_index)
                parse_index = 0
            prev_header = header
            self._data = ''

            time2 = time.time()
            diff = time2 - time1
            if diff < 0.005:
                print diff
            

            break

    def _parse_error_info(self, data):
        err_val = struct.unpack('<H', data)
        print("ERROR FOUND", err_val)
        self.error_queue.append(err_val)



    def _parse_counter_info(self, data):
        #print "!!!!!!! {:d}\t{:d}".format(struct.calcsize('<' + 'HHH'*COUNTER_INFO_LENGTH), len(data))
        derter = np.array(struct.unpack('<' + 'HHH'*COUNTER_INFO_LENGTH, data))
        self.counter_queue.append( ((derter[::3] << 16) + derter[1::3], derter[2::3]))

    def _parse_irig_info(self, data):
        unpacked_data = struct.unpack('<I' + 'H'*10 + 'I'*10, data)
        rising_edge_time = unpacked_data[0]
        irig_info = unpacked_data[1:11] # this actually encodes the time that can be converted to sec, min, hours
        #pretty_print_irig_info(irig_info)
        synch_pulse_clock_times = unpacked_data[11:21]
        self.irig_queue.append( (rising_edge_time, irig_info, synch_pulse_clock_times))

    def __del__(self):
        self._s.close()


def account_for_wrapping( input_array, nbits):
    """
    Inputs:
       input_array: array over data taking period where each entry is the encoder count
       nbits: number of bits used for storage of encoder count
    Outputs:
       out_arr: array as the same length of input_array, but wraps back to the start to account for reaching 2pi radians
    """

    add_val = int(2.0**nbits)
    out_arr = copy.copy(input_array)
    arr = copy.copy(input_array)
    arr[0] = 1
    arr[1:] = arr[1:] - arr[:-1]
    for ind in np.where(arr < 0)[0]:
        out_arr[ind:] += add_val
    return out_arr


def check_for_bad_encoder_data(enc_parser):
    enc_parser.irig_queue.clear()
    enc_parser.error_queue.clear()
    prev_v = 2**32

    notification_index = 0
    while (1):
        if len(enc_parser.counter_queue) > 0:
            vs = (enc_parser.counter_queue.popleft())[0]
            notification_index += 1
            if notification_index % 1000 == 0:
                print v
            for v in vs:
                if (v - prev_v > 4000):
                    print 'WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO BAD SHIT IS GOING DOWN'
                    print v, prev_v, v - prev_v
                prev_v = v

if __name__ == '__main__':
    """
    Test body: set up an encoder parser, take data for some time, then make some related plots
    """

    ep = EncoderParser() 

    if True:
        try:
            index = 0
            indexSec = 5
            while ep._index_sec < indexSec:
                ep.grab_and_parse_data()

        finally:
            print 'in finally'
            clk_cnts, encoder_cnts = zip(*ep.counter_queue)
            clk_cnts = np.array(clk_cnts).flatten()
            print "Length of clk_cnts:", len(clk_cnts)
            encoder_cnts = account_for_wrapping(np.array(encoder_cnts).flatten(), 16)
            print "Length of encoder_cnts:", len(encoder_cnts)

            # the clock also seems to only have 32 bits, so account for overflow on that too
            clk_cnts = account_for_wrapping(np.array(clk_cnts).flatten(), 32)
            
            rising_edge_times, irig_info, synch_clk_cnts =  zip(*ep.irig_queue)
            print "Length of rising_edge_times:", len(rising_edge_times)
            print "Length of irig_info:", len(irig_info)

            synch_clk_cnts = np.array(synch_clk_cnts).flatten()
            print "Length of synch_clk_cnts:", len(synch_clk_cnts)

            errs = zip(*ep.error_queue)
            print "Length of error queue:", len(errs)
#            print errs
            
            #Make a plot of the encoder cnts vs the clock counts
            #print len(rising_edge_times)
            #print len(encoder_cnts)
            #clkDivFactor = 139
            clkDivFactor = 1

            #Convert counts to meaningful values
            encoderTicks = 4*1e4 #Hz
            arduinoTicks = 16e6 #Hz
            hwpAngle = 2*np.pi*np.array(encoder_cnts)/float(encoderTicks)
            print hwpAngle
            hwpCos = np.cos(hwpAngle)
            timeSec = clk_cnts/float(arduinoTicks)
            print timeSec
            """
            plt.figure(0)
            plt.plot((np.array(clk_cnts) - clk_cnts[0])/float(clkDivFactor), encoder_cnts)
            plt.title("Encoder Counts vs Clock Counts")
            plt.xlabel("Clock Counts")
            plt.ylabel("Encoder Counts")
            #plt.savefig("HWP_encoderCountsVsClockCounts.png")
            """
            #Mqke a plot of the derivatives of encoder counts vs clock times
            plt.figure(1)
            deriv = np.diff(np.array(encoder_cnts))
            plt.plot(clk_cnts[1:], deriv, 'b.')
            plt.title("Derivative of Encoder Counts vs Clock Counts")
            plt.xlabel("Clock Counts")
            plt.ylabel("Derivative of Encoder Counts")
            #plt.savefig("HWP_derivativeOfEncoderCountsVSClockCounts.png")
            """
            #Make a histogram of the derivatives of encoder counts
            plt.figure(2)
            plt.hist(deriv)
            plt.title("Histogram of Derivative of Encoder Counts")
            plt.xlabel("Derivative of Encoder Counts")
            plt.ylabel("Number of Samples")
            #plt.savefig("HWP_histogramOfDerivativeOfEncoderCounts.png")
            """
            #Make a plot of the HWP angular position vs time
            plt.figure(3)
            plt.plot(timeSec, hwpCos,'b+')
            plt.title("HWP Position vs Time")
            plt.xlabel("Time [sec]")
            plt.ylabel("HWP Position")
            #plt.savefig("HWP_positionVsTime.png")

            #Write values to a text file
            #f = open("HWPPositionVsTime.txt", "w")
            #for i in range(len(timeSec)):
            #    f.write("%-20f%-20f\n" % (timeSec[i], hwpCos[i]))
            #f.close()

            # figure out if the rising edge times also wrap
            """
            plt.figure(4)
            plt.plot(rising_edge_times)
            plt.xlabel('index')
            plt.ylabel('rising_edge_times')
            """
            plt.show()

            #Write values to pickle files
            #pickle.dump(timeSec, open("timeSec_45deg_tune2.pkl", "wb"))
            #pickle.dump(hwpCos, open("hwpCos_45deg_tune2.pkl", "wb"))

    # ????
    if False:
        import threading

        def infloop():
            while (1):
                ep.grab_and_parse_data()
        t = threading.Thread(target = infloop)
        t.start()
        check_for_bad_encoder_data(ep)
