from spt3g import core
import socket, time, collections, numpy, struct, errno

def bcd10_decode(bytes):
    outbytes = [0,]*len(bytes)
    for i in range(len(bytes)):
        b = bytes[i]
        lowb = b & 0xf
        highb = b >> 5
        outbytes[i] = lowb + highb*10
    return outbytes

@core.indexmod
class CalibratorDAQ(object):
    '''
    Read calibrator data from a calibrator Arduino board.
    '''
    def __init__(self, port=42668, output='CalibratorOn', maxqueuelength=2000):
        '''
        Listen on port for UDP data from the calibrator Arduino, placing the results into
        Timepoint frames under the key given by output. If more than maxqueuelength Timepoint
        frames accumulate here with no word from the calibrator, flush the queue without
        calibrator data.
        '''
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind(('', port))
        self.sock.setblocking(False)

        self.timedcalchanges = collections.deque()
        self.calchanges = collections.deque()
        self.irigstamps = collections.deque()

        self.output = output
        self.maxqueuelength = maxqueuelength	

        self.buffer = []

    def __call__(self, frame):
        maxsize = 120

        # Use a new frame as an excuse to drain the socket queue
        while True:
            try:
                data = self.sock.recv(maxsize)
            except socket.error as e:
                if e.errno != errno.EAGAIN:
                    raise
                break

            pkttype = struct.unpack('<H', data[:2])[0]
            if pkttype == 0xbbbb: # Calibrator state change packet
                cala = struct.unpack('<QH', data[2:12])
                calb = struct.unpack('<QH', data[12:22])
                
                core.log_debug("Got cal frame ", (cala[0], cala[1]))
                core.log_debug("Got cal frame ", (calb[0], calb[1]))
                
                self.calchanges.append( (cala[0], cala[1]))
                self.calchanges.append( (calb[0], calb[1]))
            elif pkttype == 0xcafe: # IRIG packet
                

                framestarttime = struct.unpack('<Q', data[2:10])[0]
                irig_bytes = list(struct.unpack('<HHHHHHHHHH', data[10:30]))
                irig_bytes[0] = irig_bytes[0] >> 1 # Correct for bit shift in seconds field
                irig_bytes = bcd10_decode(irig_bytes)
                year = irig_bytes[5]
                yday = irig_bytes[3] + (irig_bytes[4] & 0x3)*100

                if year == 0: # Stupid IRIG generators don't always have year
                    systime = time.gmtime()
                    year = systime.tm_year - 2000 # Or centuries...
                    # Handle New Year's: if one clock is at the end of the year
                    # and the other is at the beginning, add/subtract a year
                    # appropriately. This assumes the clock deltas are at most
                    # a day, which seems fine.
                    if systime.tm_yday <= 1 and yday >= 364:
                        year -= 1
                    elif systime.tm_yday >= 364 and yday <= 1:
                        year += 1

                t = core.G3Time(s=irig_bytes[0], m=irig_bytes[1], h=irig_bytes[2], d=yday, y=year, ss=0)
                core.log_debug("Got irig frame ", t)
                subframe_starttimes = struct.unpack('<QQQQQQQQQQ', data[30:110])
                self.irigstamps.append((framestarttime, t, subframe_starttimes))
            elif pkttype == 0xe12a:
                err_code = struct.unpack('<B', data[2:3])
                core.log_error("Got error", err_code)
            else:
                core.log_warn('Unknown calibrator packet type 0x%x' % pkttype, unit='calibration.CalibratorDAQ')

        # Calibrate timing of any collected calibration packets
        self.assigntimestocaldata()

        # Now interpret the frame we were asked to read. Add to our queue,
        # and make as much progress clearing it as we can.
        self.buffer.append(frame)

        out = []
        for fr in self.buffer:
            if fr.type != core.G3FrameType.Timepoint or self.tryclear(fr):
                out.append(fr)
            else:
                break
        self.buffer = self.buffer[len(out):]

        # Warn periodically if buffer is getting long
        if len(self.buffer) > 300 and len(self.buffer) % 100 == 0:
            core.log_error('Calibrator readout frame queue at %d frames since last data. Has the calibrator encoder stopped transmitting?' % len(self.buffer))

        # Bail if the queue has gotten *really* long and we aren't making progress
        if len(self.buffer) > self.maxqueuelength and len(out) == 0:
            out = self.buffer
            self.buffer = []

        # If we got an end processing frame, we're out of luck and out of time
        if frame.type == core.G3FrameType.EndProcessing:
            out += self.buffer
            self.buffer = []

        return out

    def assigntimestocaldata(self):
        if len(self.irigstamps) < 2:
            return

        clks = [irig[0] for irig in self.irigstamps]
        times = [irig[1].time for irig in self.irigstamps]

        while len(self.calchanges) > 0:
            if self.calchanges[0][0] < clks[0]:
                # Drop archaic data
                self.calchanges.popleft()
                continue

            if self.calchanges[0][0] > clks[-1]:
                # Wait for information on data from the future
                break

            # Interpolate time otherwise
            self.timedcalchanges.append((core.G3Time(numpy.interp(self.calchanges[0][0], clks, times)), self.calchanges[0][1], self.calchanges[0][0]))
            self.calchanges.popleft()

        # Clear old IRIG calibration points
        while len(self.irigstamps) > 2: 
            if len(self.calchanges) == 0 or self.irigstamps[1][0] < self.calchanges[0][0]:
                self.irigstamps.popleft()

    def tryclear(self, frame):
        '''
        Return True if this frame has been processed successfully. False if
        we are still missing data and need to wait.
        '''
        frame_t = frame['EventHeader']

        # Bail if no data yet
        if len(self.timedcalchanges) == 0:
            return False

        if self.timedcalchanges[0][0].time > frame_t.time:
            # Not good: we have only information after this sample
            # Report state as NAN -- we'll never know
            frame[self.output] = numpy.nan
            return True
        
        # Drop brackets before this frame: time won't go backward
        while len(self.timedcalchanges) >= 2 and self.timedcalchanges[1][0].time < frame_t.time:
            self.timedcalchanges.popleft()

        # At this point, self.timedcalchanges[0] is guaranteed before or at
        # the frame time by the second if block and self.timedcalchanges[1]
        # is guaranteed greater by the while loop -- or we are short on data.
        if len(self.timedcalchanges) < 2:
            return False

        # Now the guarantee is hard and the state at this frame is what it
        # moved to in self.timedcalchanges[0]
        frame[self.output] = float(self.timedcalchanges[0][1])
        return True

