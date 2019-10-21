from spt3g import core
from spt3g import dfmux
import numpy as np


class CalAnalysis():
    def __init__(self):
        self.wiring_map = None
        return

    def __call__(self,frame):
        if frame.type == core.G3FrameType.Wiring:
            self.wiring_map=frame['WiringMap']
        elif frame.type == core.G3FrameType.Scan:
            hkmap=frame['DfMuxHousekeeping']
            kys=frame['RawTimestreams_I'].keys()
        #amplitudes=np.zeros(len(kys))
            tmp=core.G3MapDouble()
            for ind,ky in enumerate(kys):
                if dfmux.HousekeepingForBolo(hkmap,self.wiring_map,ky).carrier_amplitude >0 and dfmux.HousekeepingForBolo(hkmap,self.wiring_map,ky).state=='tuned':
                    sample_rate=frame['RawTimestreams_I'][ky].sample_rate/core.G3Units.Hz   #this is in Hz now
                    nbins=len(frame['RawTimestreams_I'][ky])
                    window=np.hanning(nbins)
                    fft_tmp=np.fft.fft(frame['RawTimestreams_I'][ky]*window)
                    freq=np.array(range(0,(nbins-1)/2+2))*(sample_rate/nbins)
                    psd_tmp=np.zeros(len(freq))
                    psd_tmp[0]=np.real(fft_tmp[0]*np.conj(fft_tmp[0]))
                    psd_tmp[nbins/2]=np.real(fft_tmp[nbins/2]*np.conj(fft_tmp[nbins/2]))
                    psd_tmp[1:nbins/2]=2*np.real(fft_tmp[1:nbins/2]*np.conj(fft_tmp[1:nbins/2]))
                    gind=np.where(np.abs(freq-6) == np.min(np.abs(freq-6)))[0]
            #amplitudes[ind]=psd_tmp[gind]
                    tmp[ky]=psd_tmp[gind[0]]
                #        frame['amplitudes']=tmp
            #frame
                else: 
                    tmp[ky]=0
            frame['cal_response_amy']=tmp
            return 



if __name__=='__main__':
    import argparse as ap
    parser=ap.ArgumentParser(description =' Process the input/output files')
    parser.add_argument('input',action='store',help='input file')
    parser.add_argument('output',action='store',help='output file')
    args=parser.parse_args()
    print 'starting pipe'
    pipe=core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=args.input)
    print 'data read'
    pipe.Add(CalAnalysis)
    print 'cal anal performed'
    pipe.Add(core.G3Writer, filename=args.output)
    print 'data wrote out'
    pipe.Run()

'''
outframe=core.G3File('/home/benderan/sandbox/tmp_frame_output.g3')
tmp_data=outframe.next()
tmp_data=outframe.next()
tmp_data=outframe.next()
kys=tmp_data['cal_response_amy'].keys()
cal_response=np.zeros(len(kys))
for ind,ky in enumerate(kys):
    cal_response[ind]=tmp_data['cal_response_amy'][ky]
'''

    
