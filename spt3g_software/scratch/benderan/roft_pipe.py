from spt3g import core
from spt3g import dfmux
from spt3g import gcp
from spt3g import coordinateutils
import numpy as np


class SelectData():
    def __init__(self, crate = 5,slot=16):
        self.wiring_map=None
        self.crate=int(crate)
        self.slot=int(slot)
    def __call__(self,frame):
        if frame.type == core.G3FrameType.Wiring:
            self.wiring_map = frame['WiringMap']
        elif frame.type == core.G3FrameType.Scan:
            hkmap=frame['DfMuxHousekeeping']
            bolo_kys=frame['CalTimestreams_I'].keys()
            kys_select=[]
            tmp = core.G3TimestreamMap()
            for bb in bolo_kys:
                if (self.wiring_map[bb].crate_serial == self.crate) and (self.wiring_map[bb].board_slot==self.slot):
                    kys_select.append(bb)
                    tmp[bb]=frame['CalTimestreams_I'][bb]
            frame['CalTimestreams_I_subset']=tmp
            for key in frame.keys():                
                if key == 'CalTimestreams_I_subset':
                    continue
                frame.pop(key, None)
            return


if __name__=='__main__':
    import argparse as ap
    parser=ap.ArgumentParser(description =' Process the input/output files')
    parser.add_argument('input',action='store',nargs='+',help='input file list',default=[])
    parser.add_argument('-o', '--output',action='store',help='output file',default='output.g3')
    parser.add_argument('-c','--crate',  action='store',help='crate', default=5)
    parser.add_argument('-s','--slot',  action='store',help='slot', default=16)

    args=parser.parse_args()
    print 'starting pipe'

    #/poleanalysis/sptdaq/calresult/calibration/boloproperties/62679601.g3'
    #/poleanalysis/benderan (output)
    pipe=core.G3Pipeline()
    #force the first file to be the latest bolometer properties file, so that it can get the transfer function
#    pipe.Add(core.G3Reader,filename=['/poleanalysis/sptdaq/calresult/calibration/boloproperties/62679601.g3', args.input])
    pipe.Add(core.G3Reader,filename=args.input)
    print 'data read'
    #pipe.Add(tod.Downsampler,ds_factor=2)
    #print 'data downsampled'
    pipe.Add(dfmux.ConvertTimestreamUnits,Input='RawTimestreams_I', Output='CalTimestreams_I', Units=core.G3TimestreamUnits.Resistance)
    print 'converted'
    pipe.Add(SelectData,crate=args.crate,slot=args.slot)
    print 'data selected'
    pipe.Add(core.G3Writer, filename=args.output)
    print 'data written'
    pipe.Run()


'''
Syntax for calling.  Note, must include bolo properties first to do the housekeepign stuff.

python /home/benderan/spt3g_software/scratch/benderan/roft_pipe.py /poleanalysis/sptdaq/calresult/calibration/boloproperties/62679601.g3 /spt_data/bolodata/downsampled/debug-forced-scanify/62754248/0000.g3 /spt_data/bolodata/downsampled/debug-forced-scanify/62754248/0001.g3 /spt_data/bolodata/downsampled/debug-forced-scanify/62754248/0002.g3 /spt_data/bolodata/downsampled/debug-forced-scanify/62754248/0003.g3 /spt_data/bolodata/downsampled/debug-forced-scanify/62754248/0004.g3 -o /poleanalysis/benderan/roft_62754248_c5_s1.g3 -c 5 -s 1

python /home/benderan/spt3g_software/scratch/benderan/roft_pipe.py /poleanalysis/sptdaq/calresult/calibration/boloproperties/62679601.g3 /spt_data/bolodata/downsampled/userdata/62886453/0000.g3 

'''
