import argparse
from spt3g import frbutils 
from spt3g import core
#/hwm_20140318


parser = argparse.ArgumentParser()

parser.add_argument('--input_file', default = 'dummy.g3')
parser.add_argument('--output_file', default = 'dummy_filtered.g3')
parser.add_argument('--hwm_fn', default='/home/nlharr/spt_code/HardwareMaps/hwm_20131029/')
parser.add_argument('--tes_pos_pkl',
                    default='/home/nlharr/data/frb_mid/TES_positions_150ghz.pkl')
parser.add_argument('--flag_type', default = 'regular')
args = parser.parse_args()

print args

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, 
         filename = args.input_file)
if args.flag_type == 'regular':
    print args.hwm_fn
    pipe.Add(frbutils.frbhunting.GeneralFrbEventFilteringSegment,
             hwm_fn = args.hwm_fn,
             tes_pos_pkl_fn = args.tes_pos_pkl )
else:
    raise RuntimeError("unknown flag type")

pipe.Add(core.G3Writer,
         filename = args.output_file)

pipe.Run()
