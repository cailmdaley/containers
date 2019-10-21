import argparse
from spt3g import core, frbutils 

parser = argparse.ArgumentParser()
parser.add_argument('--input_files', nargs='*')
parser.add_argument('--output_file')
parser.add_argument('--ev_thresh', type=float, default = 10.0)
parser.add_argument('--other_thresh', type=float, default = 8.0)

args = parser.parse_args()

core.G3Logger.global_logger.set_level(core.G3LogLevel.LOG_INFO)
frbutils.frbhunting.make_other_thresholds(in_file_list = args.input_files,
                                          out_file = args.output_file,
                                          event_thresh = args.ev_thresh,
                                          other_thresh = args.other_thresh)
