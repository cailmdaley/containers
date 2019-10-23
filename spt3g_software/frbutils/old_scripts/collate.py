import argparse
from spt3g import core, frbutils 

parser = argparse.ArgumentParser()
parser.add_argument('--input_files', nargs='*')
parser.add_argument('--output_file')

args = parser.parse_args()

core.G3Logger.global_logger.set_level(core.G3LogLevel.LOG_INFO)
frbutils.frbhunting.collate_g3_files(in_file_list = args.input_files,
                                    out_file = args.output_file)

