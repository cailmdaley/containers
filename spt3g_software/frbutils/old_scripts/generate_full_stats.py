import argparse, pickle, glob
from spt3g import core, frbutils 

parser = argparse.ArgumentParser()

parser.add_argument('--frb_keys', nargs = '*', required = True)
parser.add_argument('--input_files', nargs = '*', required = True)
parser.add_argument('--output_file', required = True)

args = parser.parse_args()

histogram_bins = range(165)

object_dic = {}
for k in args.frb_keys:
    object_dic[k] = (   
        frbutils.frbanalysis.LookForInjectedEvents(frb_event_key = k),

    )

input_files = []
for fn in args.input_files:
    input_files += glob.glob(fn)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename = input_files)
for k in args.frb_keys:
    pipe.Add(object_dic[k][0])
pipe.Run()


for k in object_dic:
    print k, object_dic[k][0].n_events_found,  object_dic[k][0].n_scans_looked_through
