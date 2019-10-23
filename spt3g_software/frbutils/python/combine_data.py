from multiprocessing import Pool
split_factor = 10

import argparse

from spt3g import core, todfilter, calibration, coordinateutils, mapmaker, todfilter
from spt3g import frbutils
from spt3g.frbutils import PolyLikelihoodFiller, DeltaEventHunter, FrbEventInfo, FrbDetInfo
from spt3g.frbutils import get_ts_poly_delta_heavi_ll
from spt3g.frbutils.impulseinjections import InjectFrbSignal, FastTransientSignalBeam
from spt3g.frbutils.impulseinjections import LookForInjectedEvents
from spt3g.frbutils.frbhunting import add_frb_files

import tempfile, os

tmp_fns = []
for i in range(split_factor):
    tmp_fns.append(tempfile.NamedTemporaryFile(delete = False).name)

parser = argparse.ArgumentParser()
parser.add_argument('--out_file')
parser.add_argument('--in_files', nargs='*')
args = parser.parse_args()

#args.in_files = args.in_files[:20]

def help_f(i):
    add_frb_files(args.in_files[i::split_factor], tmp_fns[i])
p = Pool(processes = split_factor)
p.map(help_f, range(split_factor))

#    add_frb_files(args.in_files[i::split_factor], tmp_fns[i])
print("FINAL ADD")
add_frb_files(tmp_fns, args.out_file)

for tmp_fn in tmp_fns:
    os.remove(tmp_fn)


