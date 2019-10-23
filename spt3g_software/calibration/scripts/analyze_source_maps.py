#!/usr/bin/env python

import numpy as np
from spt3g import core, calibration
from spt3g.calibration import net_and_opteff_from_source_obs as net
import argparse as ap

# Usage: analyze_source_maps.py --outdir <output directory> <input obsid>
#
# calls routines in net_and_opteff_from_source_obs on maps and TOD
# from <obsid> and writes result to .g3 file. file will have name
# net_and_opteff_<obsid>_out.g3 and be located in <output directory>.

P = ap.ArgumentParser(description='Analyze source_maps data',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default='',
               help='Input obsid')
P.add_argument('-o', '--outdir', default='./', action='store',
               help='Output directory')
args = P.parse_args()

obsid = np.int(args.input[0])
outfile = args.outdir + '/net_and_opteff_' + str(obsid) + '_out.g3'
print('Output file: ' + outfile + '\n')
oframe = net.net_and_opteff_wrapper(obsid)
writer = core.G3Writer(filename = outfile)
writer(oframe)
