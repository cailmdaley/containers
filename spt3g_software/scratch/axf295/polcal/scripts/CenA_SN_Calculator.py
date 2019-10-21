import argparse as ap
from spt3g import core
## when running scripts in same path (or else merged into spt3g_software)
#import sys
#sys.path.append('./spt3g_software/polcal/python/')
import CenA_Map_Utils as CMU
import General_Utils  as GU
import os

# Usage: CenA_SN_Calculator <input files.g3> -o /path/to/outputdir/ -m /path/to/sourcemask/ -s CenA
P = ap.ArgumentParser(description='Calc Signal to Noise of Single Bolo Source Maps',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)

P.add_argument('-obs','--ObsIDs', action='store', nargs='+', default=[],
help='CenA Observation IDs')

## Static locations added for testing, see below
#P.add_argument('input_files', action='store', nargs='+', default=[],
#help='Input files; calframe, indbolomaps')

P.add_argument('-o', '--output_dir', action='store', default='./', help='Output Directory')

P.add_argument('-m', '--mask', action='store', 
default='/home/axf295/2019/code/spt3g_software/scratch/axf295/polcal/Masks/0p25arcmin/', help='source mask location')

P.add_argument('-s', '--source', action='store', default='CenA', help='name of source')

P.add_argument('-P', '--PLOT', action = 'store', default=False, help='Plotting switch')

P.add_argument('-d', '--data_file', action='store', default='./singlebolomaps.g3', help='Location of SingleBolomaps')

P.add_argument('-c', '--cal_file', action='store', default='./cal.g3', help='Location of CalFrame')

## Static data locations on scott
#singlebolomap_loc = '/spt/user/production/calibration/CenA-pixelraster/singlebolomaps/'
#calframe_loc   = '/spt/user/production/calibration/calframe/CenA-pixelraster/'

args = P.parse_args()

source = args.source
obsids = args.ObsIDs 
PLOT   = args.PLOT

singlebolomap_loc = args.data_file
calframe_loc      = args.cal_file


## Function to get all obsIDs from 2018-2019
if obsids == [] or obsids == 'all':
    obsids = GU.get_cena_obsids(years=[2019])

        
for obs in obsids:
    ##if static_loc used, need to include obs here.
    cal_file = calframe_loc# + obs+'.g3'
    indbolomaps = singlebolomap_loc# +obs+'.g3'
    
    if not os.path.exists(cal_file):
        print('Cal File, %s, couldn\'t be found'%cal_file)
        continue
    if not os.path.exists(indbolomaps):
        print('Individual Bolo Map File, %s, couldn\'t be found'%indbolomaps)
        continue
    GU.make_directory(args.output_dir+'%s/'%obs)
    output_dir = args.output_dir + '%s/'%obs
    input_file_list = [cal_file] + [indbolomaps]

    log = open(args.output_dir+'%s/%s_logfile.txt'%(obs,obs),'w')

    log.write('Running Analysis Script on observation %s\n'%obs)
    for key in vars(args):
        log.write(key+' : '+str(vars(args)[key])+'\n')
    
    if args.mask != None:
        log.write('Calculating Signal to Noise Using Masks\n')
        CMU.calc_cena_sn_using_mask(input_file_list,obs,source,log=log,mask_location=args.mask,
                                    output_dir=output_dir,PLOT=PLOT)
        log.write('Done SN Calc of Masked Region')
    else:
        log.write('Calculating Signal to Noise\n')
        CMU.calc_signal_to_noise(output_dir,obs,log,source=Source,Save_pngs = 0)
        log.write('Done SN Calc Using Central Region')
    log.close()

    
