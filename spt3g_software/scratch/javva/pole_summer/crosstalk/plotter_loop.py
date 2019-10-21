import os
import glob
from spt3g.cluster.condor_tools import condor_submit
import argparse
from spt3g import core
import pickle 
#import make_rcw38_template

parser = argparse.ArgumentParser()
#parser.add_argument('obsids', nargs = '+', type = str)
parser.add_argument('-s','--source', action ='store', 
                    default='ra0hdec-57.5')
obsnum = '63643116'
pargs = parser.parse_args()
#obsids = pargs.obsids

main_bolos = pickle.load(open('/home/javva/spt3g_software/scratch/javva/pole_summer/crosstalk/allbolos.pkl', 'rb'))

job_root = 'xtalk_one_map_lite_subset_abs_%s/*'%obsnum
out_root = os.path.join('/spt/user/javva/crosstalk/rcw38_singlebolo_coadds_analysis', job_root)
done_files = glob.glob(out_root)

def make_plots(obsnum, df1):
    main_b = df1.split('_')[-1].split('.pkl')[0]
    os.system('python check_rcw38_answer.py /spt/user/javva/crosstalk/templates/%s.pkl /spt/user/javva/crosstalk/ind_maps/%s_%s.pkl %s'%(obsnum,obsnum, main_b, df1))


for df in done_files:
    print(df)
    make_plots(obsnum, df)
