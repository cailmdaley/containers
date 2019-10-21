import tarfile
import numpy
from glob import glob
import argparse as ap
import os.path
P = ap.ArgumentParser(description='Taring files for condor',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)


P.add_argument('source', action='store', 
           default=None, help='name of source')
P.add_argument('bp', action='store', 
           default=None, help='True False, whether to tar bolometer property files')
P.add_argument('-o','--observation', action='store', default=None, help='observation ID')
P.add_argument('-d','--daterange', action='store', default=None, help='date range for timeseries plot like 2017-08')


args = P.parse_args()
obstype = str(args.source)
tar = tarfile.open("/scratch/dyang18/" +obstype +".tar.gz", "w:gz")
root_folder ='/spt/user/production/calibration/' 
input_file_lst = glob(root_folder +obstype+'/'+'*')

for fil in input_file_lst:
    tar.add(fil,arcname=fil.split('/')[-1],recursive =False)
    
tar.close()



if args.bp:
    print("Taring the corresponding bolometer properties")
    tar = tarfile.open("/scratch/dyang18/" +obstype +"BolometerProperties" +".tar.gz", "w:gz")
    root_folder ='/spt/data/bolodata/downsampled/'+obstype+'/'
    input_file_lst = glob(root_folder+'*')
    for fil in input_file_lst:
        obsid= int((fil.split('/'))[-1])
        if os.path.exists(fil+'/nominal_online_cal.g3'):
            tar.add(fil+'/nominal_online_cal.g3',arcname=str(obsid)+'BolometerProperties.g3')


tar.close()
