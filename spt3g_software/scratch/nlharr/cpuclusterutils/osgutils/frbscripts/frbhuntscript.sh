#!/bin/bash 

#So the parts you will need to customize have been noted with #** on the line above



# Printing hostname and environment 
# Good to have for debugging, in case a site causes issues 
# and needs to blacklisted or to test whether to 
# whitelist a site again. 

echo $HOSTNAME 
env 

start_dir=$PWD 
if [ "${OSG_WN_TMP}" == "" ]; 
then
    OSG_WN_TMP=$PWD 
fi 

unset PYTHONPATH
# Set up the code in the cvmfs repo
eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh` 


#sleep $(( `date +%N` % 2400 ))

set -e

# Grab the 3g code
code_dir=$PWD/code/
mkdir ${code_dir}

#** you will want to change this to your user name
wget http://stash.osgconnect.net/spt/public/nlharr/spt3g.tgz
tar xzf spt3g.tgz -C ${code_dir}

#set up the 3g environment
export SPT3G_SOFTWARE_BUILD_PATH=${code_dir}
export PATH=${code_dir}/bin:$PATH
export LD_LIBRARY_PATH=${code_dir}/bin:$LD_LIBRARY_PATH
export PYTHONPATH=${code_dir}:$PYTHONPATH

#set up the certificate we transferred over with the submit script
export X509_USER_PROXY=$PWD/user.cert 

mkdir $PWD/tmp/
work_dir=`mktemp -d --tmpdir=$PWD/tmp/` 

cd ${work_dir} 

#** Everything below here you will want to customize for your analysis

#transfer the frb side products
side_products_dir=${work_dir}/frb_data/
mkdir ${side_products_dir}

wget http://stash.osgconnect.net/spt/public/nlharr/frb_side_products.tar.gz
tar xzf frb_side_products.tar.gz -C ${side_products_dir}

cd ${code_dir}
wget http://stash.osgconnect.net/spt/public/nlharr/pywtl.tgz
tar xvzf pywtl.tgz 

#copy input data
#globus-url-copy -vb gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/nlharr/idfs/data/${1} file://${work_dir}/${1}

#globus-url-copy -b -rst gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/${1}.g3 file://${work_dir}/${1}.g3
globus-url-copy -b -rst gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/srm/spt3g/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/${1}.g3 file://${work_dir}/${1}.g3

mkdir ${work_dir}/${1}/
#globus-url-copy -b -r -rst gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/nwhitehorn/sptpol/converted/ra0hdec-57.5/${1}/ file://${work_dir}/${1}/
globus-url-copy -b -r -rst gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/srm/spt3g/user/nwhitehorn/sptpol/converted/ra0hdec-57.5/${1}/ file://${work_dir}/${1}/

#globus-url-copy -vb gsiftp://ceph-gridftp3.grid.uchicago.edu:2811/cephfs/spt/user/nlharr/idfs/data/${1} file://${work_dir}/${1}

mkdir ${work_dir}/output 

echo "SIDE"
ls ${side_products_dir}

#actually run the script
python ${start_dir}/frbhunting.py --fn_lst ${work_dir}/${1}.g3 ${work_dir}/${1}/*.g3 --tes_pos_pkl_fn ${side_products_dir}/frb_side_products/TES_positions_150ghz.pkl --xtalk_fn ${side_products_dir}/frb_side_products/xtalk-from-venus-2013 --ftsb_filename ${side_products_dir}/frb_side_products/fast_transient_signal_conversion.pkl --out_fn ${work_dir}/output/${2} --find_ll_thresh ${4} --other_ll_thresh ${5} --inject_fake_signal ${6} --fluence ${7} --time_scale=${8} --random_seed ${9} --glitchy_prob_cutoff ${10} --simulate_timestreams ${11}

#copy output data
#globus-url-copy -vb file://${work_dir}/output/${2} gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/nlharr/output/${5}/${2}
#globus-url-copy -vb file://${work_dir}/output/${3} gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/nlharr/output/${5}/${3}


globus-url-copy -vb -rst file://${work_dir}/output/${2} gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/nlharr/output/${3}/${2}
#globus-url-copy -vb file://${work_dir}/output/${4} gsiftp://ceph-gridftp3.grid.uchicago.edu:2811/cephfs/spt/user/nlharr/output/${5}/${4}

