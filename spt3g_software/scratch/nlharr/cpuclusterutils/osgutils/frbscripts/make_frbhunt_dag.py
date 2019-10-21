import argparse
from glob import glob
from string import Template
from os import path

parser = argparse.ArgumentParser()

input_directory = '/spt/user/nlharr/idfs/data/'
parser.add_argument("--output_directory_sub_name", required=True)
parser.add_argument("--job_label", required=True)
parser.add_argument("--find_ll_thresh", required=True, type = float)
parser.add_argument("--other_ll_thresh", required=True, type = float)
parser.add_argument("--inject_fake_signal", type = int)
parser.add_argument("--fluence", type = float)
parser.add_argument("--time_scale", type = float)

parser.add_argument("--glitchy_prob", type = float, default = 0.99)
parser.add_argument("--simulate_timestreams", type = int, default = 0)


args = parser.parse_args()
stemp = Template('''JOB ${job_name} /home/nlharr/spt3g_software/scratch/nlharr/cpuclusterutils/osgutils/frbscripts/frbhunt.condor
VARS ${job_name} input_file="${input_file}" out_fn="${out_fn}" output_directory="${output_sub_folder}" find_ll_thresh="${ll_thresh}" other_ll_thresh="${other_ll_thresh}" inject_fake_signal="${inject_fake_signal}" fluence="${fluence}" time_scale="${time_scale}" glitchy_prob="${glitchy_prob}" simulate_timestreams="${simulate_timestreams}"
Retry ${job_name} 5
''')

input_files = map(path.basename, sorted(glob('/spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/*.g3')))
input_files = map(lambda s: s[:s.rfind('.')], input_files)[::10]

out_string = ''

for i, in_fn in enumerate(input_files):
    base_in_name = path.basename(in_fn)
    base_name = args.job_label + "_" + base_in_name[:base_in_name.rfind('.')]

    job_name = args.job_label+'.%d'%i
    out_fn = base_name+'_processed.g3'
    
    s = stemp.substitute( input_file = base_in_name, 
                          job_name = job_name, out_fn=out_fn,
                          output_sub_folder=args.output_directory_sub_name,
                          ll_thresh = args.find_ll_thresh, other_ll_thresh = args.other_ll_thresh,
                          inject_fake_signal= args.inject_fake_signal, fluence=args.fluence, 
                          time_scale= args.time_scale, 
                          glitchy_prob=args.glitchy_prob,
                          simulate_timestreams = args.simulate_timestreams
                          )
    print s
