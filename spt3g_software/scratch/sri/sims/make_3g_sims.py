#!/usr/bin/env python
import numpy as np, os, sys, glob, argparse

#pending - remove this later
##sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/2018_10/spt3g_sims/')

#load class now
import sims_master
sims = sims_master.sims_master()

########################################################################
"""
example run:
python make_simple_sims.py -paramfile params_for_sims.txt -op_folder None -sim_index 1 -make_cmb_skies 1 -make_foreground_skies 1
will make the simulation with sim_index = 1 and dump the CMB and foreground alms
"""
########################################################################

"""
#define argparse stuffs here
"""

parser = argparse.ArgumentParser(description='')
parser.add_argument('--paramfile', action='store', type=str, required=True,
	help='paramfile - file containing sim parameters')
parser.add_argument('--op-folder', action='store', type=str, required=True,
	help='op_folder: where all the sims products will be stored. make sure you have permissions to write in the folder')
parser.add_argument('--sim-index', action='store', type=int, default=-1,
	help='sim_index. if -1, then this code will make N=noofsims simulations. Else this code will create a specific CMB sky')
parser.add_argument('--noofsims', action='store', type=int, default=2,
	help='noofsims')
parser.add_argument('--nu', action='store', type=int, default=-1,
	help='frequency channel. important for foregrounds and instrumental beam/noise. -1 generates all frequencies [90, 150, 220]')
parser.add_argument('--random-seed', action='store', type=int, default=-1,
	help='random_seed for a specific. if -1, then the code will use create a random seed based on sim_index')
parser.add_argument('--make-cmb-skies', action='store', type=int, default=1,
	help='make_cmb_skies. if 0, then we are either making foreground or noise sims')
parser.add_argument('--make-foreground-skies', type=int, default=1,
	help='make_foreground_skies. if 0, then we are either making CMB or noise sims')

args = parser.parse_args()
'''
### this piece ensures argprase variables can be accessed as variable_name rather than args.variable_name
args_keys = args.__dict__
for kargs in args_keys:
	param_value = args_keys[kargs]
	if isinstance(param_value, str):
		cmd = '%s = "%s"' %(kargs, param_value)
	else:
		cmd = '%s = %s' %(kargs, param_value)
	exec(cmd)
'''
########################################################################
#read param file now and store it in a dic
print(sims.separator)
print('\n\tread param file now and store it in a dic')
#copy_to_class = 1 ensures that the dic is stored in sims_master.sims_master() class
param_dict = sims.fn_get_param_dict_from_paramfile(args.paramfile, copy_to_class = 1)
sims.param_dict = param_dict
##print(param_dict)

########################################################################
#set up sims. uses the param_dict (like executing / reading CAMB outputs, etc.)
print(sims.separator)
print('\n\tset up sims')
sims.fn_sim_setup(param_dict)
########################################################################
#get the specified frequencies
if args.nu == -1:
	nu_list = [95, 150, 220]
else:
	nu_list = [args.nu]

########################################################################
#now make CMB sims

if args.make_cmb_skies:
	print(sims.separator)
	print('\n\tcreating CMB skies')
	CMB_SIM_MAPS = sims.fn_make_Gaussian_simulated_skies(param_dict, args.sim_index, nu = args.nu, noofsims = args.noofsims, random_seed = args.random_seed)
	##print CMB_SIM_MAPS.shape##;sys.exit()
	subfolder = 'CMB'
	sims.fn_dump_files(CMB_SIM_MAPS, args.op_folder, subfolder)

########################################################################

#make foregorund sims
if args.make_foreground_skies:

	#now make foregrounds sims
	for nu in nu_list:

		print(sims.separator)
		print('\n\tcreating Gaussian foregrounds skies for %sGHz channel' %nu)
		FG_SIM_MAPS_DIC = sims.fn_make_Gaussian_simulated_skies(param_dict, args.sim_index, nu = nu, cmb_sky = 0, fg_sky = 1, noofsims = args.noofsims, random_seed = args.random_seed)
		##print FG_SIM_MAPS_DIC.keys()
		
		for key in FG_SIM_MAPS_DIC:
			sims.fn_dump_files(CMB_SIM_MAPS, args.op_folder, key)

		#now make Poisson realisations of point sources
		if (1):##not param_dict['Gaussian_psources']:

			"""
			#radio point sources 
			if (0):##param_dict['make_RG']:
				RG_SIMS = sims.fn_make_Poisson_source_realisation('RG', param_dict, nu, sim_index, noofsims = noofsims, random_seed = random_seed)
			"""

			#dusty point sources 
			if param_dict['make_DG']:

				print('\n\tcreating dusty Poisson sims of bright sources for %sGHz channel' %nu)
				DG_POISSON_SIMS = sims.fn_make_Poisson_source_realisation('DG', param_dict, nu, args.sim_index, noofsims = args.noofsims, random_seed = args.random_seed)
				##print DG_POISSON_SIMS.shape
				subfolder = 'DG'
				pref = 'poisson_sim'
				sims.fn_dump_files(CMB_SIM_MAPS, args.op_folder, subfolder, sim_pref = pref)

########################################################################
print(sims.separator)
