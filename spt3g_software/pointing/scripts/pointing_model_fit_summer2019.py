# used to hand Radio Day and EHT pointing data to 
# fitter. 

from __future__ import print_function
import os
import pointing_model_fit_tools as pmf
import pointing_model_fit_tools_nopmodel as pmf_nopm

data_dir = "pointing_model_fit_data/summer2019"
spt_files = ['radio_day_summary_for_fitting_01feb19.csv']
#spt_files = ['radio_day_summary_for_fitting_04feb19_newboloprops_nopm.csv']
eht_files = ['pointing190130.txt']
#eht_files = []

spt_files = [os.path.join(data_dir, f) for f in spt_files]
eht_files = [os.path.join(data_dir, f) for f in eht_files]

data, fit = pmf.run('spt_summer2019_final',
#data, fit = pmf_nopm.run('spt_summer2019_fullmodel_final_nopm_noeht',
                    spt_files,
                    eht_files,
                    fixed_params=['az_tilts', 'az_encoder_offset'],
                    output_dir=data_dir,
                    process=True,
                    verbose=True,
                    plot=True)
