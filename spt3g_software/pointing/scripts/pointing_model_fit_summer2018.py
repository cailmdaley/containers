from __future__ import print_function
import os
import pointing_model_fit_tools as pmf

data_dir = "pointing_model_fit_data/summer2018"
spt_files = ['radio_day_summary_for_jh_26jan18.csv']
eht_files = ['pointing180128.txt',
             'pointing180128b.txt',
             'pointing180129.txt']

spt_files = [os.path.join(data_dir, f) for f in spt_files]
eht_files = [os.path.join(data_dir, f) for f in eht_files]

data, fit = pmf.run('spt_eht_summer2018_final',
                    spt_files,
                    eht_files,
                    cut_hi_el=4,
                    fixed_params=['az_tilts', 'az_encoder_offset'],
                    output_dir=data_dir,
                    process=True,
                    verbose=True,
                    plot=True)
