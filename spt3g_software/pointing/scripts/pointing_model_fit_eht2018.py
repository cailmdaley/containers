from __future__ import print_function
import os
import pointing_model_fit_tools as pmf

data_dir = "pointing_model_fit_data/eht2018"
eht_files = [
    'pointing180418.txt',
    'pointing180419.txt',
    'pointing180421b.txt',
    'pointing180422.txt',
]

eht_files = [os.path.join(data_dir, f) for f in eht_files]

common_opts = dict(
    spt_files=[],
    el_range = [0, 65],
    fit_params={
        'a2': -0.00915776815744, # az tilts measured 04/24
        'a3': 0.00285415390982,
        'a4': -0.022445333333333334, # el tilt, summer2018
        'a0': 0.06319997222222222, # flexure terms, summer 2018
        'a1': 0.05921644444444445
    },
    fixed_params=['az_tilts', 'el_tilt', 'az_encoder_offset'],
    output_dir=data_dir,
    process=True,
    verbose=True,
    plot=True)

data, fit = pmf.run('eht_201804_all', eht_files=eht_files, **common_opts)
