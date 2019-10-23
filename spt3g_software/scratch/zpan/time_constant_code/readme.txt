This folder includes all the codes for fitting time constants.

Run gather_calibrator_info.py first
You will need to define the location to save the calibrator info,
including the frequency, elevation, azimuth, and time.
This code doesn't overwrite existing info if a file already exists.
It will update the existing file with the missing info. 
Running this code takes a while, so you might start from
/home/panz/spt3g/calibrator/calibrator.pkl

Then run find_and_analyze_calsweeps.py
You can customize the obsid range you want to analyze. 
It will find all cal sweeps within that obsid range. 
You need to define the location of calibrator information
generated from gather_calibrator.info.py.
Define a folder for saving calibrator analysis output files
and another folder for saving time constants.
Fit plots and tau histograms will also be saved in the second folder.
This code calls tau_analysis.py, which requires a 
physical name to logical name matching file (for plotting). 
This file exists in /home/panz/useful_data/
The 2019 file is logical_to_physical_name_19.pkl
and the 2018 file is logical_to_physical_name_post_event.pkl
