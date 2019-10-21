This document details how to use the tools in this folder to make a k-space weight for power spectrum estimation.

One trap -- note that in alot of places the number 18 is hard coded in. This is the number of bundles I'm using in this analysis. This is not *all* the bundles in the directory, which is why I'm making sure to only use 1-18 for consistency/speed.

First you must make two data products:

1.Take a bunch of sims (I used 100) and mock observe them, and then mock observe them with filtering. 
     The sim files can be generated using the sims in /spt/user/arahlin/synfast/ - these are sims cut to the 1500 d field for 3g, however it is bigger than the sptpol field so pointing and going through the mock-observing  pipe will make that OK. 
     To mock observe, I used the condor scripts /home/javva/spt3g_software/scratch/javva/condor/submt_maps_sims.py and /home/javva/spt3g_software/scratch/javva/condor/submit_maps_sims_noprocessing.py 
     state-of-the-art at time of writing for this analysis--
     input mock observed sims: /spt/user/javva/lowell/maps/ra0hdec-57.5/sim_input_maps_sasha_noninterp_lpf 
     filtered mock observed sims: /spt/user/javva/lowell/maps/ra0hdec-57.5/sims_maps_filtered_nonsense_p4_sasha_nointerp_npf
     Caveat - if you want to make sims that aren't TEB (i.e. only have E or T), use the make_lenspix_mapsXX.py in /home/javva/spt3g_software/scratch/javva/simulations. 

2. Take a bunch of processed maps, and put them into equal-weight bundles. 

   Best maps --
   raw maps: /spt/user/nwhitehorn/sptpol/lowell/maps
   made using: make_chronological_bundles.py
   bundle definitions:'/home/javva/spt3g_software/scratch/javva/pf_k_space_maps/bundle_defs_chron_remove_noisy_maps_stringent_cuts_poly4_correct_pickle.pkl'
   bundles:'/big_scratch/javva/filtering_bundles/chron_weights_remove_noisy_maps_stringent_cuts_poly4_correct/bundles_3gpipe_*' (used range(1,19))

3. Make 2-d fourier transforms of all of the maps
   made using: make_noise_k_space_maps.py
   current 2dffts: /spt/user/javva/lowell/2dpsd/p4_k_arrs*

4. For each bundle, for each bin in 2-d fourier space, for T,Q and U, take the inverse variance of the list of the points in that bin for all the maps in the  bundle
   made using: makevarbundles.py
   current noise bundles: /spt/user/javva/lowell/2dpsd/bun*

5. Now make the sim transfer function:
 First run make_sim_k_space_maps.py. This will make many pkl files with the 2d fourier transforms of the various sims. However, you want to smooth this into one overall transfer function. Do this usin
g make_smotth_sim_tf.py,  It saves them to a pkl file in the directory that you are in.
  current best: /home/javva/spt3g_software/scratch/javva/pf_k_space_maps/official_process_dir/simkspacetf_100_correctapod_median_pluscor.pkl

6. Correct for the transfer function
The k-space weighting effects the transfer function, so you wnt to correct yourpower spectrum with a correction funciton. Make it using make_tf_correction.py

7. Make a noise power spectrum comparison - make a coadd of all your maps (using /home/javva/spt3g_software/scratch/javva/pf_k_space_maps/code_on_deck/make_total_coadd.py or anything else, this is simple) and then run psd_comparison_plot.py

TO DO:
-evaluate the effect on sims of the scan pattern aka ratio between mock observed and input sim
-redo step 3 with correct apod and step 4 with correct 0 indexing
- evaluate leakage beam using T/E only sims and the synfast 
-averaging function in 7
- make 6 looping and 7 account for that too
