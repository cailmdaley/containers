pro breakup_2008_wmap7ps

filein = '/data11/cr/hrsim_spt_run3_2009_map2008_150_gauss.dat'
fileout = '/data11/cr/rotfhrsim_spt_run3_2009_map2008_150_gauss_set1.dat'
fileout2 = '/data11/cr/rotfhrsim_spt_run3_2009_map2008_150_gauss_set2.dat'
rotate_partial_sim_map,filein,fileout,3228,2988,600,on3=0,un3=300
rotate_partial_sim_map,filein,fileout2,3228,2988,600,on3=300,un3=300

filein = '/data11/cr/hrsim_spt_run3_2009_map2008_220_gauss.dat'
fileout = '/data11/cr/rotfhrsim_spt_run3_2009_map2008_220_gauss_set1.dat'
fileout2 = '/data11/cr/rotfhrsim_spt_run3_2009_map2008_220_gauss_set2.dat'
rotate_partial_sim_map,filein,fileout,3228,2988,600,on3=0,un3=300
rotate_partial_sim_map,filein,fileout2,3228,2988,600,on3=300,un3=300

end
