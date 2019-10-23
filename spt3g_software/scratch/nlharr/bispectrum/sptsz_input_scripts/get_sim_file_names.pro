

function get_sim_file_names, field, band


sim_files = file_search('/data23/lizinvt/highel2500/sim_coadds/coadds_20130616_jittertimeordered_good/'+field+'/simcoadd_'+field+'_'+strtrim(string(band), 2)+'*dat')
sim_file = sim_files[0]

return, sim_file


;return, '/data23/lizinvt/highel2500/sim_coadds/coadds_20130616_jittertimeordered_good/'+field+'/simcoadd_'+field+'_'+strtrim(string(band), 2)+'_pairwt_jittertimeordered_20130616_good.dat'

end
